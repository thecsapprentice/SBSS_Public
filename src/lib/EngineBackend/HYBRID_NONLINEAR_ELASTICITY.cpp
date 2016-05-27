//#####################################################################
// Copyright 2010-2012, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class HYBRID_NONLINEAR_ELASTICITY
//#####################################################################
#include <PhysBAM_Tools/Data_Structures/HASHTABLE.h>
#include <PhysBAM_Tools/Matrices/SYMMETRIC_MATRIX_2X2.h>
#include <PhysBAM_Tools/Matrices/SYMMETRIC_MATRIX_3X3.h>

#include "HYBRID_NONLINEAR_ELASTICITY.h"
#include <Common/RANGE_ITERATOR.h>

struct NEOHOOKEAN_TAG;
struct COROTATED_TAG;
#include <SIMD_Optimized_Kernels/Kernel_Wrappers/Isotropic_Stress_Derivative/Isotropic_Stress_Derivative.h>
#include <SIMD_Optimized_Kernels/Kernel_Wrappers/Penalty_Measure_Gradient/Penalty_Measure_Gradient.h>
#include <SIMD_Optimized_Kernels/Kernel_Wrappers/Singular_Value_Decomposition/Singular_Value_Decomposition.h>
#include <SIMD_Optimized_Kernels/Kernel_Wrappers/Central_Gradient/Central_Gradient.h>

#ifdef LOG_DETAILED_PERFORMANCE
#define LOG_DETAILED_PERFORMANCE_NE
#endif

using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T_NONLINEAR_ELASTICITY> HYBRID_NONLINEAR_ELASTICITY<T_NONLINEAR_ELASTICITY>::
HYBRID_NONLINEAR_ELASTICITY(const T_INDEX& n_input,const T h_input,const TV origin)
    :BASE(n_input,h_input,origin),number_of_mesh_cells(0),number_of_mesh_nodes(0),mesh_initialized(false)

{}
//#####################################################################
// Function Initialize_Domain_Structures
//#####################################################################
template<class T_NONLINEAR_ELASTICITY> void HYBRID_NONLINEAR_ELASTICITY<T_NONLINEAR_ELASTICITY>::
Initialize_Domain_Structures()
{
    LOG::SCOPE scope("HYBRID_NONLINEAR_ELASTICITY::Initializing domain structures");
 
    node_is_active.Fill(false);
    node_is_dirichlet.Fill(false);
    node_is_active_mesh.Fill(false);
    node_is_dirichlet_mesh.Fill(false);

    // Step 1 : Identify active nodes originating from the grid

    for(RANGE_ITERATOR<d> cell_iterator(padded_cell_domain);cell_iterator.Valid();cell_iterator.Next()){
        const T_INDEX& cell_index=cell_iterator.Index();
        if(cell_type(cell_index)==INTERIOR_CELL_TYPE || cell_type(cell_index)==BOUNDARY_CELL_TYPE){
            PHYSBAM_ASSERT(unpadded_cell_domain.Lazy_Inside(cell_index));
            for(RANGE_ITERATOR<d> node_iterator(RANGE<T_INDEX>(cell_index,cell_index+1));node_iterator.Valid();node_iterator.Next()){
                const T_INDEX& node_index=node_iterator.Index();
                PHYSBAM_ASSERT(unpadded_node_domain.Lazy_Inside(node_index));
                node_is_active(node_index)=true;}}}

    // Step 2 : Identify active nodes originating from the mesh

    HASHTABLE<int,T_INDEX> node_hash;
    for(int cell=1;cell<=number_of_mesh_cells;cell++)
        if(cell_type_mesh(cell)==INTERIOR_CELL_TYPE || cell_type_mesh(cell)==BOUNDARY_CELL_TYPE){
            const T_INDEX& cell_index=cell_indices_mesh(cell);
            PHYSBAM_ASSERT(unpadded_cell_domain.Lazy_Inside(cell_index));
            int vertex=1;
            for(RANGE_ITERATOR<d> node_iterator(RANGE<T_INDEX>(cell_index,cell_index+1));node_iterator.Valid();node_iterator.Next(),vertex++){
                const T_INDEX& node_index=node_iterator.Index();
                PHYSBAM_ASSERT(unpadded_node_domain.Lazy_Inside(node_index));
                const int node=cells_mesh(cell)(vertex);
                PHYSBAM_ASSERT(node>=0 && node<=number_of_mesh_nodes);
                if(node){
                    PHYSBAM_ASSERT(node_hash.Get_Or_Insert(node,node_index)==node_index);
                    node_is_active_mesh(node)=true;}
                else
                    node_is_active(node_index)=true;}}

    // Step 3: Identify Dirichlet nodes originating from the grid

    for(RANGE_ITERATOR<d> cell_iterator(padded_cell_domain);cell_iterator.Valid();cell_iterator.Next()){
        const T_INDEX& cell_index=cell_iterator.Index();
        if(cell_type(cell_index)==DIRICHLET_CELL_TYPE)
            for(RANGE_ITERATOR<d> node_iterator(RANGE<T_INDEX>(cell_index,cell_index+1));node_iterator.Valid();node_iterator.Next()){
                const T_INDEX& node_index=node_iterator.Index();
                if(unpadded_node_domain.Lazy_Inside(node_index)){
                    node_is_active(node_index)=false;
                    node_is_dirichlet(node_index)=true;}}}

    // Step 4 : Identify Dirichlet nodes originating from the mesh

    for(int cell=1;cell<=number_of_mesh_cells;cell++)
        if(cell_type_mesh(cell)==DIRICHLET_CELL_TYPE){
            const T_INDEX& cell_index=cell_indices_mesh(cell);
            PHYSBAM_ASSERT(padded_cell_domain.Lazy_Inside(cell_index));
            int vertex=1;
            for(RANGE_ITERATOR<d> node_iterator(RANGE<T_INDEX>(cell_index,cell_index+1));node_iterator.Valid();node_iterator.Next(),vertex++){
                const T_INDEX& node_index=node_iterator.Index();
                PHYSBAM_ASSERT(padded_node_domain.Lazy_Inside(node_index)); // We really don't need this, we just don't expect any situation where this will fail
                const int node=cells_mesh(cell)(vertex);
                PHYSBAM_ASSERT(node>=0 && node<=number_of_mesh_nodes);
                if(node){
                    PHYSBAM_ASSERT(node_hash.Get_Or_Insert(node,node_index)==node_index);
                    node_is_active_mesh(node)=false;
                    node_is_dirichlet_mesh(node)=true;}
                else{
                    node_is_active(node_index)=false;
                    node_is_dirichlet(node_index)=true;}}}
}
//#####################################################################
// Function Initialize_Mesh
//#####################################################################
template<class T_NONLINEAR_ELASTICITY> void HYBRID_NONLINEAR_ELASTICITY<T_NONLINEAR_ELASTICITY>::
Initialize_Mesh(const int number_of_mesh_cells_input,const int number_of_mesh_nodes_input)
{
    LOG::SCOPE scope("HYBRID_NONLINEAR_ELASTICITY::Initialize_Mesh()");

    // Asserts ensure this is only called once
    PHYSBAM_ASSERT(!mesh_initialized);
    number_of_mesh_cells=number_of_mesh_cells_input;
    number_of_mesh_nodes=number_of_mesh_nodes_input;
    mesh_initialized=true;

    cell_type_mesh.Resize(number_of_mesh_cells);
    cell_indices_mesh.Resize(number_of_mesh_cells);
    cells_mesh.Resize(number_of_mesh_cells);

    node_is_active_mesh.Resize(number_of_mesh_nodes);
    node_is_dirichlet_mesh.Resize(number_of_mesh_nodes);

    U_mesh.Resize(number_of_mesh_cells);
    V_mesh.Resize(number_of_mesh_cells);
    Sigma_mesh.Resize(number_of_mesh_cells);
    Q_hat_mesh.Resize(number_of_mesh_cells);
    dP_dF_mesh.Resize(number_of_mesh_cells);    
}
//#####################################################################
// Function Initialize_Blocks
//#####################################################################
template<class T_NONLINEAR_ELASTICITY> void HYBRID_NONLINEAR_ELASTICITY<T_NONLINEAR_ELASTICITY>::
Initialize_Blocks(const int number_of_partitions)
{
    Initialize_Blocks_Elasticity(number_of_partitions);
    Initialize_Blocks_Muscles(number_of_partitions);
    Initialize_Blocks_Constraints(number_of_partitions);
}

//#####################################################################
// Function Update_Position_Based_State
//#####################################################################
template<class T_NONLINEAR_ELASTICITY> void HYBRID_NONLINEAR_ELASTICITY<T_NONLINEAR_ELASTICITY>::
Update_Position_Based_State(const T_STATE& state_u, T_STATE& state_d)
{
    LOG::SCOPE scope("HYBRID_NONLINEAR_ELASTICITY::Update_Position_Based_State()");
#ifdef USE_SPECIALIZED_KERNELS
    Update_Position_Based_State_Specialized(BASE::View_Convert(state_u.x),state_u.p,state_d.x,View_Convert(state_u.x_mesh),state_u.p_mesh,state_d.x_mesh);
    Update_Position_Based_State_Muscles_And_Constraints(state_u);
#else
    BASE::Update_Position_Based_State_Elasticity(BASE::View_Convert(state_u.x),state_u.p,state_d.x);
    Update_Position_Based_State_Elasticity_Mesh(BASE::View_Convert(state_u.x),state_u.p,state_d.x,View_Convert(state_u.x_mesh),state_u.p_mesh,state_d.x_mesh);
    Update_Position_Based_State_Muscles_And_Constraints(state_u);
#endif
}
//#####################################################################
// Function Update_Position_Based_State_Elasticity_Mesh
//#####################################################################
template<class T_NONLINEAR_ELASTICITY> void HYBRID_NONLINEAR_ELASTICITY<T_NONLINEAR_ELASTICITY>::
Update_Position_Based_State_Elasticity_Mesh(T_VECTOR_VARIABLE_VIEW_CONST u,T_SCALAR_VARIABLE_VIEW_CONST p,T_VECTOR_VARIABLE_VIEW diag, T_VECTOR_VARIABLE_MESH_VIEW_CONST u_mesh,T_SCALAR_VARIABLE_MESH_VIEW_CONST p_mesh, T_VECTOR_VARIABLE_MESH_VIEW diag_mesh)
{
    LOG::SCOPE scope("HYBRID_NONLINEAR_ELASTICITY::Update_Position_Based_State_Elasticity_Mesh()");

    for(int cell=1;cell<=number_of_mesh_cells;cell++){
        const T_INDEX& cell_index=cell_indices_mesh(cell);
        
        switch(cell_type_mesh(cell)){
            case BOUNDARY_CELL_TYPE:
                PHYSBAM_ASSERT(allow_boundary_cells);
            case INTERIOR_CELL_TYPE:
                {
                    MATRIX_MXN<T> Du;Du.Resize(d,vertices_per_cell);
                    {int vertex=1;
                    for(RANGE_ITERATOR<d> iterator(RANGE<T_INDEX>(cell_index,cell_index+1));iterator.Valid();iterator.Next(),vertex++){
                        if( cells_mesh(cell)(vertex) != 0 ){
                            // Collect from Mesh
                            for(int v=1;v<=d;v++)
                                Du(v,vertex)=u_mesh(v)(cells_mesh(cell)(vertex));
                        }
                        else{
                            const T_INDEX& index=iterator.Index();
                            for(int v=1;v<=d;v++)
                                Du(v,vertex)=u(v)(index);}
                    }}
                    MATRIX<T,d> Fe;
                    Central_Gradient(Du,Fe,1.f/BASE::h);
                    Fe+=1.f;
                    
                    MATRIX<T,d>& Ue=U_mesh(cell);
                    MATRIX<T,d>& Ve=V_mesh(cell);
                    DIAGONAL_MATRIX<T,d>& Sigmae=Sigma_mesh(cell);
                    Singular_Value_Decomposition(Fe,Ue,Sigmae,Ve);            
                    
                    Sigmae=Sigmae.Clamp_Min(BASE::cutoff_value);
                    Penalty_Measure_Gradient<typename BASE::T_CONSTITUTIVE_MODEL>(Sigmae,Q_hat_mesh(cell));
                    ROTATED_STRESS_DERIVATIVE<T,d>& dP_dFe=dP_dF_mesh(cell);
                    if(BASE::T_MATERIAL::Is_Mooney_Rivlin())
                        Isotropic_Stress_Derivative<typename BASE::T_CONSTITUTIVE_MODEL>(dP_dFe,Sigmae,p_mesh(cell),BASE::constant_mu_10,BASE::constant_mu_01,BASE::constant_kappa_max,BASE::constant_alpha,true);
                    else
                        Isotropic_Stress_Derivative<typename BASE::T_CONSTITUTIVE_MODEL>(dP_dFe,Sigmae,p_mesh(cell),BASE::constant_mu,BASE::constant_kappa_max,BASE::constant_alpha,true);
                }
                break;
            case EXTERIOR_CELL_TYPE:
            case DIRICHLET_CELL_TYPE:
                break;
            default:
                PHYSBAM_FATAL_ERROR("Improper cell type in NONLINEAR_ELASTICITY::Update_Position_Based_State");
        }
    }
}
//#####################################################################
// Function Add_Force
//#####################################################################
template<class T_NONLINEAR_ELASTICITY> void HYBRID_NONLINEAR_ELASTICITY<T_NONLINEAR_ELASTICITY>::
Add_Force(const T_STATE& state_u, T_STATE& state_f)
{
    LOG::SCOPE scope("HYBRID_NONLINEAR_ELASTICITY::Add_Force()");

    #ifdef USE_SPECIALIZED_KERNELS
    if(first_order){
        Add_Force_First_Order_Elasticity_Specialized(BASE::View_Convert(state_u.x),state_u.p,View_Convert(state_u.x_mesh),state_u.p_mesh,
                                                     state_f.x,state_f.p,state_f.x_mesh,state_f.p_mesh);
        Add_Force_Muscles_And_Constraints(state_u, state_f);
    }
    else
        PHYSBAM_NOT_IMPLEMENTED();
    #else
    if(first_order){
        BASE::Add_Force_First_Order_Elasticity(BASE::View_Convert(state_u.x),state_u.p,state_f.x,state_f.p);
        Add_Force_Muscles_And_Constraints(state_u,state_f);
        Add_Force_First_Order_Elasticity_Mesh(BASE::View_Convert(state_u.x),state_u.p,View_Convert(state_u.x_mesh),state_u.p_mesh,
                                              state_f.x,state_f.p,state_f.x_mesh,state_f.p_mesh);
    }
    else
        PHYSBAM_NOT_IMPLEMENTED();
    #endif

}
//#####################################################################
// Function Add_Force_First_Order
//#####################################################################
template<class T_NONLINEAR_ELASTICITY> void HYBRID_NONLINEAR_ELASTICITY<T_NONLINEAR_ELASTICITY>::
Add_Force_First_Order_Elasticity_Mesh(T_VECTOR_VARIABLE_VIEW_CONST u,T_SCALAR_VARIABLE_VIEW_CONST p,T_VECTOR_VARIABLE_MESH_VIEW_CONST u_mesh,T_SCALAR_VARIABLE_MESH_VIEW_CONST p_mesh, T_VECTOR_VARIABLE_VIEW f,T_SCALAR_VARIABLE_VIEW q,T_VECTOR_VARIABLE_MESH_VIEW f_mesh,T_SCALAR_VARIABLE_MESH_VIEW q_mesh) const
{
    LOG::SCOPE scope("HYBRID_NONLINEAR_ELASTICITY::Add_Force_First_Order_Elasticity_Mesh()");

    Add_Force_Stab_Mesh(u,f,u_mesh,f_mesh);
    
    const T cell_volume=pow<d>(h);
    
    for(int cell=1;cell<=number_of_mesh_cells;cell++){
        const T_INDEX& cell_index=cell_indices_mesh(cell);
        
        if(cell_type_mesh(cell)==INTERIOR_CELL_TYPE || cell_type_mesh(cell)==BOUNDARY_CELL_TYPE){
            
            DIAGONAL_MATRIX<T,d> P_hat;
            if(T_MATERIAL::Is_Mooney_Rivlin())
                P_hat=T_MATERIAL::P_hat(Sigma_mesh(cell),Q_hat_mesh(cell),p_mesh(cell),BASE::constant_mu_10,BASE::constant_mu_01,BASE::constant_alpha,BASE::constant_kappa);
            else
                P_hat=T_MATERIAL::P_hat(Sigma_mesh(cell),Q_hat_mesh(cell),p_mesh(cell),BASE::constant_mu,BASE::constant_alpha,BASE::constant_kappa);

            MATRIX<T,d> P=U_mesh(cell)*P_hat.Times_Transpose(V_mesh(cell));
            MATRIX_MXN<T> H=-P*G_One*cell_volume;
            {int vertex=1;
                for(RANGE_ITERATOR<d> iterator(RANGE<T_INDEX>(cell_index,cell_index+1));iterator.Valid();iterator.Next(),vertex++){
                    if( cells_mesh(cell)(vertex) != 0 ){
                        // Accumulate to Mesh
                        for(int v=1;v<=d;v++)
                            f_mesh(v)(cells_mesh(cell)(vertex))+=H(v,vertex);
                    }
                    else{
                        const T_INDEX& index=iterator.Index();
                        for(int v=1;v<=d;v++)
                            f(v)(index)+=H(v,vertex);
                    }
                }
            }
            
            q_mesh(cell)+=T_MATERIAL::q(Sigma_mesh(cell),p_mesh(cell),BASE::constant_alpha,BASE::constant_alpha_squared_over_kappa)*cell_volume;
        }
    }

}
//#####################################################################
// Function Add_Force_Stab_Mesh
//#####################################################################
template<class T_NONLINEAR_ELASTICITY> void HYBRID_NONLINEAR_ELASTICITY<T_NONLINEAR_ELASTICITY>::
Add_Force_Stab_Mesh(T_VECTOR_VARIABLE_VIEW_CONST dx,T_VECTOR_VARIABLE_VIEW df,T_VECTOR_VARIABLE_MESH_VIEW_CONST dx_mesh,T_VECTOR_VARIABLE_MESH_VIEW df_mesh) const
{
    for(int cell=1;cell<=number_of_mesh_cells;cell++){
        const T_INDEX& cell_index=cell_indices_mesh(cell);
        
        if(cell_type_mesh(cell)==INTERIOR_CELL_TYPE || (first_order && cell_type_mesh(cell)==BOUNDARY_CELL_TYPE)){

	    MATRIX_MXN<T> dDs;dDs.Resize(d,vertices_per_cell);

            {int vertex=1;
            for(RANGE_ITERATOR<d> iterator(RANGE<T_INDEX>(cell_index,cell_index+1));iterator.Valid();iterator.Next(),vertex++){
                if(cells_mesh(cell)(vertex)!=0){
                    // Collect from Mesh
                    for(int v=1;v<=d;v++)
                        dDs(v,vertex)=dx_mesh(v)(cells_mesh(cell)(vertex));
                }
                else{
                    const T_INDEX& index=iterator.Index();
                    for(int v=1;v<=d;v++)
                        dDs(v,vertex)=dx(v)(index);}
            }}

            MATRIX_MXN<T> dH=dDs*BASE::K_stab*(BASE::stabilization_factor*BASE::constant_mu);

            {int vertex=1;
            for(RANGE_ITERATOR<d> iterator(RANGE<T_INDEX>(cell_index,cell_index+1));iterator.Valid();iterator.Next(),vertex++){
                if(cells_mesh(cell)(vertex)!=0){
                    // Accumulate to Mesh
                    for(int v=1;v<=d;v++)
                        df_mesh(v)(cells_mesh(cell)(vertex))+=dH(v,vertex);
                }
                else{
                    const T_INDEX& index=iterator.Index();
                    for(int v=1;v<=d;v++)
                        df(v)(index)+=dH(v,vertex);}
            }}
        }
    }
}
//#####################################################################
// Function Add_Force_Differential
//#####################################################################
template<class T_NONLINEAR_ELASTICITY> void HYBRID_NONLINEAR_ELASTICITY<T_NONLINEAR_ELASTICITY>::
Add_Force_Differential(const T_STATE& state_du, T_STATE& state_df) const
{
    // LOG::SCOPE scope("HYBRID_NONLINEAR_ELASTICITY::Add_Force_Differential()");
    #ifdef USE_SPECIALIZED_KERNELS
    Add_Force_Differential_Elasticity_Specialized(BASE::View_Convert(state_du.x),state_du.p,View_Convert(state_du.x_mesh),state_du.p_mesh,
                                                  state_df.x,state_df.p,state_df.x_mesh,state_df.p_mesh);
    Add_Force_Differential_Muscles_And_Constraints(state_du,state_df);
    #else
    BASE::Add_Force_Differential_Elasticity(BASE::View_Convert(state_du.x),state_du.p,state_df.x,state_df.p);
    Add_Force_Differential_Muscles_And_Constraints(state_du,state_df);
    Add_Force_Differential_Elasticity_Mesh(BASE::View_Convert(state_du.x),state_du.p,View_Convert(state_du.x_mesh),state_du.p_mesh,
                                           state_df.x,state_df.p,state_df.x_mesh,state_df.p_mesh);
    #endif

}
//#####################################################################
// Function Add_Force_Differential
//#####################################################################
template<class T_NONLINEAR_ELASTICITY> void HYBRID_NONLINEAR_ELASTICITY<T_NONLINEAR_ELASTICITY>::
Add_Force_Differential(const T_STATE& state_du, T_STATE& state_df, const int subdomain) const
{
    // LOG::SCOPE scope("HYBRID_NONLINEAR_ELASTICITY::Add_Force_Differential()");
    #ifdef USE_SPECIALIZED_KERNELS
    Add_Force_Differential_Elasticity_Specialized(BASE::View_Convert(state_du.x),state_du.p,View_Convert(state_du.x_mesh),state_du.p_mesh,
                                                  state_df.x,state_df.p,state_df.x_mesh,state_df.p_mesh, subdomain);
    Add_Force_Differential_Muscles_And_Constraints(state_du,state_df); //TODO: Fix this!
    #else
    PHYSBAM_NOT_IMPLEMENTED();
    #endif

}
//#####################################################################
// Function Add_Force_Differential
//#####################################################################
template<class T_NONLINEAR_ELASTICITY> void HYBRID_NONLINEAR_ELASTICITY<T_NONLINEAR_ELASTICITY>::
Add_Force_Differential_Elasticity_Mesh(T_VECTOR_VARIABLE_VIEW_CONST du,T_SCALAR_VARIABLE_VIEW_CONST dp,T_VECTOR_VARIABLE_MESH_VIEW_CONST du_mesh,T_SCALAR_VARIABLE_MESH_VIEW_CONST dp_mesh,
                                       T_VECTOR_VARIABLE_VIEW df,T_SCALAR_VARIABLE_VIEW dq,T_VECTOR_VARIABLE_MESH_VIEW df_mesh,T_SCALAR_VARIABLE_MESH_VIEW dq_mesh) const
{
    Add_Force_Stab_Mesh(du,df,du_mesh,df_mesh);
    
    const T cell_volume=pow<d>(BASE::h);
    
    for(int cell=1;cell<=number_of_mesh_cells;cell++){
        const T_INDEX& cell_index=cell_indices_mesh(cell);
        
        if(cell_type_mesh(cell)==INTERIOR_CELL_TYPE || cell_type_mesh(cell)==BOUNDARY_CELL_TYPE){
            
            MATRIX_MXN<T> Du;Du.Resize(d,vertices_per_cell);
            {int vertex=1;
            for(RANGE_ITERATOR<d> iterator(RANGE<T_INDEX>(cell_index,cell_index+1));iterator.Valid();iterator.Next(),vertex++){
                if(cells_mesh(cell)(vertex)!=0){
                    // Collect from Mesh
                    for(int v=1;v<=d;v++)
                        Du(v,vertex)=du_mesh(v)(cells_mesh(cell)(vertex));
                }
                else{
                    const T_INDEX& index=iterator.Index();
                    for(int v=1;v<=d;v++)
                        Du(v,vertex)=du(v)(index);}
            }}
            MATRIX<T,d> dF;
            Central_Gradient(Du,dF,1.f/BASE::h);
                                
            MATRIX<T,d> dF_hat=U_mesh(cell).Transpose_Times(dF)*V_mesh(cell);
            dq_mesh(cell)+=BASE::T_MATERIAL::dq(Q_hat_mesh(cell),dF_hat,dp_mesh(cell),BASE::constant_alpha,BASE::constant_alpha_squared_over_kappa)*cell_volume;         
            MATRIX<T,d> dP_hat=BASE::T_MATERIAL::dP_hat(dP_dF_mesh(cell),Q_hat_mesh(cell),dF_hat,dp_mesh(cell),BASE::constant_alpha);
            MATRIX<T,d> dP=U_mesh(cell)*dP_hat.Times_Transpose(V_mesh(cell));
            

            MATRIX_MXN<T> H=-dP*G_One*cell_volume;
            {int vertex=1;
            for(RANGE_ITERATOR<d> iterator(RANGE<T_INDEX>(cell_index,cell_index+1));iterator.Valid();iterator.Next(),vertex++){
                if(cells_mesh(cell)(vertex)!=0){
                    // Accumulate to Mesh
                    for(int v=1;v<=d;v++)
                        df_mesh(v)(cells_mesh(cell)(vertex))+=H(v,vertex);
                }
                else{
                    const T_INDEX& index=iterator.Index();
                    for(int v=1;v<=d;v++)
                        df(v)(index)+=H(v,vertex);}
            }}
        }
    }
}

//#####################################################################
// Function Initialize_Undeformed_Configuration
//#####################################################################
template<class T_NONLINEAR_ELASTICITY> void HYBRID_NONLINEAR_ELASTICITY<T_NONLINEAR_ELASTICITY>::
Initialize_Undeformed_Configuration(T_STATE& state)
{
    for(RANGE_ITERATOR<d> iterator(unpadded_node_domain);iterator.Valid();iterator.Next()){
        const T_INDEX& index=iterator.Index();
        if(node_is_active(index) || node_is_dirichlet(index))
            for(int v=1;v<=d;v++)
                state.x(v)(index)=T(0);
    }
    for(int iterator = 1;iterator <= Number_Of_Mesh_Nodes();iterator++){
        if(node_is_active_mesh(iterator) || node_is_dirichlet_mesh(iterator))
            for(int v=1;v<=d;v++)
                state.x_mesh(v)(iterator)=T(0);
    }
/*
    for(RANGE_ITERATOR<d> iterator(unpadded_cell_domain);iterator.Valid();iterator.Next()){
        const T_INDEX& index=iterator.Index();
        state.p(index)=T(0);
    }
    for(int iterator = 1;iterator <= Number_Of_Mesh_Cells();iterator++){
        state.p_mesh(iterator)=T(0);
    }
*/
}

//#####################################################################
// Function Deformation
//#####################################################################
template<class T_NONLINEAR_ELASTICITY> typename HYBRID_NONLINEAR_ELASTICITY<T_NONLINEAR_ELASTICITY>::TV HYBRID_NONLINEAR_ELASTICITY<T_NONLINEAR_ELASTICITY>::
Deformation(const T_INDEX& cell_index,const TV& multilinear_coordinates, const T_STATE& state) const
{
   for(int m=1;m<=number_of_mesh_cells;m++)
       if(cell_indices_mesh(m) == cell_index)
           return Deformation_Mesh(m, multilinear_coordinates, BASE::View_Convert(state.x), View_Convert(state.x_mesh));
   return Deformation_Grid(cell_index,multilinear_coordinates, BASE::View_Convert(state.x));
}
//#####################################################################
// Function Multilinear_Interpolation_Stencil (GRID)
//#####################################################################
template<class T_NONLINEAR_ELASTICITY> typename HYBRID_NONLINEAR_ELASTICITY<T_NONLINEAR_ELASTICITY>::T_STENCIL HYBRID_NONLINEAR_ELASTICITY<T_NONLINEAR_ELASTICITY>::
Multilinear_Interpolation_Stencil_Grid(const T_INDEX& cell_index,const TV& multilinear_coordinates)
{
    T_STENCIL interpolation_stencil;
    for(RANGE_ITERATOR<d> iterator(RANGE<T_INDEX>(cell_index,cell_index+1));iterator.Valid();iterator.Next()){
        T_INDEX node_index=iterator.Index();
        T weight=(T)1.;
        for(int v=1;v<=d;v++) weight*=node_index(v)==cell_index(v)?(T)1.-multilinear_coordinates(v):multilinear_coordinates(v);
        interpolation_stencil.Insert(node_index,weight);}
    return interpolation_stencil;
}
//#####################################################################
// Function Multilinear_Interpolation_Stencil (MESH)
//#####################################################################
template<class T_NONLINEAR_ELASTICITY> typename HYBRID_NONLINEAR_ELASTICITY<T_NONLINEAR_ELASTICITY>::T_STENCIL HYBRID_NONLINEAR_ELASTICITY<T_NONLINEAR_ELASTICITY>::
Multilinear_Interpolation_Stencil_Mesh(const TV& multilinear_coordinates)
{
    return Multilinear_Interpolation_Stencil_Grid(T_INDEX(), multilinear_coordinates);
}
//#####################################################################
// Function Deformation (GRID)
//#####################################################################
template<class T_NONLINEAR_ELASTICITY> typename HYBRID_NONLINEAR_ELASTICITY<T_NONLINEAR_ELASTICITY>::TV HYBRID_NONLINEAR_ELASTICITY<T_NONLINEAR_ELASTICITY>::
Deformation_Grid(const T_INDEX& cell_index,const TV& multilinear_coordinates, T_VECTOR_VARIABLE_VIEW_CONST u) const
{
#if 0
    T_STENCIL interpolation_stencil=Multilinear_Interpolation_Stencil_Grid(T_INDEX(),multilinear_coordinates);
    TV result;
    T_VECTOR_VARIABLE_RAW temp_u;
    for(int v=1;v<=d;v++) temp_u(v).Resize(RANGE<T_INDEX>(T_INDEX(), T_INDEX()+1));
    for(RANGE_ITERATOR<d> iterator(RANGE<T_INDEX>(T_INDEX(),T_INDEX()+1));iterator.Valid();iterator.Next()){
        const T_INDEX& node_index = iterator.Index();
        for(int v=1;v<=d;v++) temp_u(v)( node_index ) = BASE::internal_state.x(v)(cell_index+node_index);
    }
    for(int v=1;v<=d;v++) result(v)=interpolation_stencil*temp_u(v);
    return result+BASE::grid.Node(cell_index)+BASE::grid.dX*multilinear_coordinates;
#else
    return Displacement_Grid(cell_index, multilinear_coordinates, u) + BASE::grid.Node(cell_index)+BASE::grid.dX*multilinear_coordinates;
#endif
}
//#####################################################################
// Function Displacement (GRID)
//#####################################################################
template<class T_NONLINEAR_ELASTICITY> typename HYBRID_NONLINEAR_ELASTICITY<T_NONLINEAR_ELASTICITY>::TV HYBRID_NONLINEAR_ELASTICITY<T_NONLINEAR_ELASTICITY>::
Displacement_Grid(const T_INDEX& cell_index,const TV& multilinear_coordinates, T_VECTOR_VARIABLE_VIEW_CONST du)
{
    T_STENCIL interpolation_stencil=Multilinear_Interpolation_Stencil_Grid(T_INDEX(),multilinear_coordinates);
    TV result;
    T_VECTOR_VARIABLE_RAW temp_du;
    for(int v=1;v<=d;v++) temp_du(v).Resize(RANGE<T_INDEX>(T_INDEX(), T_INDEX()+1));
    for(RANGE_ITERATOR<d> iterator(RANGE<T_INDEX>(T_INDEX(),T_INDEX()+1));iterator.Valid();iterator.Next()){
        const T_INDEX& node_index = iterator.Index();
        for(int v=1;v<=d;v++) temp_du(v)( node_index ) = du(v)(cell_index+node_index);
    }
    for(int v=1;v<=d;v++) result(v)=interpolation_stencil*temp_du(v);
    return result;
}
//#####################################################################
// Function Deformation (MESH)
//#####################################################################
template<class T_NONLINEAR_ELASTICITY> typename HYBRID_NONLINEAR_ELASTICITY<T_NONLINEAR_ELASTICITY>::TV HYBRID_NONLINEAR_ELASTICITY<T_NONLINEAR_ELASTICITY>::
Deformation_Mesh(const int& mesh_element,const TV& multilinear_coordinates,T_VECTOR_VARIABLE_VIEW_CONST u,T_VECTOR_VARIABLE_MESH_VIEW_CONST u_mesh) const
{
#if 0
    T_STENCIL interpolation_stencil=Multilinear_Interpolation_Stencil_Mesh(multilinear_coordinates);
    TV result;

    T_VECTOR_VARIABLE_RAW temp_u;
    for(int v=1;v<=d;v++) temp_u(v).Resize(RANGE<T_INDEX>(T_INDEX(), T_INDEX()+1));
    int i = 1;
    const T_INDEX grid_index = cell_indices_mesh(mesh_element);
    for(RANGE_ITERATOR<d> iterator(RANGE<T_INDEX>(T_INDEX(),T_INDEX()+1));iterator.Valid();iterator.Next(),i++){
        const T_INDEX& node_index = iterator.Index();
        if( cells_mesh(mesh_element)(i) == 0 )
            for(int v=1;v<=d;v++) temp_u(v)( node_index ) = BASE::internal_state.x(v)(grid_index+node_index);
        else{
            for(int v=1;v<=d;v++) temp_u(v)( node_index ) = internal_state_mesh.x_mesh(v)( cells_mesh(mesh_element)(i) );
        }
    }
    for(int v=1;v<=d;v++) result(v)=interpolation_stencil*temp_u(v);
    return result+BASE::grid.Node(cell_indices_mesh(mesh_element))+BASE::grid.dX*multilinear_coordinates;
#else
    return Displacement_Mesh(mesh_element, multilinear_coordinates, u, u_mesh) + BASE::grid.Node(cell_indices_mesh(mesh_element))+BASE::grid.dX*multilinear_coordinates;
#endif
}
//#####################################################################
// Function Displacement (MESH)
//#####################################################################
template<class T_NONLINEAR_ELASTICITY> typename HYBRID_NONLINEAR_ELASTICITY<T_NONLINEAR_ELASTICITY>::TV HYBRID_NONLINEAR_ELASTICITY<T_NONLINEAR_ELASTICITY>::
Displacement_Mesh(const int& mesh_element,const TV& multilinear_coordinates,T_VECTOR_VARIABLE_VIEW_CONST du,T_VECTOR_VARIABLE_MESH_VIEW_CONST du_mesh) const
{
    T_STENCIL interpolation_stencil=Multilinear_Interpolation_Stencil_Mesh(multilinear_coordinates);
    TV result;

    T_VECTOR_VARIABLE_RAW temp_du;
    for(int v=1;v<=d;v++) temp_du(v).Resize(RANGE<T_INDEX>(T_INDEX(), T_INDEX()+1));
    int i = 1;
    const T_INDEX grid_index = cell_indices_mesh(mesh_element);
    for(RANGE_ITERATOR<d> iterator(RANGE<T_INDEX>(T_INDEX(),T_INDEX()+1));iterator.Valid();iterator.Next(),i++){
        const T_INDEX& node_index = iterator.Index();
        if( cells_mesh(mesh_element)(i) == 0 )
            for(int v=1;v<=d;v++) temp_du(v)( node_index ) = du(v)(grid_index+node_index);
        else
            for(int v=1;v<=d;v++) temp_du(v)( node_index ) = du_mesh(v)( cells_mesh(mesh_element)(i) );
    }
    for(int v=1;v<=d;v++) result(v)=interpolation_stencil*temp_du(v);
    return result;
}
//#####################################################################
template class HYBRID_NONLINEAR_ELASTICITY<SKINNING_NONLINEAR_ELASTICITY<float,2,true,true> >;
template class HYBRID_NONLINEAR_ELASTICITY<SKINNING_NONLINEAR_ELASTICITY<float,2,true,false> >;
template class HYBRID_NONLINEAR_ELASTICITY<SKINNING_NONLINEAR_ELASTICITY<float,3,true,true> >;
template class HYBRID_NONLINEAR_ELASTICITY<SKINNING_NONLINEAR_ELASTICITY<float,3,true,false> >;
template class HYBRID_NONLINEAR_ELASTICITY<NONLINEAR_ELASTICITY<float,3> >;
template class HYBRID_NONLINEAR_ELASTICITY<NONLINEAR_ELASTICITY<float,2> >;
#ifndef USE_SPECIALIZED_KERNELS
template class HYBRID_NONLINEAR_ELASTICITY<SKINNING_NONLINEAR_ELASTICITY<double,2,true,true> >;
template class HYBRID_NONLINEAR_ELASTICITY<SKINNING_NONLINEAR_ELASTICITY<double,2,true,false> >;
template class HYBRID_NONLINEAR_ELASTICITY<SKINNING_NONLINEAR_ELASTICITY<double,3,true,true> >;
template class HYBRID_NONLINEAR_ELASTICITY<SKINNING_NONLINEAR_ELASTICITY<double,3,true,false> >;
template class HYBRID_NONLINEAR_ELASTICITY<NONLINEAR_ELASTICITY<double,3> >;
template class HYBRID_NONLINEAR_ELASTICITY<NONLINEAR_ELASTICITY<double,2> >;
#endif

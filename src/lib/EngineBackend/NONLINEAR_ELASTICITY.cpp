//#####################################################################
// Copyright 2011-2013, Nathan Mitchell, Taylor Patterson, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class NONLINEAR_ELASTICITY
//#####################################################################
#include <PhysBAM_Tools/Matrices/SYMMETRIC_MATRIX_2X2.h>
#include <PhysBAM_Tools/Matrices/SYMMETRIC_MATRIX_3X3.h>

#include "NONLINEAR_ELASTICITY.h"
#include <Common_Tools/Math_Tools/RANGE_ITERATOR.h>

#ifdef LOG_DETAILED_PERFORMANCE
#define LOG_DETAILED_PERFORMANCE_NE
#endif

// #define NO_CUT_CELLS

struct NEOHOOKEAN_TAG;
struct COROTATED_TAG;
#include <SIMD_Optimized_Kernels/Kernel_Wrappers/Isotropic_Stress_Derivative/Isotropic_Stress_Derivative.h>
#include <SIMD_Optimized_Kernels/Kernel_Wrappers/Penalty_Measure_Gradient/Penalty_Measure_Gradient.h>
#include <SIMD_Optimized_Kernels/Kernel_Wrappers/Singular_Value_Decomposition/Singular_Value_Decomposition.h>
#include <SIMD_Optimized_Kernels/Kernel_Wrappers/Central_Gradient/Central_Gradient.h>

using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
namespace{
template<int d>
void Check_Array_Pointer_Conventions()
{
    // Ensure that our array conventions (of what Get_Array_Pointer() returns) match what is in PhysBAM
    ARRAY<int> dummy_array(1);
    PHYSBAM_ASSERT(&dummy_array(1)==dummy_array.Get_Array_Pointer());
}
}

template<class T,int d> NONLINEAR_ELASTICITY<T,d>::
NONLINEAR_ELASTICITY(const T_INDEX& n_input,const T h_input,const TV origin)
    :h(h_input),n(n_input),grid(n+1,RANGE<TV>(origin,origin+TV(n)*h)),
     unpadded_cell_domain(grid.Cell_Indices()),unpadded_node_domain(grid.Node_Indices()),
     padded_cell_domain(unpadded_cell_domain.min_corner-lower_cell_padding,unpadded_cell_domain.max_corner+upper_cell_padding),
     padded_node_domain(unpadded_node_domain.min_corner-lower_node_padding,unpadded_node_domain.max_corner+upper_node_padding),
     allow_boundary_cells(false), first_order(true), constant_partitions(1)
{
    Check_Array_Pointer_Conventions<d>();

    cell_type.Resize(padded_cell_domain,true,false);
    node_is_active.Resize(padded_node_domain,true,false);
    node_is_dirichlet.Resize(padded_node_domain,true,false);

    // This has been moved to Initialize_Internal_State
    //for(int v=1;v<=d;v++) u(v).Resize(padded_node_domain,true,false);
    //p.Resize(padded_cell_domain,true,false);
    //for(int v=1;v<=d;v++) f_external(v).Resize(padded_node_domain,true,false);
    // Note: q_external should be resized, but left null here to prevent accidental use


#ifdef USE_SPECIALIZED_KERNELS
    if(d==2)
#endif
    {
        U.Resize(padded_cell_domain,true,false);
        V.Resize(padded_cell_domain,true,false);
        Sigma.Resize(padded_cell_domain,true,false);
        Q_hat.Resize(padded_cell_domain,true,false);
        dP_dF.Resize(padded_cell_domain,true,false);
    }

    // Run initialization routines that only depend on initialization-time constants 

    Build_Cell_Centered_Derivative_Operator();
    Initialize_Cell_Centered_Quadrature();
    Initialize_Stabilization_Kernel();

}
//#####################################################################
// Function Build_Cell_Centered_Derivative_Operator
//#####################################################################
template<class T,int d> void NONLINEAR_ELASTICITY<T,d>::
Build_Cell_Centered_Derivative_Operator()
{
    LOG::SCOPE scope("Building cell-centered derivative operator");

    for(int w=1;w<=d;w++)
        cell_centered_derivative_operator(w).Remove_All();

    T factor=(T)1./h;for(int v=1;v<=d-1;v++) factor*=(T).5;

    for(RANGE_ITERATOR<d> iterator(RANGE<T_INDEX>(T_INDEX(),T_INDEX::All_Ones_Vector()));iterator.Valid();iterator.Next())
        for(int w=1;w<=d;w++)
            cell_centered_derivative_operator(w).Insert(iterator.Index(),iterator.Index()(w)?factor:-factor);
}
//#####################################################################
// Function Initialize_Cell_Centered_Quadrature
//#####################################################################
template<class T,int d> void NONLINEAR_ELASTICITY<T,d>::
Initialize_Cell_Centered_Quadrature()
{
    // Construct G_One (gradient) matrix for one-point quadrature

    G_One.Resize(d,vertices_per_cell);

    int vertex=1;
    for(RANGE_ITERATOR<d> iterator(RANGE<T_INDEX>(T_INDEX(),T_INDEX::All_Ones_Vector()));iterator.Valid();iterator.Next(),vertex++){
        const T_INDEX& index=iterator.Index();

	for(int v=1;v<=d;v++)
	    if(d==2)
	        if(index(v)==0) G_One(v,vertex)=-1./(2.*h); else G_One(v,vertex)=1./(2.*h);
	    else
	        if(index(v)==0) G_One(v,vertex)=-1./(4.*h); else G_One(v,vertex)=1./(4.*h);
    }
}
//#####################################################################
// Function Initialize_Gauss_Quadrature (local)
// Function Initialize_Stabilization_Kernel
//#####################################################################
namespace{
template<int d,class T> void
Initialize_Gauss_Quadrature(ARRAY<MATRIX_MXN<T> >& gradient_matrices,ARRAY<T>& cell_volume_weights,const T h)
{
    typedef VECTOR<int,d> T_INDEX;
    enum {vertices_per_cell=POWER<2,d>::value};

    // Set up locations of quadrature points and cell volume corresponding to each quadrature point

    const int number_of_quadrature_points=POWER<2,d>::value;
    MATRIX_MXN<T> quadrature_weights(d,number_of_quadrature_points);
    cell_volume_weights.Resize(number_of_quadrature_points);

    int quadrature_point=1;
    for(RANGE_ITERATOR<d> iterator(RANGE<T_INDEX>(T_INDEX(),T_INDEX::All_Ones_Vector()));iterator.Valid();iterator.Next(),quadrature_point++){
        cell_volume_weights(quadrature_point)=pow<d>(h)/((T)number_of_quadrature_points);
        const T_INDEX& index=iterator.Index();
        VECTOR<T,d> weights;
        for(int i=1;i<=d;i++)
	    if(index(i)==0) quadrature_weights(i,quadrature_point)=(1-one_over_root_three)*.5;
	    else quadrature_weights(i,quadrature_point)=(1+one_over_root_three)*.5;}

    // Construct gradient_matrices (gradient) matrix from Gauss quadrature points

    gradient_matrices.Resize(number_of_quadrature_points);

    for(int q=1;q<=number_of_quadrature_points;q++){
        gradient_matrices(q).Resize(d,vertices_per_cell);
        VECTOR<T,d> weights;
        for(int v=1;v<=d;v++)
            weights(v)=quadrature_weights(v,q);

        int vertex=1;
	for(RANGE_ITERATOR<d> iterator(RANGE<T_INDEX>(T_INDEX(),T_INDEX::All_Ones_Vector()));iterator.Valid();iterator.Next(),vertex++){
            const T_INDEX& index=iterator.Index();

            for(int v=1;v<=d;v++){
                T& coefficient=gradient_matrices(q)(v,vertex);
                coefficient=1.;

                for(int w=1;w<=d;w++)
                    if(w==v)
                        if(index(w)==0) coefficient*=-1./h; else coefficient*=1./h;
                    else
		        if(index(w)==0) coefficient*=(1.-weights(w)); else coefficient*=weights(w);}}
    }
}
}

template<class T,int d> void NONLINEAR_ELASTICITY<T,d>::
Initialize_Stabilization_Kernel()
{
    const T cell_volume=pow<d>(h);

    ARRAY<MATRIX_MXN<T> > gradient_matrices;
    ARRAY<T> cell_volume_weights;

    Initialize_Gauss_Quadrature<d>(gradient_matrices,cell_volume_weights,h);

    K_stab=G_One.Normal_Equations_Matrix()*cell_volume;
    for(int q=1;q<=gradient_matrices.m;q++)
        K_stab-=gradient_matrices(q).Normal_Equations_Matrix()*cell_volume_weights(q);

    K_stab*=2.f;
}
//#####################################################################
// Initialize_Parameters
//#####################################################################
template<class T,int d> void NONLINEAR_ELASTICITY<T,d>::
Initialize_Parameters(const T mu,const T kappa,const T alpha,const T cutoff_value_input,const T stabilization_factor_input)
{
    constant_mu=mu;
    constant_kappa=kappa;
    constant_alpha=alpha;
    cutoff_value=cutoff_value_input;
    stabilization_factor=stabilization_factor_input;

    constant_mu_stab=mu*(h/24.f)*stabilization_factor;

    //constant_kappa_max=std::min((T)5.*mu,kappa);
    constant_kappa_max=constant_kappa;

    constant_alpha_squared_over_kappa=sqr(alpha)/kappa;
    constant_lambda=kappa;
}
template<class T,int d> void NONLINEAR_ELASTICITY<T,d>::
Initialize_Parameters(const T mu_10,const T mu_01,const T kappa,const T alpha,const T cutoff_value_input,const T stabilization_factor_input)
{
    constant_mu_10=mu_10;
    constant_mu_01=mu_01;
    constant_kappa=kappa;
    constant_alpha=alpha;
    cutoff_value=cutoff_value_input;
    stabilization_factor=stabilization_factor_input;

    constant_mu=2*(mu_10+mu_01);
    constant_mu_stab=constant_mu*(h/24.f)*stabilization_factor;

    //constant_kappa_max=std::min((T)10.*(mu_10+mu_01),kappa);
    constant_kappa_max=constant_kappa;

    constant_alpha_squared_over_kappa=sqr(alpha)/kappa;
    constant_lambda=(T)one_third*(8*mu_01-4*mu_10)+kappa;
 }
//#####################################################################
// Function Initialize_Domain_Structures
//#####################################################################
template<class T,int d> void NONLINEAR_ELASTICITY<T,d>::
Initialize_Domain_Structures()
{
    LOG::SCOPE scope("Initializing domain structures");

    node_is_active.Fill(false);
    node_is_dirichlet.Fill(false);

    for(RANGE_ITERATOR<d> iterator(padded_cell_domain);iterator.Valid();iterator.Next()){
        const T_INDEX& index=iterator.Index();
        if(cell_type(index)==INTERIOR_CELL_TYPE || cell_type(index)==BOUNDARY_CELL_TYPE){
            PHYSBAM_ASSERT(unpadded_cell_domain.Lazy_Inside(index));
            for(RANGE_ITERATOR<d> dof_iterator(RANGE<T_INDEX>(index,index+1));dof_iterator.Valid();dof_iterator.Next()){
                const T_INDEX& dof_index=dof_iterator.Index();
                PHYSBAM_ASSERT(unpadded_node_domain.Lazy_Inside(dof_index));
                node_is_active(dof_index)=true;}}}

    for(RANGE_ITERATOR<d> iterator(padded_cell_domain);iterator.Valid();iterator.Next()){
        const T_INDEX& index=iterator.Index();
        if(cell_type(index)==DIRICHLET_CELL_TYPE)
            for(RANGE_ITERATOR<d> dof_iterator(RANGE<T_INDEX>(index,index+1));dof_iterator.Valid();dof_iterator.Next()){
                const T_INDEX& dof_index=dof_iterator.Index();
                if(unpadded_node_domain.Lazy_Inside(dof_index)){
                    node_is_active(dof_index)=false;
                    node_is_dirichlet(dof_index)=true;}}}
}
//#####################################################################
// Function Initialize_Blocks
//#####################################################################
namespace{
template<class T> void
Initialize_Blocks_Helper(NONLINEAR_ELASTICITY<T,2>& nonlinear_elasticity,const int number_of_partitions)
{PHYSBAM_NOT_IMPLEMENTED();}
template<class T> void
Initialize_Blocks_Helper(NONLINEAR_ELASTICITY<T,3>& nonlinear_elasticity,const int number_of_partitions)
{nonlinear_elasticity.Initialize_Blocks_Specialized(number_of_partitions);}
}

template<class T,int d> void NONLINEAR_ELASTICITY<T,d>::
Initialize_Blocks(const int number_of_partitions)
{Initialize_Blocks_Helper(*this,number_of_partitions);}

//#####################################################################
// Function Update_Position_Based_State
//#####################################################################
template<class T,int d> void NONLINEAR_ELASTICITY<T,d>::
Update_Position_Based_State(const T_STATE& state_u, T_STATE& state_d)
{
    LOG::SCOPE scope("NONLINEAR_ELASTICITY::Update_Position_Based_State()");
    
#ifdef USE_SPECIALIZED_KERNELS
    Update_Position_Based_State_Specialized<false,false>(View_Convert(state_u.x), state_u.p, state_d.x);
#else
    Update_Position_Based_State_Elasticity(View_Convert(state_u.x), state_u.p, View_Convert(state_d.x));
#endif
}
//#####################################################################
// Function Update_Position_Based_State_Elasticity
//#####################################################################
template<class T,int d> void NONLINEAR_ELASTICITY<T,d>::
Update_Position_Based_State_Elasticity(T_VECTOR_VARIABLE_VIEW_CONST u,T_SCALAR_VARIABLE_VIEW_CONST p, T_VECTOR_VARIABLE_VIEW diag)
{
    LOG::SCOPE scope("NONLINEAR_ELASTICITY::Update_Position_Based_State_Elasticity()");
        for(RANGE_ITERATOR<d> cell_iterator(unpadded_cell_domain);cell_iterator.Valid();cell_iterator.Next()){
            const T_INDEX& cell_index=cell_iterator.Index();

            switch(cell_type(cell_index))
            {
                case BOUNDARY_CELL_TYPE:
                    PHYSBAM_ASSERT(allow_boundary_cells);
                case INTERIOR_CELL_TYPE:
                    {
                        MATRIX_MXN<T> Du;Du.Resize(d,vertices_per_cell);
                        {int vertex=1;
                            for(RANGE_ITERATOR<d> iterator(RANGE<T_INDEX>(cell_index,cell_index+1));iterator.Valid();iterator.Next(),vertex++){
                                const T_INDEX& index=iterator.Index();
                                for(int v=1;v<=d;v++)
                                    Du(v,vertex)=u(v)(index);}}
                        MATRIX<T,d> Fe;
                        Central_Gradient(Du,Fe,1.f/h);
                        Fe+=1.f;

                        MATRIX<T,d>& Ue=U(cell_index);
                        MATRIX<T,d>& Ve=V(cell_index);
                        DIAGONAL_MATRIX<T,d>& Sigmae=Sigma(cell_index);
                        Singular_Value_Decomposition(Fe,Ue,Sigmae,Ve);            

                        Sigmae=Sigmae.Clamp_Min(cutoff_value);
                        Penalty_Measure_Gradient<T_CONSTITUTIVE_MODEL>(Sigmae,Q_hat(cell_index));
                        ROTATED_STRESS_DERIVATIVE<T,d>& dP_dFe=dP_dF(cell_index);
                        if(T_MATERIAL::Is_Mooney_Rivlin())
                            Isotropic_Stress_Derivative<T_CONSTITUTIVE_MODEL>(dP_dFe,Sigmae,p(cell_index),constant_mu_10,constant_mu_01,constant_kappa_max,constant_alpha,true);
                        else
                            Isotropic_Stress_Derivative<T_CONSTITUTIVE_MODEL>(dP_dFe,Sigmae,p(cell_index),constant_mu,constant_kappa_max,constant_alpha,true);
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
template<class T,int d> void NONLINEAR_ELASTICITY<T,d>::
Add_Force(const T_STATE& state_u, T_STATE& state_f) 
{
    LOG::SCOPE scope("NONLINEAR_ELASTICITY::Add_Force");
#ifdef USE_SPECIALIZED_KERNELS  
    if(first_order)
        Add_Force_First_Order_Elasticity_Specialized<false,false>(View_Convert(state_u.x),state_u.p,state_f.x,state_f.p);
    else
        Add_Force_Second_Order_Elasticity(View_Convert(state_u.x),state_u.p,state_f.x,state_f.p);
#else
    if(first_order)
        Add_Force_First_Order_Elasticity(View_Convert(state_u.x),state_u.p,state_f.x,state_f.p);
    else
        Add_Force_Second_Order_Elasticity(View_Convert(state_u.x),state_u.p,state_f.x,state_f.p);
#endif        
}
//#####################################################################
// Function Add_Force_First_Order_Elasticity
//#####################################################################
template<class T,int d> void NONLINEAR_ELASTICITY<T,d>::
Add_Force_First_Order_Elasticity(T_VECTOR_VARIABLE_VIEW_CONST u,T_SCALAR_VARIABLE_VIEW_CONST p, T_VECTOR_VARIABLE_VIEW f,T_SCALAR_VARIABLE_VIEW q) const
{
    LOG::SCOPE scope("Adding force (1-pt quadrature) No Stab");
    
    Add_Force_Stab(u,f,true);
    
    const T cell_volume=pow<d>(h);
    
    for(RANGE_ITERATOR<d> cell_iterator(unpadded_cell_domain);cell_iterator.Valid();cell_iterator.Next()){
        const T_INDEX& cell_index=cell_iterator.Index();
        
        if(cell_type(cell_index)==INTERIOR_CELL_TYPE || cell_type(cell_index)==BOUNDARY_CELL_TYPE){
            
            DIAGONAL_MATRIX<T,d> P_hat;
            if(T_MATERIAL::Is_Mooney_Rivlin())
                P_hat=T_MATERIAL::P_hat(Sigma(cell_index),Q_hat(cell_index),p(cell_index),constant_mu_10,constant_mu_01,constant_alpha,constant_kappa);
            else
                P_hat=T_MATERIAL::P_hat(Sigma(cell_index),Q_hat(cell_index),p(cell_index),constant_mu,constant_alpha,constant_kappa);
            
            MATRIX<T,d> P=U(cell_index)*P_hat.Times_Transpose(V(cell_index));
            
            MATRIX_MXN<T> H=-P*G_One*cell_volume;
            {int vertex=1;
                for(RANGE_ITERATOR<d> iterator(RANGE<T_INDEX>(cell_index,cell_index+1));iterator.Valid();iterator.Next(),vertex++){
                    const T_INDEX& index=iterator.Index();
                    for(int v=1;v<=d;v++)
                        f(v)(index)+=H(v,vertex);}
            }
            
            q(cell_index)+=T_MATERIAL::q(Sigma(cell_index),p(cell_index),constant_alpha,constant_alpha_squared_over_kappa)*cell_volume;
        }
    }
}
//#####################################################################
// Function Gradient_Matrix (local)
// Function Add_Force_Second_Order_Elasticity
//#####################################################################
namespace{
template<class T,int d> void
Gradient_Matrix(MATRIX_MXN<T>& G,const VECTOR<T,d>& weights,const T h)
{
    typedef VECTOR<int,d> T_INDEX;
    enum {vertices_per_cell=POWER<2,d>::value};

    G.Resize(d,vertices_per_cell);
    PHYSBAM_ASSERT(weights.Min()>=0. && weights.Max()<=1.);

    int vertex=1;
    for(RANGE_ITERATOR<d> iterator(RANGE<T_INDEX>(T_INDEX(),T_INDEX::All_Ones_Vector()));iterator.Valid();iterator.Next(),vertex++){
        const T_INDEX& index=iterator.Index();
        
        for(int v=1;v<=d;v++){
            T& coefficient=G(v,vertex);
            coefficient=1.;
            
            for(int w=1;w<=d;w++)
                if(w==v)
                    if(index(w)==0) coefficient*=-1./h; else coefficient*=1./h;
                else
                    if(index(w)==0) coefficient*=(1.-weights(w)); else coefficient*=weights(w);}}
}
}

template<class T,int d> void NONLINEAR_ELASTICITY<T,d>::
Add_Force_Second_Order_Elasticity(T_VECTOR_VARIABLE_VIEW_CONST u,T_SCALAR_VARIABLE_VIEW_CONST p, T_VECTOR_VARIABLE_VIEW f,T_SCALAR_VARIABLE_VIEW q) const
{
    LOG::SCOPE scope("Adding force (1-pt quadrature) No Stab");

    Add_Force_Stab(u,f,false);

    const T cell_volume=pow<d>(h);

    for(RANGE_ITERATOR<d> cell_iterator(unpadded_cell_domain);cell_iterator.Valid();cell_iterator.Next()){
        const T_INDEX& cell_index=cell_iterator.Index();

        switch(cell_type(cell_index))
        {
            case INTERIOR_CELL_TYPE:
                {
                    MATRIX_MXN<T> Du;Du.Resize(d,vertices_per_cell);
                    {int vertex=1;
                        for(RANGE_ITERATOR<d> iterator(RANGE<T_INDEX>(cell_index,cell_index+1));iterator.Valid();iterator.Next(),vertex++){
                            const T_INDEX& index=iterator.Index();
                            for(int v=1;v<=d;v++)
                                Du(v,vertex)=u(v)(index);}}

                    MATRIX<T,d> Fe,Ue,Ve;
                    DIAGONAL_MATRIX<T,d> Sigmae,Qe_hat,Pe_hat;
                    Central_Gradient(Du,Fe,1.f/h);
                    Fe+=1.f;

                    Singular_Value_Decomposition(Fe,Ue,Sigmae,Ve);            

                    Sigmae=Sigmae.Clamp_Min(cutoff_value);
                    Penalty_Measure_Gradient<T_CONSTITUTIVE_MODEL>(Sigmae,Qe_hat);

                    if(T_MATERIAL::Is_Mooney_Rivlin())
                        Pe_hat=T_MATERIAL::P_hat(Sigmae,Qe_hat,p(cell_index),constant_mu_10,constant_mu_01,constant_alpha,constant_kappa);
                    else
                        Pe_hat=T_MATERIAL::P_hat(Sigmae,Qe_hat,p(cell_index),constant_mu,constant_alpha,constant_kappa);

                    MATRIX<T,d> Pe=Ue*Pe_hat.Times_Transpose(Ve);

                    MATRIX_MXN<T> H=-Pe*G_One*cell_volume;
                    {int vertex=1;
                        for(RANGE_ITERATOR<d> iterator(RANGE<T_INDEX>(cell_index,cell_index+1));iterator.Valid();iterator.Next(),vertex++){
                            const T_INDEX& index=iterator.Index();
                            for(int v=1;v<=d;v++)
                                f(v)(index)+=H(v,vertex);}
                    }

                    q(cell_index)+=T_MATERIAL::q(Sigmae,p(cell_index),constant_alpha,constant_alpha_squared_over_kappa)*cell_volume;
                }
                break;
            case BOUNDARY_CELL_TYPE:
                {
                    PHYSBAM_FATAL_ERROR("Boundary cell support not fully implemented");
                    if(d==2){
                        PHYSBAM_ASSERT(quadrature_point_coordinates(cell_index).m==3);
                        PHYSBAM_ASSERT(quadrature_point_weights(cell_index).m==3);}
                    else{
                        PHYSBAM_ASSERT(quadrature_point_coordinates(cell_index).m==4);
                        PHYSBAM_ASSERT(quadrature_point_weights(cell_index).m==4);}

                    for(int i=1;i<=quadrature_point_coordinates(cell_index).m;i++){
                        const TV& qp_coordinates=quadrature_point_coordinates(cell_index)(i);
                        const T qp_weight=quadrature_point_weights(cell_index)(i);

                        MATRIX_MXN<T> Du;Du.Resize(d,vertices_per_cell);
                        {int vertex=1;
                            for(RANGE_ITERATOR<d> iterator(RANGE<T_INDEX>(cell_index,cell_index+1));iterator.Valid();iterator.Next(),vertex++){
                                const T_INDEX& index=iterator.Index();
                                for(int v=1;v<=d;v++)
                                    Du(v,vertex)=u(v)(index);}}
                        MATRIX_MXN<T> G;Gradient_Matrix(G,qp_coordinates,h);
                        MATRIX<T,d> Fe=MATRIX<T,d>(Du*G.Transposed())+1.;

                        MATRIX<T,d> Ue,Ve;DIAGONAL_MATRIX<T,d> Sigmae;
                        Fe.Fast_Singular_Value_Decomposition(Ue,Sigmae,Ve);
                        Sigmae=Sigmae.Clamp_Min(cutoff_value);

                        DIAGONAL_MATRIX<T,d> P_hat;
                        if(T_MATERIAL::Is_Mooney_Rivlin())
                            P_hat=T_MATERIAL::P_hat(Sigmae,T_MATERIAL::Q_hat(Sigmae),p(cell_index),constant_mu_10,constant_mu_01,constant_alpha,constant_kappa);
                        else
                            P_hat=T_MATERIAL::P_hat(Sigmae,T_MATERIAL::Q_hat(Sigmae),p(cell_index),constant_mu,constant_alpha,constant_kappa);
                        MATRIX<T,d> P=Ue*P_hat.Times_Transpose(Ve);

                        MATRIX_MXN<T> H=-P*G*cell_volume*qp_weight;

                        {int vertex=1;
                            for(RANGE_ITERATOR<d> iterator(RANGE<T_INDEX>(cell_index,cell_index+1));iterator.Valid();iterator.Next(),vertex++){
                                const T_INDEX& index=iterator.Index();
                                for(int v=1;v<=d;v++)
                                    f(v)(index)+=H(v,vertex);}}

                        q(cell_index)+=T_MATERIAL::q(Sigmae,p(cell_index),constant_alpha,constant_alpha_squared_over_kappa)*cell_volume*qp_weight;
                    }
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
// Function Add_Force_Differential_Stab
//#####################################################################
template<class T,int d> void NONLINEAR_ELASTICITY<T,d>::
Add_Force_Stab(T_VECTOR_VARIABLE_VIEW_CONST dx,T_VECTOR_VARIABLE_VIEW df,const bool first_order) const
{
#ifdef LOG_DETAILED_PERFORMANCE_NE
    LOG::SCOPE scope("NONLINEAR_ELASTICITY::Add_Force_Differential_Stab");
#endif

    for(RANGE_ITERATOR<d> cell_iterator(unpadded_cell_domain);cell_iterator.Valid();cell_iterator.Next()){
        const T_INDEX& cell_index=cell_iterator.Index();

        if(cell_type(cell_index)==INTERIOR_CELL_TYPE || (first_order && cell_type(cell_index)==BOUNDARY_CELL_TYPE)){

	    MATRIX_MXN<T> dDs;dDs.Resize(d,vertices_per_cell);

            {int vertex=1;
            for(RANGE_ITERATOR<d> iterator(RANGE<T_INDEX>(cell_index,cell_index+1));iterator.Valid();iterator.Next(),vertex++){
                const T_INDEX& index=iterator.Index();
                for(int v=1;v<=d;v++)
		    dDs(v,vertex)=dx(v)(index);}}

            MATRIX_MXN<T> dH=dDs*K_stab*(stabilization_factor*constant_mu);

	    {int vertex=1;
	    for(RANGE_ITERATOR<d> iterator(RANGE<T_INDEX>(cell_index,cell_index+1));iterator.Valid();iterator.Next(),vertex++){
                const T_INDEX& index=iterator.Index();
                for(int v=1;v<=d;v++)
                    df(v)(index)+=dH(v,vertex);}}
        }
    }
}
//#####################################################################
// Function Add_Force_Differential
//#####################################################################
template<class T,int d> void NONLINEAR_ELASTICITY<T,d>::
Add_Force_Differential(const T_STATE& state_u, T_STATE& state_f ) const
{
#ifdef LOG_DETAILED_PERFORMANCE_NE
    LOG::SCOPE scope("NONLINEAR_ELASTICITY::Add_Force_Differential");
#endif
#ifdef USE_SPECIALIZED_KERNELS
    Add_Force_Differential_Elasticity_Specialized<false,false>(View_Convert(state_u.x),state_u.p,state_f.x,state_f.p);
#else
    Add_Force_Differential_Elasticity(View_Convert(state_u.x),state_u.p,state_f.x,state_f.p);
#endif
}
//#####################################################################
// Function Add_Force_Differential
//#####################################################################
template<class T,int d> void NONLINEAR_ELASTICITY<T,d>::
Add_Force_Differential(const T_STATE& state_u, T_STATE& state_f, const int subdomain ) const
{
#ifdef LOG_DETAILED_PERFORMANCE_NE
    LOG::SCOPE scope("NONLINEAR_ELASTICITY::Add_Force_Differential");
#endif
#ifdef USE_SPECIALIZED_KERNELS
    Add_Force_Differential_Elasticity_Specialized<false,false>(View_Convert(state_u.x),state_u.p,state_f.x,state_f.p, subdomain);
#else
    PHYSBAM_NOT_IMPLEMENTED();
#endif
}
//#####################################################################
// Function Add_Force_Differential_Nostab
//#####################################################################
template<class T,int d> void NONLINEAR_ELASTICITY<T,d>::
Add_Force_Differential_Elasticity(T_VECTOR_VARIABLE_VIEW_CONST du,T_SCALAR_VARIABLE_VIEW_CONST dp,T_VECTOR_VARIABLE_VIEW df,T_SCALAR_VARIABLE_VIEW dq) const
{
#ifdef LOG_DETAILED_PERFORMANCE_NE
    LOG::SCOPE scope("NONLINEAR_ELASTICITY::Add_Force_Differential_Elasticity");
#endif
    
    Add_Force_Stab(du,df,true);
    
    const T cell_volume=pow<d>(h);
    
    for(RANGE_ITERATOR<d> cell_iterator(unpadded_cell_domain);cell_iterator.Valid();cell_iterator.Next()){
        const T_INDEX& cell_index=cell_iterator.Index();
        
        if(cell_type(cell_index)==INTERIOR_CELL_TYPE || cell_type(cell_index)==BOUNDARY_CELL_TYPE){
            
            MATRIX_MXN<T> Du;Du.Resize(d,vertices_per_cell);
            {int vertex=1;
                for(RANGE_ITERATOR<d> iterator(RANGE<T_INDEX>(cell_index,cell_index+1));iterator.Valid();iterator.Next(),vertex++){
                    const T_INDEX& index=iterator.Index();
                    for(int v=1;v<=d;v++)
                        Du(v,vertex)=du(v)(index);}
            }
            MATRIX<T,d> dF;
            Central_Gradient(Du,dF,1.f/h);
            
            MATRIX<T,d> dF_hat=U(cell_index).Transpose_Times(dF)*V(cell_index);
            dq(cell_index)+=T_MATERIAL::dq(Q_hat(cell_index),dF_hat,dp(cell_index),constant_alpha,constant_alpha_squared_over_kappa)*cell_volume;         
            MATRIX<T,d> dP_hat=T_MATERIAL::dP_hat(dP_dF(cell_index),Q_hat(cell_index),dF_hat,dp(cell_index),constant_alpha);
            MATRIX<T,d> dP=U(cell_index)*dP_hat.Times_Transpose(V(cell_index));
            
            MATRIX_MXN<T> H=-dP*G_One*cell_volume;
            {int vertex=1;
                for(RANGE_ITERATOR<d> iterator(RANGE<T_INDEX>(cell_index,cell_index+1));iterator.Valid();iterator.Next(),vertex++){
                    const T_INDEX& index=iterator.Index();
                    for(int v=1;v<=d;v++)
                        df(v)(index)+=H(v,vertex);}
            }
        }
        else if(cell_type(cell_index)==BOUNDARY_CELL_TYPE)
            PHYSBAM_FATAL_ERROR();
    }
}
//#####################################################################
template class NONLINEAR_ELASTICITY<float,2>;
template class NONLINEAR_ELASTICITY<float,3>;
#ifndef USE_SPECIALIZED_KERNELS
template class NONLINEAR_ELASTICITY<double,2>;
template class NONLINEAR_ELASTICITY<double,3>;
#endif

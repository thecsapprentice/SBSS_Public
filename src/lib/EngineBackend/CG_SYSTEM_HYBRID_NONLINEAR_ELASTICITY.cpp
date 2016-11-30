//#####################################################################
// Copyright 2011, Taylor Patterson, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class CG_SYSTEM (Specialization for HYBRID_NONLINEAR_ELASTICITY)
//#####################################################################
#include <PhysBAM_Tools/Arrays_Computations/DOT_PRODUCT.h>
#include <PhysBAM_Tools/Arrays_Computations/MAGNITUDE.h>

#include "CG_SYSTEM.h"
#include "CG_VECTOR.h"
#include "HYBRID_NONLINEAR_ELASTICITY.h"
#include <Common_Tools/Math_Tools/RANGE_ITERATOR.h>


#include <CG_Optimized_Kernels/Dot_Product/Dot_Product_Helper.h>
#include <CG_Optimized_Kernels/Convergence_Norm/Convergence_Norm_Helper.h>
#include <CG_Optimized_Kernels/Project/Project_Helper.h>
#include <CG_Optimized_Kernels/Vector_Times/Vector_Times_Helper.h>

#ifdef LOG_DETAILED_PERFORMANCE
#define LOG_DETAILED_PERFORMANCE_CG_SYSTEM
#endif


#define CG_SYSTEM_THREADS 12
#define CG_SYSTEM_PARTITIONS 12
#define USE_FAST_CG_SYSTEM_KERNELS

//#define TEST_FOR_ZERO_WHILE_PROJECTING
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T_ELASTICITY> CG_SYSTEM<T_ELASTICITY >::
CG_SYSTEM(const T_ELASTICITY& elasticity_input)
:BASE(false,false),elasticity(elasticity_input),D(NULL),D_primal(NULL),D_dual(NULL)
{
    VECTOR<int,3> node_grid_dims = elasticity.padded_node_domain.max_corner - elasticity.padded_node_domain.min_corner +1 ;
    X_StrideNode = node_grid_dims(2)*node_grid_dims(3);
    Y_StrideNode =                   node_grid_dims(3);

    VECTOR<int,3> cell_grid_dims = elasticity.padded_cell_domain.max_corner - elasticity.padded_cell_domain.min_corner +1;
    X_StrideCell = cell_grid_dims(2)*cell_grid_dims(3);
    Y_StrideCell =                   cell_grid_dims(3);
}
//#####################################################################
// Function Multiply
//#####################################################################
template<class T_ELASTICITY> void CG_SYSTEM<T_ELASTICITY >::
Multiply(const VECTOR_BASE& v,VECTOR_BASE& result) const
{
#ifdef LOG_DETAILED_PERFORMANCE_CG_SYSTEM
    LOG::SCOPE scope("CG_SYSTEM::Multiply");
#endif

    typedef CG_VECTOR<CG_POLICY<T_ELASTICITY> > T_CG_VECTOR;
    T_CG_VECTOR& cgv_result=dynamic_cast<T_CG_VECTOR&>(result);
    const T_STATE& state_u = CG_VECTOR<CG_POLICY<T_ELASTICITY> >::State(v);
    T_STATE& state_f = CG_VECTOR<CG_POLICY<T_ELASTICITY> >::State(result);
    
    cgv_result.Clear();
    elasticity.Add_Force_Differential(state_u,state_f);
    cgv_result*=-1.;
}
//#####################################################################
// Function Inner_Product
//#####################################################################
template<class T_ELASTICITY> double CG_SYSTEM<T_ELASTICITY >::
Inner_Product(const VECTOR_BASE& v1,const VECTOR_BASE& v2) const
{
#ifdef LOG_DETAILED_PERFORMANCE_CG_SYSTEM
    LOG::SCOPE scope("CG_SYSTEM::Inner_Product");
#endif
    typedef ARRAY_VIEW<T> T_SCALAR_VARIABLE_MESH_VIEW;
    typedef VECTOR<T_SCALAR_VARIABLE_MESH_VIEW,d> T_VECTOR_VARIABLE_MESH_VIEW;
    typedef ARRAY_VIEW<const T> T_SCALAR_VARIABLE_MESH_VIEW_CONST;
    typedef VECTOR<T_SCALAR_VARIABLE_MESH_VIEW_CONST,d> T_VECTOR_VARIABLE_MESH_VIEW_CONST;

    const T_STATE& state_1 = CG_VECTOR<CG_POLICY<T_ELASTICITY> >::State(v1);
    const T_STATE& state_2 = CG_VECTOR<CG_POLICY<T_ELASTICITY> >::State(v2);
    T_VECTOR_VARIABLE_VIEW_CONST x_array1 = T_ELASTICITY::BASE::View_Convert(state_1.x);
    T_VECTOR_VARIABLE_VIEW_CONST x_array2 = T_ELASTICITY::BASE::View_Convert(state_2.x);
    T_SCALAR_VARIABLE_VIEW_CONST p_array1=state_1.p;
    T_SCALAR_VARIABLE_VIEW_CONST p_array2=state_2.p;
    T_VECTOR_VARIABLE_MESH_VIEW_CONST x_array1_mesh=T_ELASTICITY::View_Convert(state_1.x_mesh);
    T_VECTOR_VARIABLE_MESH_VIEW_CONST x_array2_mesh=T_ELASTICITY::View_Convert(state_2.x_mesh);
    T_SCALAR_VARIABLE_MESH_VIEW_CONST p_array1_mesh=state_1.p_mesh;
    T_SCALAR_VARIABLE_MESH_VIEW_CONST p_array2_mesh=state_2.p_mesh;

    double result=0.;

    int number_of_blocks = elasticity.cgblock_base_offsets.m;
    ARRAY<int> partition_offsets;
    for(int partition=0;partition<elasticity.constant_partitions;partition++)
        partition_offsets.Append((partition)*(number_of_blocks/ elasticity.constant_partitions )+min(partition,number_of_blocks% elasticity.constant_partitions ));
    
#if !defined(USE_FAST_CG_SYSTEM_KERNELS)
    const RANGE<T_INDEX>& unpadded_node_domain=elasticity.unpadded_node_domain;
    const T_FLAG& node_is_active=elasticity.node_is_active;
    for(RANGE_ITERATOR<d> node_iterator(unpadded_node_domain);node_iterator.Valid();node_iterator.Next()){
        const T_INDEX& node_index=node_iterator.Index();
        if(node_is_active(node_index))
            for(int v=1;v<=d;v++)
	        result+=x_array1(v)(node_index)*x_array2(v)(node_index);}

    const RANGE<T_INDEX>& unpadded_cell_domain=elasticity.unpadded_cell_domain;
    const T_CELL_TYPE_FIELD& cell_type=elasticity.cell_type;

    for(RANGE_ITERATOR<d> cell_iterator(unpadded_cell_domain);cell_iterator.Valid();cell_iterator.Next()){
        const T_INDEX& cell_index=cell_iterator.Index();
        if(cell_type(cell_index)==INTERIOR_CELL_TYPE || cell_type(cell_index)==BOUNDARY_CELL_TYPE)
            result+=p_array1(cell_index)*p_array2(cell_index);}

    for(int v=1;v<=d;v++) result+=ARRAYS_COMPUTATIONS::Dot_Product_Double_Precision(x_array1_mesh(v),x_array2_mesh(v));

    result+=ARRAYS_COMPUTATIONS::Dot_Product_Double_Precision(p_array1_mesh,p_array2_mesh);
#else
    const RANGE<T_INDEX>& padded_node_domain=elasticity.padded_node_domain;
    {
        for(int i=1;i<=d;i++)
            {


                result += MT_Streaming_Kernels::
                    Vector_Dot_Product_Helper<T,16>(&x_array1(i)(padded_node_domain.min_corner),
                                                    &x_array2(i)(padded_node_domain.min_corner),
                                                    elasticity.cg_offsets.Get_Array_Pointer(),
                                                    elasticity.cg_offsets.m,
                                                    elasticity.constant_partitions);


            }
    }

    const RANGE<T_INDEX>& padded_cell_domain=elasticity.padded_cell_domain;
    {

/*
        result += MT_Streaming_Kernels::
            Vector_Dot_Product_Helper<T,16>(&p_array1(padded_cell_domain.min_corner),
                                            &p_array2(padded_cell_domain.min_corner),
                                            elasticity.cg_offsets.Get_Array_Pointer(),
                                            elasticity.cg_offsets.m,
                                            elasticity.constant_partitions);
*/
    }

    for(int i=1;i<=d;i++)
        {
            result += MT_Streaming_Kernels::
                Vector_Dot_Product_Helper<T,16>(x_array1_mesh(i).Get_Array_Pointer(),
                                                x_array2_mesh(i).Get_Array_Pointer(),
                                                elasticity.cg_offsets_node_mesh.Get_Array_Pointer(),
                                                elasticity.cg_offsets_node_mesh.m,
                                                elasticity.constant_partitions);

        }
/*
    result += MT_Streaming_Kernels::
        Vector_Dot_Product_Helper<T,16>(p_array1_mesh.Get_Array_Pointer(),
                                        p_array2_mesh.Get_Array_Pointer(),
                                        elasticity.cg_offsets_cell_mesh.Get_Array_Pointer(),
                                        elasticity.cg_offsets_cell_mesh.m,
                                        elasticity.constant_partitions);
*/
#endif
    return result;
}
//#####################################################################
// Function Convergence_Norm
//#####################################################################
template<class T_ELASTICITY>  typename T_ELASTICITY::SCALAR CG_SYSTEM<T_ELASTICITY >::
Convergence_Norm(const VECTOR_BASE& v) const
{
#ifdef LOG_DETAILED_PERFORMANCE_CG_SYSTEM
    LOG::SCOPE scope("CG_SYSTEM::Convergence_Norm");
#endif
    typedef ARRAY_VIEW<T> T_SCALAR_VARIABLE_MESH_VIEW;
    typedef VECTOR<T_SCALAR_VARIABLE_MESH_VIEW,d> T_VECTOR_VARIABLE_MESH_VIEW;
    typedef ARRAY_VIEW<const T> T_SCALAR_VARIABLE_MESH_VIEW_CONST;
    typedef VECTOR<T_SCALAR_VARIABLE_MESH_VIEW_CONST,d> T_VECTOR_VARIABLE_MESH_VIEW_CONST;

    const T_STATE& state = CG_VECTOR<CG_POLICY<T_ELASTICITY> >::State(v);
    T_VECTOR_VARIABLE_VIEW_CONST x_array= T_ELASTICITY::BASE::View_Convert(state.x);
    T_SCALAR_VARIABLE_VIEW_CONST p_array= state.p;
    T_VECTOR_VARIABLE_MESH_VIEW_CONST x_array_mesh=T_ELASTICITY::View_Convert(state.x_mesh);
    T_SCALAR_VARIABLE_MESH_VIEW_CONST p_array_mesh=state.p_mesh;
    const T_FLAG& node_is_active=elasticity.node_is_active;

    T maximum=0.;

    int number_of_blocks = elasticity.cgblock_base_offsets.m;
    ARRAY<int> partition_offsets;
    for(int partition=0;partition<elasticity.constant_partitions;partition++)
        partition_offsets.Append((partition)*(number_of_blocks/ elasticity.constant_partitions )+min(partition,number_of_blocks% elasticity.constant_partitions ));   

#if !defined(USE_FAST_CG_SYSTEM_KERNELS) 
    const RANGE<T_INDEX>& unpadded_node_domain=elasticity.unpadded_node_domain;
    for(RANGE_ITERATOR<d> node_iterator(unpadded_node_domain);node_iterator.Valid();node_iterator.Next()){
        const T_INDEX& node_index=node_iterator.Index();
        if(node_is_active(node_index))
            for(int v=1;v<=d;v++)
	        maximum=std::max(maximum,(T)abs(x_array(v)(node_index)));}

    const RANGE<T_INDEX>& unpadded_cell_domain=elasticity.unpadded_cell_domain;
    const T_CELL_TYPE_FIELD& cell_type=elasticity.cell_type;

    for(RANGE_ITERATOR<d> cell_iterator(unpadded_cell_domain);cell_iterator.Valid();cell_iterator.Next()){
        const T_INDEX& cell_index=cell_iterator.Index();
        if(cell_type(cell_index)==INTERIOR_CELL_TYPE || cell_type(cell_index)==BOUNDARY_CELL_TYPE)
            maximum=std::max(maximum,(T)abs(p_array(cell_index)));}

    for(int v=1;v<=d;v++) maximum=std::max(maximum,ARRAYS_COMPUTATIONS::Maximum_Magnitude(x_array_mesh(v)));

    maximum=std::max(maximum,ARRAYS_COMPUTATIONS::Maximum_Magnitude(p_array_mesh));
#else
    {
        for(int i=1;i<=d;i++)
            {
                maximum=std::max(maximum,(T)MT_Streaming_Kernels::
                                 Vector_Convergence_Norm_Helper<T,16,false>(x_array(i).array.Get_Array_Pointer(),
                                                                            NULL,
                                                                            elasticity.cg_offsets.Get_Array_Pointer(),
                                                                            elasticity.cg_offsets.m,
                                                                            elasticity.constant_partitions));

            }
    }
    /*{
        maximum = std::max( maximum, (T)MT_Streaming_Kernels::
                            Vector_Convergence_Norm_Helper<T,16,false>(p_array.array.Get_Array_Pointer(),
                                                                 NULL,
                                                                 elasticity.cg_offsets.Get_Array_Pointer(),
                                                                 elasticity.cg_offsets.m,
                                                                 elasticity.constant_partitions));

                                                                 }*/

    for(int i=1;i<=d;i++)
        {
            maximum=std::max(maximum,(T)MT_Streaming_Kernels::
                             Vector_Convergence_Norm_Helper<T,16,false>(x_array_mesh(i).Get_Array_Pointer(),
                                                               NULL, //TODO Make sure this okay all the time
                                                               elasticity.cg_offsets_node_mesh.Get_Array_Pointer(),
                                                               elasticity.cg_offsets_node_mesh.m,
                                                               elasticity.constant_partitions));
        }
    /*{
        maximum=std::max(maximum,(T)MT_Streaming_Kernels::
                         Vector_Convergence_Norm_Helper<T,16,false>(p_array_mesh.Get_Array_Pointer(),
                                                       NULL, //TODO Make sure this is okay all the time
                                                       elasticity.cg_offsets_cell_mesh.Get_Array_Pointer(),
                                                       elasticity.cg_offsets_cell_mesh.m,
                                                       elasticity.constant_partitions));
                                                       }*/
#endif
    return maximum;
}
//#####################################################################
// Function Project
//#####################################################################
template<class T_ELASTICITY>  void CG_SYSTEM<T_ELASTICITY >::
Project(VECTOR_BASE& v) const
{
#ifdef LOG_DETAILED_PERFORMANCE_CG_SYSTEM
    LOG::SCOPE scope("CG_SYSTEM::Project");
#endif
    typedef ARRAY_VIEW<T> T_SCALAR_VARIABLE_MESH_VIEW;
    typedef VECTOR<T_SCALAR_VARIABLE_MESH_VIEW,d> T_VECTOR_VARIABLE_MESH_VIEW;
    typedef ARRAY_VIEW<const T> T_SCALAR_VARIABLE_MESH_VIEW_CONST;
    typedef VECTOR<T_SCALAR_VARIABLE_MESH_VIEW_CONST,d> T_VECTOR_VARIABLE_MESH_VIEW_CONST;
    typedef ARRAY<int> T_FLAG_MESH;

    T_STATE& state = CG_VECTOR<CG_POLICY<T_ELASTICITY> >::State(v);
    T_VECTOR_VARIABLE_VIEW x_array= state.x;
    T_VECTOR_VARIABLE_MESH_VIEW x_array_mesh=state.x_mesh;
    const T_FLAG& node_is_dirichlet=elasticity.node_is_dirichlet;
    const T_FLAG_MESH& node_is_dirichlet_mesh=elasticity.node_is_dirichlet_mesh;

    int number_of_blocks = elasticity.cgblock_base_offsets.m;
    ARRAY<int> partition_offsets;
    for(int partition=0;partition<elasticity.constant_partitions;partition++)
        partition_offsets.Append((partition)*(number_of_blocks/ elasticity.constant_partitions )+min(partition,number_of_blocks% elasticity.constant_partitions ));
 
#ifdef TEST_FOR_ZERO_WHILE_PROJECTING    
    const T_FLAG& node_is_active=elasticity.node_is_active;
#endif

#if !defined(USE_FAST_CG_SYSTEM_KERNELS) 
    const RANGE<T_INDEX>& unpadded_node_domain=elasticity.unpadded_node_domain;
    for(RANGE_ITERATOR<d> node_iterator(unpadded_node_domain);node_iterator.Valid();node_iterator.Next()){
        const T_INDEX& node_index=node_iterator.Index();
        if(node_is_dirichlet(node_index)){
            for(int v=1;v<=d;v++)
	        x_array(v)(node_index)=T();}
#ifdef TEST_FOR_ZERO_WHILE_PROJECTING    
        else if(!node_is_active(node_index))
            for(int v=1;v<=d;v++)
                PHYSBAM_ASSERT(x_array(v)(node_index)==T());
#endif
    }
#else
    {
        for(int i=1;i<=d;i++)
            {

                MT_Streaming_Kernels::
                    Vector_Project_Helper<T,16>(x_array(i).array.Get_Array_Pointer(),
                                             node_is_dirichlet.array.Get_Array_Pointer(),
                                                elasticity.cg_offsets.Get_Array_Pointer(),
                                                elasticity.cg_offsets.m,
                                                elasticity.constant_partitions);

            }

        for(int i=1;i<=d;i++)
            {

                MT_Streaming_Kernels::
                    Vector_Project_Helper<T,16>(x_array_mesh(i).Get_Array_Pointer(),
                                                node_is_dirichlet_mesh.Get_Array_Pointer(),
                                                elasticity.cg_offsets_node_mesh.Get_Array_Pointer(),
                                                elasticity.cg_offsets_node_mesh.m,
                                                elasticity.constant_partitions);
            }
    }
#ifdef TEST_FOR_ZERO_WHILE_PROJECTING   
    const RANGE<T_INDEX>& unpadded_node_domain=elasticity.unpadded_node_domain;
    for(RANGE_ITERATOR<d> node_iterator(unpadded_node_domain);node_iterator.Valid();node_iterator.Next()){
        const T_INDEX& node_index=node_iterator.Index();
        if(!node_is_active(node_index))
            for(int v=1;v<=d;v++)
                PHYSBAM_ASSERT(x_array(v)(node_index)==T());
    }
#endif   
#endif
#ifdef TEST_FOR_ZERO_WHILE_PROJECTING    
    T_SCALAR_VARIABLE_VIEW p_array= state.p;

    const RANGE<T_INDEX>& unpadded_cell_domain=elasticity.unpadded_cell_domain;
    const T_CELL_TYPE_FIELD& cell_type=elasticity.cell_type;

    for(RANGE_ITERATOR<d> cell_iterator(unpadded_cell_domain);cell_iterator.Valid();cell_iterator.Next()){
        const T_INDEX& cell_index=cell_iterator.Index();
        if(cell_type(cell_index)!=INTERIOR_CELL_TYPE && cell_type(cell_index)!=BOUNDARY_CELL_TYPE){
	    PHYSBAM_ASSERT(p_array(cell_index)==T());}}
#endif

    // TODO : Figure out what happens with dirichlet mesh stuff
}
//#####################################################################
// Function Set_Boundary_Conditions
//#####################################################################
template<class T_ELASTICITY>  void CG_SYSTEM<T_ELASTICITY >::
Set_Boundary_Conditions(VECTOR_BASE& x) const
{
    Project(x);
}
//#####################################################################
// Function Project_Nullspace
//#####################################################################
template<class T_ELASTICITY>  void CG_SYSTEM<T_ELASTICITY >::
Project_Nullspace(VECTOR_BASE& x) const {}
//#####################################################################
// Function Apply_Preconditioner
//#####################################################################
template<class T_ELASTICITY>  void CG_SYSTEM<T_ELASTICITY >::
Apply_Preconditioner(const VECTOR_BASE& y, VECTOR_BASE& x) const
{
#if 0
    //LOG::SCOPE scope("CG_SYSTEM::Apply_Preconditioner");

    if( D == NULL )
        PHYSBAM_FATAL_ERROR("Preconditioner failed. No diagonal matrix provided." );

    if( D_primal == NULL )
        PHYSBAM_FATAL_ERROR("Preconditioner failed. No diagonal primal matrix provided." );
    
    if( D_dual == NULL )
        PHYSBAM_FATAL_ERROR("Preconditioner failed. No diagonal dual matrix provided." );

    typedef CG_VECTOR<CG_POLICY<T_ELASTICITY> > T_CG_VECTOR;
    typedef typename CG_POLICY<T_ELASTICITY>::STATE_BINDER T_STATE_BINDER;
    const T_CG_VECTOR& cgv_y=dynamic_cast<const T_CG_VECTOR&>(y);
    T_CG_VECTOR& cgv_x=dynamic_cast<T_CG_VECTOR&>(x);

    const T_CG_VECTOR& cgv_DIAG=dynamic_cast<const T_CG_VECTOR&>(*D);
    const T_CG_VECTOR& cgv_DIAG_P=dynamic_cast<const T_CG_VECTOR&>(*D_primal);
    const T_CG_VECTOR& cgv_DIAG_D=dynamic_cast<const T_CG_VECTOR&>(*D_dual);

    const T_STATE& state_y = CG_VECTOR<CG_POLICY<T_ELASTICITY> >::State(y);
    T_STATE& state_x = CG_VECTOR<CG_POLICY<T_ELASTICITY> >::State(x);

    const T_STATE& DIAGONAL = CG_VECTOR<CG_POLICY<T_ELASTICITY> >::State(*D);
    const T_STATE& DIAGONAL_P = CG_VECTOR<CG_POLICY<T_ELASTICITY> >::State(*D_primal);
    const T_STATE& DIAGONAL_D = CG_VECTOR<CG_POLICY<T_ELASTICITY> >::State(*D_dual);

    T_STATE_BINDER r_b, d_b;
    T_STATE r, d;
    elasticity.Initialize_State(r_b, r);
    elasticity.Initialize_State(d_b, d);
    T_CG_VECTOR krylov_r(elasticity,r);
    T_CG_VECTOR krylov_d(elasticity,d);

#define GLOBAL_ITERATIONS 0
#define PRIMAL_ITERATIONS 2
#define DUAL_ITERATIONS 3

#if 0
#define ENABLE_NORM(X) X
#else
#define ENABLE_NORM(X)
#endif
           
    int iterations = 0;
    int global_iterations = GLOBAL_ITERATIONS;
    int primal_iterations = PRIMAL_ITERATIONS;
    int dual_iterations = DUAL_ITERATIONS;
    
    T norm;

    cgv_x.Clear();
    krylov_r.Clear();
    norm = Convergence_Norm(krylov_r);
    ENABLE_NORM(LOG::cout << "    " << norm << std::endl);

    r += state_y;

    norm = Convergence_Norm(krylov_r);

    ENABLE_NORM(LOG::cout << "    " << norm << std::endl);

    for(iterations = 0; iterations < global_iterations; iterations++){
        krylov_d.Clear();
        krylov_d += krylov_r;
        krylov_d *= cgv_DIAG;
        cgv_x += krylov_d;
        elasticity.Add_Force_Differential( d, r );
        Project( krylov_r );
        norm = Convergence_Norm(krylov_r); 
        ENABLE_NORM(LOG::cout << "  G " << norm << std::endl);
    }
    for(iterations = 0; iterations < primal_iterations; iterations++){
        krylov_d.Clear();
        krylov_d += krylov_r;
        krylov_d *= cgv_DIAG_P;
        cgv_x += krylov_d;
        elasticity.Add_Force_Differential( d, r );
        Project( krylov_r );
        norm = Convergence_Norm(krylov_r); 
        ENABLE_NORM(LOG::cout << "  P " << norm << std::endl);
    }
    for(iterations = 0; iterations < dual_iterations; iterations++){
        krylov_d.Clear();
        krylov_d += krylov_r;
        krylov_d *= cgv_DIAG_D;
        cgv_x += krylov_d;
        elasticity.Add_Force_Differential( d, r );
        Project( krylov_r );
        norm = Convergence_Norm(krylov_r); 
        ENABLE_NORM(LOG::cout << "  D " << norm << std::endl);
    }
    for(iterations = 0; iterations < primal_iterations; iterations++){
        krylov_d.Clear();
        krylov_d += krylov_r;
        krylov_d *= cgv_DIAG_P;
        cgv_x += krylov_d;
        elasticity.Add_Force_Differential( d, r );
        Project( krylov_r );
        norm = Convergence_Norm(krylov_r); 
        ENABLE_NORM(LOG::cout << "  P " <<  norm << std::endl);
    }
    for(iterations = 0; iterations < global_iterations; iterations++){
        krylov_d.Clear();
        krylov_d += krylov_r;
        krylov_d *= cgv_DIAG;
        cgv_x += krylov_d;
        elasticity.Add_Force_Differential( d, r );
        Project( krylov_r );
        norm = Convergence_Norm(krylov_r); 
        ENABLE_NORM(LOG::cout << "  G " << norm << std::endl);
    }
#endif
}
//#####################################################################
//template class CG_SYSTEM<float,2>;
template class CG_SYSTEM<HYBRID_NONLINEAR_ELASTICITY<SKINNING_NONLINEAR_ELASTICITY<float,3,true,true> > >;
template class CG_SYSTEM<HYBRID_NONLINEAR_ELASTICITY<SKINNING_NONLINEAR_ELASTICITY<float,3,true,false> > >;
template class CG_SYSTEM<HYBRID_NONLINEAR_ELASTICITY<NONLINEAR_ELASTICITY<float,3> > >;
#ifndef USE_SPECIALIZED_KERNELS
template class CG_SYSTEM<HYBRID_NONLINEAR_ELASTICITY<SKINNING_NONLINEAR_ELASTICITY<double,3,true,true> > >;
template class CG_SYSTEM<HYBRID_NONLINEAR_ELASTICITY<SKINNING_NONLINEAR_ELASTICITY<double,3,true,false> > >;
template class CG_SYSTEM<HYBRID_NONLINEAR_ELASTICITY<NONLINEAR_ELASTICITY<double,3> > >;
#endif



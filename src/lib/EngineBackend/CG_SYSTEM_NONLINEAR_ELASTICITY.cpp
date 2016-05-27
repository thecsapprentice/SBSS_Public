//#####################################################################
// Copyright 2011, Taylor Patterson, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class CG_SYSTEM
//#####################################################################
#include "CG_SYSTEM.h"
#include "CG_VECTOR.h"
#include "SKINNING_NONLINEAR_ELASTICITY.h"
#include "NONLINEAR_ELASTICITY.h"
#include <Common/RANGE_ITERATOR.h>

#include <CG_Optimized_Kernels/Dot_Product/Dot_Product_Helper.h>
#include <CG_Optimized_Kernels/Convergence_Norm/Convergence_Norm_Helper.h>
#include <CG_Optimized_Kernels/Project/Project_Helper.h>
#include <CG_Optimized_Kernels/Vector_Times/Vector_Times_Helper.h>

#ifdef LOG_DETAILED_PERFORMANCE
#define LOG_DETAILED_PERFORMANCE_CG_SYSTEM
#endif

#define CG_SYSTEM_THREADS 6
//#define USE_FAST_CG_SYSTEM_KERNELS

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

//    T_VECTOR_VARIABLE_VIEW_CONST u=CG_VECTOR<T,d,HYBRID_NONLINEAR_ELASTICITY<SKINNING_NONLINEAR_ELASTICITY<T,d,true,false> > >::X_Array(v);
//    T_SCALAR_VARIABLE_VIEW_CONST p=CG_VECTOR<T,d,HYBRID_NONLINEAR_ELASTICITY<SKINNING_NONLINEAR_ELASTICITY<T,d,true,false> > >::P_Array(v);
//    T_VECTOR_VARIABLE_VIEW f=CG_VECTOR<T,d,HYBRID_NONLINEAR_ELASTICITY<SKINNING_NONLINEAR_ELASTICITY<T,d,true,false> > >::X_Array(result);
//    T_SCALAR_VARIABLE_VIEW q=CG_VECTOR<T,d,HYBRID_NONLINEAR_ELASTICITY<SKINNING_NONLINEAR_ELASTICITY<T,d,true,false> > >::P_Array(result);

    const T_STATE& state_u = CG_VECTOR<CG_POLICY<T_ELASTICITY> >::State(v);
    T_STATE& state_f = CG_VECTOR<CG_POLICY<T_ELASTICITY> >::State(result);



#if !defined(USE_FAST_CG_SYSTEM_KERNELS)
    for(int v=1;v<=d;v++) state_f.x(v).Fill(T());
    state_f.p.Fill(T());
#else
    {
        for(int i=1;i<=d;i++)
            {
                MT_Streaming_Kernels::
                    Vector_Times_Helper<T>(state_f.x(i).array.Get_Array_Pointer(),
                                           T(0),
                                           elasticity.node_block_base_offsets.Get_Array_Pointer(),
                                           elasticity.node_block_base_offsets.m,
                                           elasticity.node_block_partition_offsets.Get_Array_Pointer(),
                                           elasticity.node_block_partition_offsets.m,
                                           X_StrideNode, Y_StrideNode, 3);
            }

        MT_Streaming_Kernels::
            Vector_Times_Helper<T>(state_f.p.array.Get_Array_Pointer(),
                                   T(0),
                                   elasticity.cell_block_base_offsets.Get_Array_Pointer(),
                                   elasticity.cell_block_base_offsets.m,
                                   elasticity.cell_block_partition_offsets.Get_Array_Pointer(),
                                   elasticity.cell_block_partition_offsets.m,
                                   X_StrideCell, Y_StrideCell, 2);
    }
#endif

    elasticity.Add_Force_Differential(state_u,state_f);

#if !defined(USE_FAST_CG_SYSTEM_KERNELS)
    state_f.p*=-1.;
    for(int v=1;v<=d;v++) state_f.x(v)*=-1.;
#else
    {
        for(int i=1;i<=d;i++)
            {
                MT_Streaming_Kernels::
                    Vector_Times_Helper<T>(state_f.x(i).array.Get_Array_Pointer(),
                                           T(-1),
                                           elasticity.node_block_base_offsets.Get_Array_Pointer(),
                                           elasticity.node_block_base_offsets.m,
                                           elasticity.node_block_partition_offsets.Get_Array_Pointer(),
                                           elasticity.node_block_partition_offsets.m,
                                           X_StrideNode, Y_StrideNode, 3);
            }

        MT_Streaming_Kernels::
            Vector_Times_Helper<T>(state_f.p.array.Get_Array_Pointer(),
                                   T(-1),
                                   elasticity.cell_block_base_offsets.Get_Array_Pointer(),
                                   elasticity.cell_block_base_offsets.m,
                                   elasticity.cell_block_partition_offsets.Get_Array_Pointer(),
                                   elasticity.cell_block_partition_offsets.m,
                                   X_StrideCell, Y_StrideCell, 2);
    }
#endif




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

    //T_VECTOR_VARIABLE_VIEW_CONST x_array1=CG_VECTOR<T,d,HYBRID_NONLINEAR_ELASTICITY<SKINNING_NONLINEAR_ELASTICITY<T,d,true,false> > >::X_Array(v1);
    //T_VECTOR_VARIABLE_VIEW_CONST x_array2=CG_VECTOR<T,d,HYBRID_NONLINEAR_ELASTICITY<SKINNING_NONLINEAR_ELASTICITY<T,d,true,false> > >::X_Array(v2);
    const T_STATE& state_1 = CG_VECTOR<CG_POLICY<T_ELASTICITY> >::State(v1);
    const T_STATE& state_2 = CG_VECTOR<CG_POLICY<T_ELASTICITY> >::State(v2);

    T_VECTOR_VARIABLE_VIEW_CONST x_array1 = T_ELASTICITY::View_Convert(state_1.x);
    T_VECTOR_VARIABLE_VIEW_CONST x_array2 = T_ELASTICITY::View_Convert(state_2.x);

    double result=0.;

#if !defined(USE_FAST_CG_SYSTEM_KERNELS)
    const RANGE<T_INDEX>& unpadded_node_domain=elasticity.unpadded_node_domain;
    const T_FLAG& node_is_active=elasticity.node_is_active;
    for(RANGE_ITERATOR<d> node_iterator(unpadded_node_domain);node_iterator.Valid();node_iterator.Next()){
        const T_INDEX& node_index=node_iterator.Index();
        if(node_is_active(node_index))
            for(int v=1;v<=d;v++)
	        result+=x_array1(v)(node_index)*x_array2(v)(node_index);}
#else
    const RANGE<T_INDEX>& padded_node_domain=elasticity.padded_node_domain;

/*
    double result_d=0.;
    const RANGE<T_INDEX>& unpadded_node_domain=elasticity.unpadded_node_domain;
    const T_FLAG& node_is_active=elasticity.node_is_active;
    for(int v=1;v<=d;v++){
        for(RANGE_ITERATOR<d> node_iterator(padded_node_domain);node_iterator.Valid();node_iterator.Next()){
            const T_INDEX& node_index=node_iterator.Index();
            if(node_is_active(node_index)){
                PHYSBAM_ASSERT(unpadded_node_domain.Lazy_Inside(node_index));
                result_d+=x_array1(v)(node_index)*x_array2(v)(node_index);
            }
            else{
                for(int v=1;v<=d;v++){
                    PHYSBAM_ASSERT(x_array1(v)(node_index)==0);
                    PHYSBAM_ASSERT(x_array2(v)(node_index)==0);
                }
            }
        }
    }
*/  
 
    {
/*
        int size = (padded_node_domain.max_corner-padded_node_domain.min_corner+1).Product();
        for(int i=1;i<=d;i++)
            {
                MT_Streaming_Kernels::Dot_Product_Helper<T> test(&x_array1(i)(padded_node_domain.min_corner),
                                                                 &x_array2(i)(padded_node_domain.min_corner), size);
                result += test.Run_Parallel( CG_SYSTEM_THREADS );        
            }
*/
        for(int i=1;i<=d;i++)
            {
                result += MT_Streaming_Kernels::
                    Vector_Dot_Product_Helper<T>(&x_array1(i)(padded_node_domain.min_corner),
                                                 &x_array2(i)(padded_node_domain.min_corner),
                                                 elasticity.node_block_base_offsets.Get_Array_Pointer(),
                                                 elasticity.node_block_base_offsets.m,
                                                 elasticity.node_block_partition_offsets.Get_Array_Pointer(),
                                                 elasticity.node_block_partition_offsets.m,
                                                 X_StrideNode, Y_StrideNode, 3);
            }


    }

    //LOG::cout << "Dot Product:      DEBUG:"<<result_d<<"   VECTOR:"<<result<<std::endl;

#endif


    //T_SCALAR_VARIABLE_VIEW_CONST p_array1=CG_VECTOR<T,d,HYBRID_NONLINEAR_ELASTICITY<SKINNING_NONLINEAR_ELASTICITY<T,d,true,false> > >::P_Array(v1);
    //T_SCALAR_VARIABLE_VIEW_CONST p_array2=CG_VECTOR<T,d,HYBRID_NONLINEAR_ELASTICITY<SKINNING_NONLINEAR_ELASTICITY<T,d,true,false> > >::P_Array(v2);

    T_SCALAR_VARIABLE_VIEW_CONST p_array1=state_1.p;
    T_SCALAR_VARIABLE_VIEW_CONST p_array2=state_2.p;

#if !defined(USE_FAST_CG_SYSTEM_KERNELS) 
    const RANGE<T_INDEX>& unpadded_cell_domain=elasticity.unpadded_cell_domain;
    const T_CELL_TYPE_FIELD& cell_type=elasticity.cell_type;

    for(RANGE_ITERATOR<d> cell_iterator(unpadded_cell_domain);cell_iterator.Valid();cell_iterator.Next()){
        const T_INDEX& cell_index=cell_iterator.Index();
        if(cell_type(cell_index)==INTERIOR_CELL_TYPE || cell_type(cell_index)==BOUNDARY_CELL_TYPE)
            result+=p_array1(cell_index)*p_array2(cell_index);}
#else
    const RANGE<T_INDEX>& padded_cell_domain=elasticity.padded_cell_domain;
    // double result_d=result;
    // const RANGE<T_INDEX>& unpadded_cell_domain=elasticity.unpadded_cell_domain;
    // const T_CELL_TYPE_FIELD& cell_type=elasticity.cell_type;
    // for(RANGE_ITERATOR<d> cell_iterator(padded_cell_domain);cell_iterator.Valid();cell_iterator.Next()){
    //     const T_INDEX& cell_index=cell_iterator.Index();
    //     if(cell_type(cell_index)==INTERIOR_CELL_TYPE|| cell_type(cell_index)==BOUNDARY_CELL_TYPE)
    //         {
    //             PHYSBAM_ASSERT(unpadded_cell_domain.Lazy_Inside(cell_index));
    //             result_d+=p_array1(cell_index)*p_array2(cell_index);
    //         }
    //     else
    //         {
    //             try
    //                 {
    //                     PHYSBAM_ASSERT(p_array1(cell_index)<1e-4);
    //                     PHYSBAM_ASSERT(p_array2(cell_index)<1e-4);                
    //                 }
    //             catch( ASSERTION_ERROR e)
    //                 {
    //                     LOG::cout << "P1 at " << cell_index << " is " << p_array1(cell_index) << std::endl;
    //                     LOG::cout << "P2 at " << cell_index << " is " << p_array2(cell_index) << std::endl;
    //                 }
    //         }
    // }
    
    {
        /*
        int size = (padded_cell_domain.max_corner-padded_cell_domain.min_corner+1).Product();

        MT_Streaming_Kernels::Dot_Product_Helper<T> test1(&p_array1(padded_cell_domain.min_corner),
                                                          &p_array2(padded_cell_domain.min_corner), size);

        result+=test1.Run_Parallel( CG_SYSTEM_THREADS );
        */

        result += MT_Streaming_Kernels::
            Vector_Dot_Product_Helper<T>(&p_array1(padded_cell_domain.min_corner),
                                         &p_array2(padded_cell_domain.min_corner),
                                         elasticity.cell_block_base_offsets.Get_Array_Pointer(),
                                         elasticity.cell_block_base_offsets.m,
                                         elasticity.cell_block_partition_offsets.Get_Array_Pointer(),
                                         elasticity.cell_block_partition_offsets.m,
                                         X_StrideCell, Y_StrideCell, 2);

    }


    //LOG::cout << "Dot Product:      DEBUG:"<<result_d<<"   VECTOR:"<<result<<std::endl;


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

    //T_VECTOR_VARIABLE_VIEW_CONST x_array=CG_VECTOR<T,d,HYBRID_NONLINEAR_ELASTICITY<SKINNING_NONLINEAR_ELASTICITY<T,d,true,false> > >::X_Array(v);

    const T_STATE& state = CG_VECTOR<CG_POLICY<T_ELASTICITY> >::State(v);
    T_VECTOR_VARIABLE_VIEW_CONST x_array= T_ELASTICITY::View_Convert(state.x);

    const T_FLAG& node_is_active=elasticity.node_is_active;


    T maximum=0.;

#if !defined(USE_FAST_CG_SYSTEM_KERNELS) 
    const RANGE<T_INDEX>& unpadded_node_domain=elasticity.unpadded_node_domain;
    for(RANGE_ITERATOR<d> node_iterator(unpadded_node_domain);node_iterator.Valid();node_iterator.Next()){
        const T_INDEX& node_index=node_iterator.Index();
        if(node_is_active(node_index))
            for(int v=1;v<=d;v++)
	        maximum=std::max(maximum,(T)abs(x_array(v)(node_index)));}
#else
    //const RANGE<T_INDEX>& padded_node_domain=elasticity.padded_node_domain;
    {
/*
        int size = (padded_node_domain.max_corner-padded_node_domain.min_corner+1).Product();
        
        for(int i=1;i<=d;i++)
            {
                MT_Streaming_Kernels::Convergence_Norm_Helper<T> test(x_array(i).array.Get_Array_Pointer(),
                                                                      node_is_active.array.Get_Array_Pointer(), size);
                maximum = std::max( maximum, test.Run_Parallel( CG_SYSTEM_THREADS ));
            }
*/

        for(int i=1;i<=d;i++)
            {
                maximum=std::max(maximum,(T)MT_Streaming_Kernels::
                                 Vector_Convergence_Norm_Helper<T>(x_array(i).array.Get_Array_Pointer(),
                                                                   node_is_active.array.Get_Array_Pointer(),
                                                                   elasticity.node_block_base_offsets.Get_Array_Pointer(),
                                                                   elasticity.node_block_base_offsets.m,
                                                                   elasticity.node_block_partition_offsets.Get_Array_Pointer(),
                                                                   elasticity.node_block_partition_offsets.m,
                                                                   X_StrideNode, Y_StrideNode, 3));
            }

    }
#endif

    //T_SCALAR_VARIABLE_VIEW_CONST p_array=CG_VECTOR<T,d,HYBRID_NONLINEAR_ELASTICITY<SKINNING_NONLINEAR_ELASTICITY<T,d,true,false> > >::P_Array(v);
    T_SCALAR_VARIABLE_VIEW_CONST p_array= state.p;


#if !defined(USE_FAST_CG_SYSTEM_KERNELS) 
    const RANGE<T_INDEX>& unpadded_cell_domain=elasticity.unpadded_cell_domain;
    const T_CELL_TYPE_FIELD& cell_type=elasticity.cell_type;

    for(RANGE_ITERATOR<d> cell_iterator(unpadded_cell_domain);cell_iterator.Valid();cell_iterator.Next()){
        const T_INDEX& cell_index=cell_iterator.Index();
        if(cell_type(cell_index)==INTERIOR_CELL_TYPE || cell_type(cell_index)==BOUNDARY_CELL_TYPE)
            maximum=std::max(maximum,(T)abs(p_array(cell_index)));}
#else
    //const RANGE<T_INDEX>& padded_cell_domain=elasticity.padded_cell_domain;
    
    {
/*
        int size = (padded_cell_domain.max_corner-padded_cell_domain.min_corner+1).Product();
        
        MT_Streaming_Kernels::Convergence_Norm_Helper<T> test(p_array.array.Get_Array_Pointer(), size);
        maximum = std::max( maximum, test.Run_Parallel( CG_SYSTEM_THREADS ));
*/

        maximum = std::max( maximum, (T)MT_Streaming_Kernels::
                            Vector_Convergence_Norm_Helper<T>(p_array.array.Get_Array_Pointer(),
                                                              NULL,
                                                              elasticity.cell_block_base_offsets.Get_Array_Pointer(),
                                                              elasticity.cell_block_base_offsets.m,
                                                              elasticity.cell_block_partition_offsets.Get_Array_Pointer(),
                                                              elasticity.cell_block_partition_offsets.m,
                                                              X_StrideCell, Y_StrideCell, 2));

      
    }


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

    //T_VECTOR_VARIABLE_VIEW x_array=CG_VECTOR<T,d,HYBRID_NONLINEAR_ELASTICITY<SKINNING_NONLINEAR_ELASTICITY<T,d,true,false> > >::X_Array(v);
    T_STATE& state = CG_VECTOR<CG_POLICY<T_ELASTICITY> >::State(v);
    T_VECTOR_VARIABLE_VIEW x_array= state.x;


    const T_FLAG& node_is_dirichlet=elasticity.node_is_dirichlet;

#ifdef TEST_FOR_ZERO_WHILE_PROJECTING    
    const T_FLAG& node_is_active=elasticity.node_is_active;
#endif

#if  !defined(USE_FAST_CG_SYSTEM_KERNELS) 
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
    //gconst RANGE<T_INDEX>& padded_node_domain=elasticity.padded_node_domain;
    
    {
/*
        int size = (padded_node_domain.max_corner-padded_node_domain.min_corner+1).Product();
        
        for(int i=1;i<=d;i++)
            {
                MT_Streaming_Kernels::Project_Helper<T> test(x_array(i).array.Get_Array_Pointer(),
                                                             node_is_dirichlet.array.Get_Array_Pointer(), size);
                test.Run_Parallel( CG_SYSTEM_THREADS );
            }
*/
        for(int i=1;i<=d;i++)
            {
                MT_Streaming_Kernels::
                    Vector_Project_Helper<T>(x_array(i).array.Get_Array_Pointer(),
                                             node_is_dirichlet.array.Get_Array_Pointer(),
                                             elasticity.node_block_base_offsets.Get_Array_Pointer(),
                                             elasticity.node_block_base_offsets.m,
                                             elasticity.node_block_partition_offsets.Get_Array_Pointer(),
                                             elasticity.node_block_partition_offsets.m,
                                             X_StrideNode, Y_StrideNode, 3);
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
    //T_SCALAR_VARIABLE_VIEW p_array=CG_VECTOR<T,d,HYBRID_NONLINEAR_ELASTICITY<SKINNING_NONLINEAR_ELASTICITY<T,d,true,false> > >::P_Array(v);
    T_SCALAR_VARIABLE_VIEW p_array= state.p;

    const RANGE<T_INDEX>& unpadded_cell_domain=elasticity.unpadded_cell_domain;
    const T_CELL_TYPE_FIELD& cell_type=elasticity.cell_type;

    for(RANGE_ITERATOR<d> cell_iterator(unpadded_cell_domain);cell_iterator.Valid();cell_iterator.Next()){
        const T_INDEX& cell_index=cell_iterator.Index();
        if(cell_type(cell_index)!=INTERIOR_CELL_TYPE && cell_type(cell_index)!=BOUNDARY_CELL_TYPE){
	    PHYSBAM_ASSERT(p_array(cell_index)==T());}}
#endif
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
Apply_Preconditioner(const VECTOR_BASE& r, VECTOR_BASE& z) const
{
    const T_STATE& state_r = CG_VECTOR<CG_POLICY<T_ELASTICITY> >::State(r);
    T_STATE& state_z = CG_VECTOR<CG_POLICY<T_ELASTICITY> >::State(z);

    for(int v=1;v<=d;v++) state_z.x(v).Fill(T());
    state_z.p.Fill(T());

    for(int v=1;v<=d;v++) state_z.x += state_r.x;
    state_z.p += state_r.p;
}
//#####################################################################
//template class CG_SYSTEM<float,2>;
template class CG_SYSTEM<SKINNING_NONLINEAR_ELASTICITY<float,3,true,true> >;
template class CG_SYSTEM<SKINNING_NONLINEAR_ELASTICITY<float,3,true,false> >;
template class CG_SYSTEM<NONLINEAR_ELASTICITY<float,3> >;
#ifndef USE_SPECIALIZED_KERNELS
template class CG_SYSTEM<SKINNING_NONLINEAR_ELASTICITY<double,3,true,true> >;
template class CG_SYSTEM<SKINNING_NONLINEAR_ELASTICITY<double,3,true,false> >;
template class CG_SYSTEM<NONLINEAR_ELASTICITY<double,3> >;
#endif



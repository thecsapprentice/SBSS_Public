//#####################################################################
// Copyright 2011, Taylor Patterson, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class CG_VECTOR
//#####################################################################
#include "CG_VECTOR.h"
#include "CG_POLICY.h"
#include "NONLINEAR_ELASTICITY.h"

#include <Common/RANGE_ITERATOR.h>
using namespace PhysBAM;

#ifdef LOG_DETAILED_PERFORMANCE
#define LOG_DETAILED_PERFORMANCE_CG_VECTOR
#endif


#define CG_VECTOR_THREADS 6
//#define USE_FAST_CG_VECTOR_KERNELS

#include <CG_Optimized_Kernels/Vector_Add/Vector_Add_Helper.h>
#include <CG_Optimized_Kernels/Vector_Subtract/Vector_Subtract_Helper.h>
#include <CG_Optimized_Kernels/Vector_Times/Vector_Times_Helper.h>
#include <CG_Optimized_Kernels/Vector_Set/Vector_Set_Helper.h>
#include <CG_Optimized_Kernels/Vector_SAXPY_A/Vector_SAXPY_A_Helper.h>
#include <CG_Optimized_Kernels/Vector_SAXPY_AB/Vector_SAXPY_AB_Helper.h>
#include <CG_Optimized_Kernels/Vector_Element_Times/Vector_Element_Times_Helper.h>

//#####################################################################
// operator+=
//#####################################################################
template<class POLICY> KRYLOV_VECTOR_BASE<typename POLICY::ELASTICITY::SCALAR >& CG_VECTOR<POLICY>::
operator+=(const BASE& bv)
{
#ifdef LOG_DETAILED_PERFORMANCE_CG_VECTOR
    LOG::SCOPE scope("CG_VECTOR::Addition");
#endif

#if !defined(USE_FAST_CG_VECTOR_KERNELS)
    state.x += State(bv).x;
    state.p += State(bv).p;
#else   
    {
        for(int i=1;i<=d;i++)
            {
                MT_Streaming_Kernels::
                    Vector_Add_Helper<T>(state.x(i).array.Get_Array_Pointer(),
                                         State(bv).x(i).array.Get_Array_Pointer(),
                                         ne.node_block_base_offsets.Get_Array_Pointer(),
                                         ne.node_block_base_offsets.m,
                                         ne.node_block_partition_offsets.Get_Array_Pointer(),
                                         ne.node_block_partition_offsets.m,
                                         X_StrideNode, Y_StrideNode, 3);
            }

        MT_Streaming_Kernels::
            Vector_Add_Helper<T>(state.p.array.Get_Array_Pointer(),
                                 State(bv).p.array.Get_Array_Pointer(),
                                 ne.cell_block_base_offsets.Get_Array_Pointer(),
                                 ne.cell_block_base_offsets.m,
                                 ne.cell_block_partition_offsets.Get_Array_Pointer(),
                                 ne.cell_block_partition_offsets.m,
                                 X_StrideCell, Y_StrideCell, 2);
    }
#endif
    return *this;
}
//#####################################################################
// operator-=
//#####################################################################
template<class POLICY> KRYLOV_VECTOR_BASE<typename POLICY::ELASTICITY::SCALAR >& CG_VECTOR<POLICY>::
operator-=(const BASE& bv)
{
#ifdef LOG_DETAILED_PERFORMANCE_CG_VECTOR
    LOG::SCOPE scope("CG_VECTOR::Subtraction");
#endif

#if !defined(USE_FAST_CG_VECTOR_KERNELS)
    state.x -= State(bv).x;
    state.p -= State(bv).p;
#else   
    {       
        for(int i=1;i<=d;i++)
            {
                MT_Streaming_Kernels::
                    Vector_Subtract_Helper<T>(state.x(i).array.Get_Array_Pointer(),
                                              State(bv).x(i).array.Get_Array_Pointer(),
                                              ne.node_block_base_offsets.Get_Array_Pointer(),
                                              ne.node_block_base_offsets.m,
                                              ne.node_block_partition_offsets.Get_Array_Pointer(),
                                              ne.node_block_partition_offsets.m,
                                              X_StrideNode, Y_StrideNode, 3);
            }

        MT_Streaming_Kernels::
            Vector_Subtract_Helper<T>(state.p.array.Get_Array_Pointer(),
                                      State(bv).p.array.Get_Array_Pointer(),
                                      ne.cell_block_base_offsets.Get_Array_Pointer(),
                                      ne.cell_block_base_offsets.m,
                                      ne.cell_block_partition_offsets.Get_Array_Pointer(),
                                      ne.cell_block_partition_offsets.m,
                                      X_StrideCell, Y_StrideCell, 2);
    }
#endif

    return *this;
}
//#####################################################################
// operator*=
//#####################################################################
template<class POLICY> KRYLOV_VECTOR_BASE<typename POLICY::ELASTICITY::SCALAR >& CG_VECTOR<POLICY>::
operator*=(const T a)
{
#ifdef LOG_DETAILED_PERFORMANCE_CG_VECTOR
    LOG::SCOPE scope("CG_VECTOR::Multiply");
#endif

    // VECTOR<ARRAY<T, T_INDEX>, d> x_array_debug;
    // for( int i = 0; i < d; i++){
    //     x_array_debug(i+1).Resize( ne.padded_node_domain );

    //     for( RANGE_ITERATOR<d> iterator(ne.padded_node_domain);iterator.Valid();iterator.Next()){
    //         x_array_debug(i+1)(iterator.Index()) = x_array(i+1)(iterator.Index());
    //     }
    // }

    // for(int v=1;v<=d;v++) x_array_debug(v)*=a;

#if !defined(USE_FAST_CG_VECTOR_KERNELS)
    for(int v=1;v<=d;v++) state.x(v)*=a;
    state.p*=a;
#else   
    {     
        for(int i=1;i<=d;i++)
            {
                MT_Streaming_Kernels::
                    Vector_Times_Helper<T>(state.x(i).array.Get_Array_Pointer(),
                                           a,
                                           ne.node_block_base_offsets.Get_Array_Pointer(),
                                           ne.node_block_base_offsets.m,
                                           ne.node_block_partition_offsets.Get_Array_Pointer(),
                                           ne.node_block_partition_offsets.m,
                                           X_StrideNode, Y_StrideNode, 3);
            }

        MT_Streaming_Kernels::
            Vector_Times_Helper<T>(state.p.array.Get_Array_Pointer(),
                                   a,
                                   ne.cell_block_base_offsets.Get_Array_Pointer(),
                                   ne.cell_block_base_offsets.m,
                                   ne.cell_block_partition_offsets.Get_Array_Pointer(),
                                   ne.cell_block_partition_offsets.m,
                                   X_StrideCell, Y_StrideCell, 2);        
    }
#endif

    
    //   LOG::cout << std::endl << a << std::endl;

    // for( int i = 0; i < d; i ++)
    //     for( RANGE_ITERATOR<d> iterator(ne.padded_node_domain);iterator.Valid();iterator.Next()){
    //         if( x_array_debug(i+1)(iterator.Index()) != x_array(i+1)(iterator.Index()))
                
    //             LOG::cout << "Index: " << iterator.Index() << " IS ACTIVE: " << ne.node_is_active(iterator.Index() )<< "  Debug: " << x_array_debug(i+1)(iterator.Index()) << "  Test: " << x_array(i+1)(iterator.Index()) << std::endl;
    //     } 



    return *this;
}
//#####################################################################
// operator*=
//#####################################################################
template<class POLICY> KRYLOV_VECTOR_BASE<typename POLICY::ELASTICITY::SCALAR >& CG_VECTOR<POLICY>::
operator*=(const BASE& bv)
{
#ifdef LOG_DETAILED_PERFORMANCE_CG_VECTOR
    LOG::SCOPE scope("CG_VECTOR::Element_Times");
#endif

#if !defined(USE_FAST_CG_VECTOR_KERNELS)
    state.x *= State(bv).x;
    state.p *= State(bv).p;
#else   
    {     
        int number_of_blocks = ne.cgblock_base_offsets.m;
        ARRAY<int> partition_offsets;
        for(int partition=0;partition<ne.constant_partitions;partition++)
            partition_offsets.Append((partition)*(number_of_blocks/ ne.constant_partitions )+min(partition,number_of_blocks% ne.constant_partitions ));
          
        for(int i=1;i<=d;i++)
            {
                MT_Streaming_Kernels::
                    Vector_Element_Times_Helper<T,16>(state.x(i).array.Get_Array_Pointer(),
                                              State(bv).x(i).array.Get_Array_Pointer(),
                                              ne.cg_offsets.Get_Array_Pointer(),
                                              ne.cg_offsets.m,
                                              ne.constant_partitions);
            }

        MT_Streaming_Kernels::
            Vector_Element_Times_Helper<T,16>(state.p.array.Get_Array_Pointer(),
                                      State(bv).p.array.Get_Array_Pointer(),
                                      ne.cg_offsets.Get_Array_Pointer(),
                                      ne.cg_offsets.m,
                                      ne.constant_partitions);

    }
#endif

    return *this;
}
//#####################################################################
// Function Copy
//#####################################################################
template<class POLICY> void CG_VECTOR<POLICY>::
Copy(const T c,const BASE& bv)
{
#ifdef LOG_DETAILED_PERFORMANCE_CG_VECTOR
    LOG::SCOPE scope("CG_VECTOR::SAXPY_A");
#endif

#if !defined(USE_FAST_CG_VECTOR_KERNELS)
    for(int v=1;v<=d;v++) T_STATE::T_SCALAR_VARIABLE::Copy(c,State(bv).x(v),state.x(v));
    T_STATE::T_SCALAR_VARIABLE::Copy(c,State(bv).p,state.p);
#else

    {
        for(int i=1;i<=d;i++)
            {
                MT_Streaming_Kernels::
                    Vector_SAXPY_A_Helper<T>(state.x(i).array.Get_Array_Pointer(),
                                             State(bv).x(i).array.Get_Array_Pointer(),
                                             c,
                                             ne.node_block_base_offsets.Get_Array_Pointer(),
                                             ne.node_block_base_offsets.m,
                                             ne.node_block_partition_offsets.Get_Array_Pointer(),
                                             ne.node_block_partition_offsets.m,
                                             X_StrideNode, Y_StrideNode, 3);
            }

        MT_Streaming_Kernels::
            Vector_SAXPY_A_Helper<T>(state.p.array.Get_Array_Pointer(),
                                     State(bv).p.array.Get_Array_Pointer(),
                                     c,
                                     ne.cell_block_base_offsets.Get_Array_Pointer(),
                                     ne.cell_block_base_offsets.m,
                                     ne.cell_block_partition_offsets.Get_Array_Pointer(),
                                     ne.cell_block_partition_offsets.m,
                                     X_StrideCell, Y_StrideCell, 2);


    }

#endif

}
//#####################################################################
// Function Copy
//#####################################################################
template<class POLICY> void CG_VECTOR<POLICY>::
Copy(const T c1,const BASE& bv1,const BASE& bv2)
{
#ifdef LOG_DETAILED_PERFORMANCE_CG_VECTOR
    LOG::SCOPE scope("CG_VECTOR::SAXPY_AB");
#endif

#if !defined(USE_FAST_CG_VECTOR_KERNELS)
    for(int v=1;v<=d;v++) T_STATE::T_SCALAR_VARIABLE::Copy(c1,State(bv1).x(v),State(bv2).x(v),state.x(v));
    T_STATE::T_SCALAR_VARIABLE::Copy(c1,State(bv1).p,State(bv2).p,state.p);
#else

    {
        for(int i=1;i<=d;i++)
            {
                MT_Streaming_Kernels::
                    Vector_SAXPY_AB_Helper<T>(state.x(i).array.Get_Array_Pointer(),
                                              State(bv1).x(i).array.Get_Array_Pointer(),
                                              c1,
                                              State(bv2).x(i).array.Get_Array_Pointer(),
                                              ne.node_block_base_offsets.Get_Array_Pointer(),
                                              ne.node_block_base_offsets.m,
                                              ne.node_block_partition_offsets.Get_Array_Pointer(),
                                              ne.node_block_partition_offsets.m,
                                              X_StrideNode, Y_StrideNode, 3);
            }
        
        MT_Streaming_Kernels::
            Vector_SAXPY_AB_Helper<T>(state.p.array.Get_Array_Pointer(),
                                      State(bv1).p.array.Get_Array_Pointer(),
                                      c1,
                                      State(bv2).p.array.Get_Array_Pointer(),
                                      ne.cell_block_base_offsets.Get_Array_Pointer(),
                                      ne.cell_block_base_offsets.m,
                                      ne.cell_block_partition_offsets.Get_Array_Pointer(),
                                      ne.cell_block_partition_offsets.m,
                                      X_StrideCell, Y_StrideCell, 2);
        
    }

#endif

}
//#####################################################################
// Function Raw_Size
//#####################################################################
template<class POLICY> int CG_VECTOR<POLICY>::
Raw_Size() const
{
    PHYSBAM_NOT_IMPLEMENTED();
}
//#####################################################################
// Function Raw_Get
//#####################################################################
template<class POLICY> typename POLICY::ELASTICITY::SCALAR& CG_VECTOR<POLICY>::
Raw_Get(int i)
{
    PHYSBAM_NOT_IMPLEMENTED();
}
//#####################################################################
// Function Clear
//#####################################################################
template<class POLICY> void CG_VECTOR<POLICY>::
Clear()
{
#ifdef LOG_DETAILED_PERFORMANCE_CG_VECTOR
    LOG::SCOPE scope("CG_VECTOR::Clear");
#endif

#if !defined(USE_FAST_CG_VECTOR_KERNELS)
    for(int v=1;v<=d;v++) state.x(v).Fill(T());
    state.p.Fill(T());
#else
    {
        for(int i=1;i<=d;i++)
            {
                MT_Streaming_Kernels::
                    Vector_Set_Helper<T>(state.x(i).array.Get_Array_Pointer(),
                                           T(0),
                                           ne.node_block_base_offsets.Get_Array_Pointer(),
                                           ne.node_block_base_offsets.m,
                                           ne.node_block_partition_offsets.Get_Array_Pointer(),
                                           ne.node_block_partition_offsets.m,
                                           X_StrideNode, Y_StrideNode, 3);
            }

        MT_Streaming_Kernels::
            Vector_Set_Helper<T>(state.p.array.Get_Array_Pointer(),
                                   T(0),
                                   ne.cell_block_base_offsets.Get_Array_Pointer(),
                                   ne.cell_block_base_offsets.m,
                                   ne.cell_block_partition_offsets.Get_Array_Pointer(),
                                   ne.cell_block_partition_offsets.m,
                                   X_StrideCell, Y_StrideCell, 2);
    }
#endif
}
//#####################################################################
// Function Set
//#####################################################################
template<class POLICY> void CG_VECTOR<POLICY>::
Set(const T a)
{
#ifdef LOG_DETAILED_PERFORMANCE_CG_VECTOR
    LOG::SCOPE scope("CG_VECTOR::Set");
#endif

#if !defined(USE_FAST_CG_VECTOR_KERNELS)
    for(int v=1;v<=d;v++) state.x(v).Fill(T(a));
    state.p.Fill(T(a));
#else
    {
        for(int i=1;i<=d;i++)
            {
                MT_Streaming_Kernels::
                    Vector_Set_Helper<T>(state.x(i).array.Get_Array_Pointer(),
                                           T(a),
                                           ne.node_block_base_offsets.Get_Array_Pointer(),
                                           ne.node_block_base_offsets.m,
                                           ne.node_block_partition_offsets.Get_Array_Pointer(),
                                           ne.node_block_partition_offsets.m,
                                           X_StrideNode, Y_StrideNode, 3);
            }

        MT_Streaming_Kernels::
            Vector_Set_Helper<T>(state.p.array.Get_Array_Pointer(),
                                   T(a),
                                   ne.cell_block_base_offsets.Get_Array_Pointer(),
                                   ne.cell_block_base_offsets.m,
                                   ne.cell_block_partition_offsets.Get_Array_Pointer(),
                                   ne.cell_block_partition_offsets.m,
                                   X_StrideCell, Y_StrideCell, 2);
    }
#endif
}
//#####################################################################
//template class CG_VECTOR<float,2>;
template class CG_VECTOR<CG_POLICY<SKINNING_NONLINEAR_ELASTICITY<float,3,true,true> > >;
template class CG_VECTOR<CG_POLICY<SKINNING_NONLINEAR_ELASTICITY<float,3,true,false> > >;
template class CG_VECTOR<CG_POLICY<NONLINEAR_ELASTICITY<float,3> > >;
#ifndef USE_SPECIALIZED_KERNELS
template class CG_VECTOR<CG_POLICY<SKINNING_NONLINEAR_ELASTICITY<double,3,true,true> > >;
template class CG_VECTOR<CG_POLICY<SKINNING_NONLINEAR_ELASTICITY<double,3,true,false> > >;
template class CG_VECTOR<CG_POLICY<NONLINEAR_ELASTICITY<double,3> > >;
#endif

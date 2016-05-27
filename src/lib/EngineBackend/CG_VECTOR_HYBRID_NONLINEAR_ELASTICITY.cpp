//#####################################################################
// Copyright 2011, Taylor Patterson, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class CG_VECTOR (Specialization for HYBRID_NONLINEAR_ELASTICITY)
//#####################################################################
#include "CG_VECTOR.h"
#include "HYBRID_NONLINEAR_ELASTICITY.h"

#include <CG_Optimized_Kernels/Vector_Add/Vector_Add_Helper.h>
#include <CG_Optimized_Kernels/Vector_Subtract/Vector_Subtract_Helper.h>
#include <CG_Optimized_Kernels/Vector_Times/Vector_Times_Helper.h>
#include <CG_Optimized_Kernels/Vector_Set/Vector_Set_Helper.h>
#include <CG_Optimized_Kernels/Vector_SAXPY_A/Vector_SAXPY_A_Helper.h>
#include <CG_Optimized_Kernels/Vector_SAXPY_AB/Vector_SAXPY_AB_Helper.h>
#include <CG_Optimized_Kernels/Vector_Element_Times/Vector_Element_Times_Helper.h>

#include <PhysBAM_Tools/Arrays_Computations/ARRAY_COPY.h>

#ifdef LOG_DETAILED_PERFORMANCE
#define LOG_DETAILED_PERFORMANCE_CG_VECTOR
#endif

#define CG_VECTOR_THREADS 12
#define CG_VECTOR_PARTITIONS 12
#define USE_FAST_CG_VECTOR_KERNELS

using namespace PhysBAM;
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
    state.x_mesh+=State(bv).x_mesh;
    state.p_mesh+=State(bv).p_mesh;
#else   
    {
        
        for(int i=1;i<=d;i++)
            {
                MT_Streaming_Kernels::
                    Vector_Add_Helper<T,16>(state.x(i).array.Get_Array_Pointer(),
                                         State(bv).x(i).array.Get_Array_Pointer(),
                                         ne.cg_offsets.Get_Array_Pointer(),
                                         ne.cg_offsets.m,
                                         ne.constant_partitions);


                MT_Streaming_Kernels::
                    Vector_Add_Helper<T,16>(state.x_mesh(i).Get_Array_Pointer(),
                                         State(bv).x_mesh(i).Get_Array_Pointer(),
                                         ne.cg_offsets_node_mesh.Get_Array_Pointer(),
                                         ne.cg_offsets_node_mesh.m,
                                         ne.constant_partitions );
            }
/*
        MT_Streaming_Kernels::
            Vector_Add_Helper<T,16>(state.p.array.Get_Array_Pointer(),
                                 State(bv).p.array.Get_Array_Pointer(),
                                 ne.cg_offsets.Get_Array_Pointer(),
                                 ne.cg_offsets.m,
                                 ne.constant_partitions);

        MT_Streaming_Kernels::
            Vector_Add_Helper<T,16>(state.p_mesh.Get_Array_Pointer(),
                                 State(bv).p_mesh.Get_Array_Pointer(),
                                 ne.cg_offsets_cell_mesh.Get_Array_Pointer(),
                                 ne.cg_offsets_cell_mesh.m,
                                 ne.constant_partitions );      
*/

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
    state.x_mesh-=State(bv).x_mesh;
    state.p_mesh-=State(bv).p_mesh;
#else   
    {     
          
        for(int i=1;i<=d;i++)
            {
                MT_Streaming_Kernels::
                    Vector_Subtract_Helper<T,16>(state.x(i).array.Get_Array_Pointer(),
                                              State(bv).x(i).array.Get_Array_Pointer(),
                                              ne.cg_offsets.Get_Array_Pointer(),
                                              ne.cg_offsets.m,
                                              ne.constant_partitions);

                MT_Streaming_Kernels::
                    Vector_Subtract_Helper<T,16>(state.x_mesh(i).Get_Array_Pointer(),
                                              State(bv).x_mesh(i).Get_Array_Pointer(),
                                              ne.cg_offsets_node_mesh.Get_Array_Pointer(),
                                              ne.cg_offsets_node_mesh.m,
                                              ne.constant_partitions);
            }
/*
        MT_Streaming_Kernels::
            Vector_Subtract_Helper<T,16>(state.p.array.Get_Array_Pointer(),
                                      State(bv).p.array.Get_Array_Pointer(),
                                      ne.cg_offsets.Get_Array_Pointer(),
                                      ne.cg_offsets.m,
                                      ne.constant_partitions);

        MT_Streaming_Kernels::
            Vector_Subtract_Helper<T,16>(state.p_mesh.Get_Array_Pointer(),
                                      State(bv).p_mesh.Get_Array_Pointer(),
                                      ne.cg_offsets_cell_mesh.Get_Array_Pointer(),
                                      ne.cg_offsets_cell_mesh.m,
                                      ne.constant_partitions);
*/
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

#if !defined(USE_FAST_CG_VECTOR_KERNELS)
    for(int v=1;v<=d;v++) state.x(v)*=a;
    state.p*=a;
    for(int v=1;v<=d;v++) state.x_mesh(v)*=a;
    state.p_mesh*=a;
#else   
    {     
       
        for(int i=1;i<=d;i++)
            {
                MT_Streaming_Kernels::
                    Vector_Times_Helper<T,16>(state.x(i).array.Get_Array_Pointer(),
                                           a,
                                           ne.cg_offsets.Get_Array_Pointer(),
                                           ne.cg_offsets.m,
                                           ne.constant_partitions);
                
                MT_Streaming_Kernels::
                    Vector_Times_Helper<T,16>(state.x_mesh(i).Get_Array_Pointer(),
                                           a,
                                           ne.cg_offsets_node_mesh.Get_Array_Pointer(),
                                           ne.cg_offsets_node_mesh.m,
                                           ne.constant_partitions);
            }
/*
        MT_Streaming_Kernels::
            Vector_Times_Helper<T,16>(state.p.array.Get_Array_Pointer(),
                                   a,
                                   ne.cg_offsets.Get_Array_Pointer(),
                                   ne.cg_offsets.m,
                                   ne.constant_partitions);
        
        MT_Streaming_Kernels::
            Vector_Times_Helper<T,16>(state.p_mesh.Get_Array_Pointer(),
                                   a,
                                   ne.cg_offsets_cell_mesh.Get_Array_Pointer(),
                                   ne.cg_offsets_cell_mesh.m,
                                   ne.constant_partitions);
*/
    }
#endif

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
    state.x_mesh*=State(bv).x_mesh;
    state.p_mesh*=State(bv).p_mesh;
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

                MT_Streaming_Kernels::
                    Vector_Element_Times_Helper<T,16>(state.x_mesh(i).Get_Array_Pointer(),
                                              State(bv).x_mesh(i).Get_Array_Pointer(),
                                              ne.cg_offsets_node_mesh.Get_Array_Pointer(),
                                              ne.cg_offsets_node_mesh.m,
                                              ne.constant_partitions);
            }
/*
        MT_Streaming_Kernels::
            Vector_Element_Times_Helper<T,16>(state.p.array.Get_Array_Pointer(),
                                      State(bv).p.array.Get_Array_Pointer(),
                                      ne.cg_offsets.Get_Array_Pointer(),
                                      ne.cg_offsets.m,
                                      ne.constant_partitions);

        MT_Streaming_Kernels::
            Vector_Element_Times_Helper<T,16>(state.p_mesh.Get_Array_Pointer(),
                                      State(bv).p_mesh.Get_Array_Pointer(),
                                      ne.cg_offsets_cell_mesh.Get_Array_Pointer(),
                                      ne.cg_offsets_cell_mesh.m,
                                      ne.constant_partitions);
*/
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
    for(int v=1;v<=d;v++) T_STATE::T_SCALAR_VARIABLE_MESH::Copy(c,State(bv).x_mesh(v),state.x_mesh(v));
    T_STATE::T_SCALAR_VARIABLE_MESH::Copy(c,State(bv).p_mesh,state.p_mesh);
#else

    {
      
        for(int i=1;i<=d;i++)
            {
                MT_Streaming_Kernels::
                    Vector_SAXPY_A_Helper<T,16>(state.x(i).array.Get_Array_Pointer(),
                                             State(bv).x(i).array.Get_Array_Pointer(),
                                             c,
                                             ne.cg_offsets.Get_Array_Pointer(),
                                             ne.cg_offsets.m,
                                             ne.constant_partitions);

                MT_Streaming_Kernels::
                    Vector_SAXPY_A_Helper<T,16>(state.x_mesh(i).Get_Array_Pointer(),
                                             State(bv).x_mesh(i).Get_Array_Pointer(),
                                             c,
                                             ne.cg_offsets_node_mesh.Get_Array_Pointer(),
                                             ne.cg_offsets_node_mesh.m,
                                             ne.constant_partitions);               
            }
/*
        MT_Streaming_Kernels::
            Vector_SAXPY_A_Helper<T,16>(state.p.array.Get_Array_Pointer(),
                                     State(bv).p.array.Get_Array_Pointer(),
                                     c,
                                     ne.cg_offsets.Get_Array_Pointer(),
                                     ne.cg_offsets.m,
                                     ne.constant_partitions);

        MT_Streaming_Kernels::
            Vector_SAXPY_A_Helper<T,16>(state.p_mesh.Get_Array_Pointer(),
                                     State(bv).p_mesh.Get_Array_Pointer(),
                                     c,
                                     ne.cg_offsets_cell_mesh.Get_Array_Pointer(),
                                     ne.cg_offsets_cell_mesh.m,
                                     ne.constant_partitions);
*/

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
    for(int v=1;v<=d;v++) T_STATE::T_SCALAR_VARIABLE_MESH::Copy(c1,State(bv1).x_mesh(v),State(bv2).x_mesh(v),state.x_mesh(v));
    T_STATE::T_SCALAR_VARIABLE_MESH::Copy(c1,State(bv1).p_mesh,State(bv2).p_mesh,state.p_mesh);
#else

    {
       for(int i=1;i<=d;i++)
            {
                MT_Streaming_Kernels::
                    Vector_SAXPY_AB_Helper<T,16>(state.x(i).array.Get_Array_Pointer(),
                                                 State(bv1).x(i).array.Get_Array_Pointer(),
                                                 c1,
                                                 State(bv2).x(i).array.Get_Array_Pointer(),
                                                 ne.cg_offsets.Get_Array_Pointer(),
                                                 ne.cg_offsets.m,
                                                 ne.constant_partitions);
                
                MT_Streaming_Kernels::
                    Vector_SAXPY_AB_Helper<T,16>(state.x_mesh(i).Get_Array_Pointer(),
                                              State(bv1).x_mesh(i).Get_Array_Pointer(),
                                              c1,
                                              State(bv2).x_mesh(i).Get_Array_Pointer(),
                                              ne.cg_offsets_node_mesh.Get_Array_Pointer(),
                                              ne.cg_offsets_node_mesh.m,
                                              ne.constant_partitions);
                
                
            }      
/*
        MT_Streaming_Kernels::
            Vector_SAXPY_AB_Helper<T,16>(state.p.array.Get_Array_Pointer(),
                                         State(bv1).p.array.Get_Array_Pointer(),
                                      c1,
                                      State(bv2).p.array.Get_Array_Pointer(),
                                      ne.cg_offsets.Get_Array_Pointer(),
                                      ne.cg_offsets.m,
                                      ne.constant_partitions);

        MT_Streaming_Kernels::
            Vector_SAXPY_AB_Helper<T,16>(state.p_mesh.Get_Array_Pointer(),
                                      State(bv1).p_mesh.Get_Array_Pointer(),
                                      c1,
                                      State(bv2).p_mesh.Get_Array_Pointer(),
                                      ne.cg_offsets_cell_mesh.Get_Array_Pointer(),
                                      ne.cg_offsets_cell_mesh.m,
                                      ne.constant_partitions);
*/
        
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
    //for(int v=1;v<=d;v++) state.x_mesh(v).Fill(T());
    //state.p_mesh.Fill(T());
    for(int v=1;v<=d;v++) for(int i(1);i<=state.x_mesh(v).m;i++) state.x_mesh(v)(i)=0.0;
    for(int i(1);i<=state.p_mesh.m;i++) state.p_mesh(i)=0.0;

#else
    {
       
        for(int i=1;i<=d;i++)
            {
                #if 1
                MT_Streaming_Kernels::
                    Vector_Set_Helper<T,16>(state.x(i).array.Get_Array_Pointer(),
                                            T(0),
                                            ne.cg_offsets.Get_Array_Pointer(),
                                            ne.cg_offsets.m,
                                            ne.constant_partitions);
                #else
                state.x(i).Fill(T());
                #endif              
                
                #if 1
                MT_Streaming_Kernels::
                    Vector_Set_Helper<T,16>(state.x_mesh(i).Get_Array_Pointer(),
                                           T(0),
                                            ne.cg_offsets_node_mesh.Get_Array_Pointer(),
                                            ne.cg_offsets_node_mesh.m,
                                            ne.constant_partitions);
                #else
                for(int n(1);n<=state.x_mesh(i).m;n++) state.x_mesh(i)(n)=0.0;
                #endif
            }
/*
        #if 1
        MT_Streaming_Kernels::
            Vector_Set_Helper<T,16>(state.p.array.Get_Array_Pointer(),
                                    T(0),
                                    ne.cg_offsets.Get_Array_Pointer(),
                                    ne.cg_offsets.m,
                                    ne.constant_partitions);
        #else
        state.p.Fill(T());
        #endif
        
        #if 1
        MT_Streaming_Kernels::
            Vector_Set_Helper<T,16>(state.p_mesh.Get_Array_Pointer(),
                                   T(0),
                                    ne.cg_offsets_cell_mesh.Get_Array_Pointer(),
                                    ne.cg_offsets_cell_mesh.m,
                                    ne.constant_partitions);
        #else
        for(int n(1);n<=state.p_mesh.m;n++) state.p_mesh(n)=0.0;
        #endif
*/
    }
#endif

}
//#####################################################################
// Function Set
//#####################################################################
template<class POLICY> void CG_VECTOR<POLICY>::
Set( const T a)
{
#ifdef LOG_DETAILED_PERFORMANCE_CG_VECTOR
    LOG::SCOPE scope("CG_VECTOR::Set");
#endif

#if !defined(USE_FAST_CG_VECTOR_KERNELS)
    for(int v=1;v<=d;v++) state.x(v).Fill(T(a));
    state.p.Fill(T(a));
    //for(int v=1;v<=d;v++) state.x_mesh(v).Fill(T());
    //state.p_mesh.Fill(T());
    for(int v=1;v<=d;v++) for(int i(1);i<=state.x_mesh(v).m;i++) state.x_mesh(v)(i)=a;
    for(int i(1);i<=state.p_mesh.m;i++) state.p_mesh(i)=a;

#else
    {
       
        for(int i=1;i<=d;i++)
            {
                #if 1
                MT_Streaming_Kernels::
                    Vector_Set_Helper<T,16>(state.x(i).array.Get_Array_Pointer(),
                                            T(a),
                                            ne.cg_offsets.Get_Array_Pointer(),
                                            ne.cg_offsets.m,
                                            ne.constant_partitions);
                #else
                state.x(i).Fill(T(a));
                #endif              
                
                #if 1
                MT_Streaming_Kernels::
                    Vector_Set_Helper<T,16>(state.x_mesh(i).Get_Array_Pointer(),
                                            T(a),
                                            ne.cg_offsets_node_mesh.Get_Array_Pointer(),
                                            ne.cg_offsets_node_mesh.m,
                                            ne.constant_partitions);
                #else
                for(int n(1);n<=state.x_mesh(i).m;n++) state.x_mesh(i)(n)=a;
                #endif
            }
/*
        #if 1
        MT_Streaming_Kernels::
            Vector_Set_Helper<T,16>(state.p.array.Get_Array_Pointer(),
                                    T(a),
                                    ne.cg_offsets.Get_Array_Pointer(),
                                    ne.cg_offsets.m,
                                    ne.constant_partitions);
        #else
        state.p.Fill(T(a));
        #endif
        
        #if 1
        MT_Streaming_Kernels::
            Vector_Set_Helper<T,16>(state.p_mesh.Get_Array_Pointer(),
                                    T(a),
                                    ne.cg_offsets_cell_mesh.Get_Array_Pointer(),
                                    ne.cg_offsets_cell_mesh.m,
                                    ne.constant_partitions);
        #else
        for(int n(1);n<=state.p_mesh.m;n++) state.p_mesh(n)=a;
        #endif
*/
    }
#endif

}
//#####################################################################
//template class CG_VECTOR<float,2>;
template class CG_VECTOR<CG_POLICY<HYBRID_NONLINEAR_ELASTICITY<SKINNING_NONLINEAR_ELASTICITY<float,3,true,true> > > >;
template class CG_VECTOR<CG_POLICY<HYBRID_NONLINEAR_ELASTICITY<SKINNING_NONLINEAR_ELASTICITY<float,3,true,false> > > >;
template class CG_VECTOR<CG_POLICY<HYBRID_NONLINEAR_ELASTICITY<NONLINEAR_ELASTICITY<float,3> > > >;
#ifndef USE_SPECIALIZED_KERNELS
template class CG_VECTOR<CG_POLICY<HYBRID_NONLINEAR_ELASTICITY<SKINNING_NONLINEAR_ELASTICITY<double,3,true,true> > > >;
template class CG_VECTOR<CG_POLICY<HYBRID_NONLINEAR_ELASTICITY<SKINNING_NONLINEAR_ELASTICITY<double,3,true,false> > > >;
template class CG_VECTOR<CG_POLICY<HYBRID_NONLINEAR_ELASTICITY<NONLINEAR_ELASTICITY<double,3> > > >;
#endif

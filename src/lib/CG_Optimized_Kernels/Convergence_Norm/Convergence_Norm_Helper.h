//#####################################################################
// Copyright 2010, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __Convergence_Norm_Helper__
#define __Convergence_Norm_Helper__
#include <algorithm>
#include <cmath>
#include <iostream>
#include "../Kernel_Base/Kernel_Base_Reducer_Helper.h"
#include "../Kernel_Base/Kernel_Mesh_Base_Reducer_Helper.h"
#include "../Kernel_Base/Kernel_Serial_Base_Reducer_Helper.h"
#include "../Kernel_Base/Kernel_OpenMP_Base_Reducer_Helper.h"

namespace MT_Streaming_Kernels{


template <class T, int stride=1, bool USING_D=false> 
class Vector_Convergence_Norm_Op
    {
    private:
        const T* const x;
        const int* const d;

    public:
        explicit Vector_Convergence_Norm_Op(const T* const x_input, const int* const d_input) : x(x_input), d(d_input) {}
        void Execute(int index, double& result)
        {  
            // NOT IMPLEMENTED
            assert(false);
/*
            T f[stride];
            T _d, _x;
            if(USING_D){               
                for(int i=0;i<stride;i++){
                    _d = (T)d[index+i];
                    _x = (T)x[index+i];
                    f[i] = _d * _x;
                    f[i] = f[i] > 0 ? f[i] : -f[i];
                }

                for(int i=0;i<stride;i++){
                    result = result > f[i] ? result : f[i];
                }

            }
            else{    
                for(int i=0;i<stride;i++){
                    f[i] = (T)x[index+i];
                    f[i] = f[i] > 0 ? f[i] : -f[i];
                }

                for(int i=0;i<stride;i++){
                    result = result > f[i] ? result : f[i];
                }
            }
*/
        }
        void Join(double& result, const double& partial_result)
        {
            result=std::max(result, fabs(partial_result));
        }
    };

 template <class T, bool USING_D> 
     class Vector_Convergence_Norm_Op<T,1,USING_D>
    {
    private:
        const T* const x;
        const int* const d;

    public:
        explicit Vector_Convergence_Norm_Op(const T* const x_input, const int* const d_input) : x(x_input), d(d_input) {}
        void Execute(int index, double& result)
        {
            double f;
            if(USING_D){
                f = (((T)d[index]) * ((T)x[index]));
                f = f > 0 ? f : -f;
                result = result > f ? result : f;
            }
            else{
                f = (T)x[index];
                f = f > 0 ? f : -f;
                result = result > f ? result : f;
            }
        }
        void Join(double& result, const double& partial_result)
        {
            result=std::max(result, fabs(partial_result));
        }
    };


 template <class T, bool using_d> inline
     double Vector_Convergence_Norm_Helper( const T* const x_input,const int* const d_input,
                               const int* block_offsets_input,const int number_of_blocks_input,
                               const int* partition_offsets_input, const int number_of_partitions_input,
                               const int x_stride_input, const int y_stride_input, const int block_size_input)
     {
         Vector_Convergence_Norm_Op<T,1,using_d> op(x_input, d_input);
         Kernel_Base_Reducer_Helper<Vector_Convergence_Norm_Op<T,1,using_d> > helper(op, block_offsets_input,
                                                                           number_of_blocks_input,
                                                                           partition_offsets_input,
                                                                           number_of_partitions_input,
                                                                           x_stride_input,
                                                                           y_stride_input,
                                                                           block_size_input);
         return helper.Run_Parallel();
     };


 template <class T, class LAYOUT, bool using_d> inline
     double Vector_Convergence_Norm_Helper( const T* const x_input,const int* const d_input,
                               const LAYOUT* block_layouts_input,const int number_of_blocks_input,
                               const int* partition_offsets_input, const int number_of_partitions_input,
                               const int x_stride_input, const int y_stride_input, const int block_size_input)
     {
         Vector_Convergence_Norm_Op<T,1,using_d> op(x_input, d_input);
         Kernel_Mesh_Base_Reducer_Helper<Vector_Convergence_Norm_Op<T,1,using_d>, LAYOUT > helper(op, block_layouts_input,
                                                                           number_of_blocks_input,
                                                                           partition_offsets_input,
                                                                           number_of_partitions_input,
                                                                           x_stride_input,
                                                                           y_stride_input,
                                                                           block_size_input);
         return helper.Run_Parallel();
     };


 template <class T, bool using_d> inline
     double Vector_Convergence_Norm_Helper( const T* const x_input,const int* const d_input,
                                            const int number_of_entries_input,
                                            const int number_of_partitions_input)
     {
         Vector_Convergence_Norm_Op<T,1,using_d> op(x_input, d_input);
         Kernel_Serial_Base_Reducer_Helper<Vector_Convergence_Norm_Op<T,1,using_d> > helper(op, number_of_entries_input,
                                                                                  number_of_partitions_input);
         return helper.Run_Parallel();
     };

#include <omp.h>

 template <class T, int stride, bool using_d> inline
     double Vector_Convergence_Norm_Helper(const T* const x_input,const int* const d_input,
                                           const int* offsets_input,
                                           const int number_of_offsets_input,
                                           const int number_of_partitions_input)
 {

     if(number_of_offsets_input == 0)
         return 0.0;

#if 0
     Vector_Convergence_Norm_Op<T,stride,using_d> op(x_input, d_input);
     Kernel_OpenMP_Base_Reducer_Helper<Vector_Convergence_Norm_Op<T,stride,using_d> >helper(op, offsets_input, 
                                                                            number_of_offsets_input,
                                                                            number_of_partitions_input);
     return helper.Run_Parallel();
#else

     
#ifdef __MIC__
     assert(stride==16);
     double result=1e-10;
     __m512 t_result_p[number_of_partitions_input];

     int partition_offsets[number_of_partitions_input];
     int partition;
     partition_offsets[0]=0;
     for( int i=0, partition=0; i < number_of_offsets_input && partition < number_of_partitions_input;
          i+=number_of_offsets_input / number_of_partitions_input, partition++ )
         partition_offsets[partition] = i;
     

     if(using_d){
#pragma omp parallel num_threads(number_of_partitions_input) shared(partition_offsets)
         {
             const int tid=omp_get_thread_num();
             const int begin = partition_offsets[tid];
             const int end = ((tid<number_of_partitions_input-1)?partition_offsets[tid+1]:number_of_offsets_input);
             
             t_result_p[tid] = _mm512_set1_ps(1e-10);

             for( int i=begin; i<end; i++){
                 
                 _mm_prefetch((char*)&x_input[offsets_input[i+64]],_MM_HINT_T1);
                 _mm_prefetch((char*)&x_input[offsets_input[i+8]],_MM_HINT_T0);
                 __m512 rX=_mm512_setzero_ps();
                 rX=_mm512_load_ps(&x_input[offsets_input[i]]);
                 
                 _mm_prefetch((char*)&d_input[offsets_input[i+64]],_MM_HINT_T1);
                 _mm_prefetch((char*)&d_input[offsets_input[i+8]],_MM_HINT_T0);
                 __m512i rD=_mm512_setzero_epi32();
                 rD=_mm512_load_epi32(&d_input[offsets_input[i]]);
                 
                 
                 __m512i rZ =_mm512_setzero_epi32();
                 __mmask16 rDm=_mm512_cmp_epi32_mask(rD, rZ,  _MM_CMPINT_NE);
                 t_result_p[tid] = _mm512_mask_gmaxabs_ps(t_result_p[tid], rDm, t_result_p[tid], rX);
                     
             }
         }         
     }
      else{
#pragma omp parallel num_threads(number_of_partitions_input) shared(partition_offsets)
         {
             const int tid=omp_get_thread_num();
             const int begin = partition_offsets[tid];
             const int end = ((tid<number_of_partitions_input-1)?partition_offsets[tid+1]:number_of_offsets_input);
             
             t_result_p[tid] = _mm512_set1_ps(1e-10);
             
             for( int i=begin; i<end; i++){
                 
                 _mm_prefetch((char*)&x_input[offsets_input[i+64]],_MM_HINT_T1);
                 _mm_prefetch((char*)&x_input[offsets_input[i+8]],_MM_HINT_T0);
                 __m512 rX=_mm512_setzero_ps();
                 rX=_mm512_load_ps(&x_input[offsets_input[i]]);
                                
                 t_result_p[tid] = _mm512_gmaxabs_ps(t_result_p[tid], rX);
                 
             }
         }
      }
     
     result=_mm512_reduce_max_ps(t_result_p[0]);
     for(int j=1; j<number_of_partitions_input;j++){
         double r = _mm512_reduce_max_ps(t_result_p[j]);
         result = result > r ? result : r;
     }
     
     return result;
     

#else
     double t_result[stride];
     memset(t_result, 1e-10, stride*sizeof(double));
     double result=1e-10;
     double t_result_p[number_of_partitions_input][stride];
     memset(t_result_p, 1e-10, number_of_partitions_input*stride*sizeof(double));

     int partition_offsets[number_of_partitions_input];
     int partition;
     partition_offsets[0]=0;
     for( int i=0, partition=0; i < number_of_offsets_input && partition < number_of_partitions_input;
          i+=number_of_offsets_input / number_of_partitions_input, partition++ )
         partition_offsets[partition] = i;
     

     if(using_d){
#pragma omp parallel num_threads(number_of_partitions_input) shared(partition_offsets)
         {
             const int tid=omp_get_thread_num();
             const int begin = partition_offsets[tid];
             const int end = ((tid<number_of_partitions_input-1)?partition_offsets[tid+1]:number_of_offsets_input);
             
             for( int i=begin; i<end; i++){
                 int offset = offsets_input[i];
                 for( int k=0; k<stride; k++){
                     double r = d_input[offset+k]*x_input[offset+k];
                     assert( r==r );
                     r = r > 0 ? r : -r;
                     t_result_p[tid][k] = r > t_result_p[tid][k] ? r : t_result_p[tid][k];
                 }
             }             
         }
     }
     else{
#pragma omp parallel num_threads(number_of_partitions_input) shared(partition_offsets)
         {
             const int tid=omp_get_thread_num();
             const int begin = partition_offsets[tid];
             const int end = ((tid<number_of_partitions_input-1)?partition_offsets[tid+1]:number_of_offsets_input);
             
             for( int i=begin; i<end; i++){
                 int offset = offsets_input[i];
                 for( int k=0; k<stride; k++){
                     double r = x_input[offset+k];
                     assert( r==r );
                     r = r > 0 ? r : -r;
                     t_result_p[tid][k] = r > t_result_p[tid][k] ? r : t_result_p[tid][k];
                 }
             }             
         }
     }

     result=t_result_p[0][0];
     for(int j=0; j<number_of_partitions_input;j++)
         for(int i=1;i<stride;i++)
             result = result > t_result_p[j][i] ? result : t_result_p[j][i];
     
     return result;

#endif 



#endif
 }

     
}
#endif

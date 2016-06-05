//#####################################################################
// Copyright 2010, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __Dot_Product_Helper__
#define __Dot_Product_Helper__

#include "../Kernel_Base/Kernel_Base_Reducer_Helper.h"
#include "../Kernel_Base/Kernel_Mesh_Base_Reducer_Helper.h"
#include "../Kernel_Base/Kernel_Serial_Base_Reducer_Helper.h"
#include "../Kernel_Base/Kernel_OpenMP_Base_Reducer_Helper.h"

#ifdef __MIC__
#include <immintrin.h>
#endif


namespace MT_Streaming_Kernels{

    template <class T, int stride=1> 
class Vector_Dot_Product_Op
    {
    private:
        const T* const x;
        const T* const y;


    public:
        explicit Vector_Dot_Product_Op(const T* const x_input,const T* const y_input) : x(x_input), y(y_input) {}
        void Execute(int index, double& result)
        {
            // NOT IMPLEMENTED 
            assert(false);
/*
            T f[stride];
#pragma vector always
#pragma ivdep
            for(int i=0;i<stride;i++)
                f[i] =  (double)x[index+i] * (double)y[index+i];

            T r=0;           
#ifdef __MIC__
            if(stride==16)
                {
                    


                }
            else
#endif
#pragma vector always
#pragma ivdep
                for(int i=0;i<stride;i++)
                    r += f[i];

            result +=r;
*/
        }
        void Join(double& result, const double& partial_result)
        {
            result += partial_result;
        }
    };

template <class T> 
    class Vector_Dot_Product_Op<T,1>
    {
    private:
        const T* const x;
        const T* const y;


    public:
        explicit Vector_Dot_Product_Op(const T* const x_input,const T* const y_input) : x(x_input), y(y_input) {}
        void Execute(int index, double& result)
        {
            result += (double)x[index] * (double)y[index];
        }
        void Join(double& result, const double& partial_result)
        {
            result += partial_result;
        }
    };


 template <class T> inline
     double Vector_Dot_Product_Helper( const T* const x_input,const T* const y_input,
                               const int* block_offsets_input,const int number_of_blocks_input,
                               const int* partition_offsets_input, const int number_of_partitions_input,
                               const int x_stride_input, const int y_stride_input, const int block_size_input)
     {
         Vector_Dot_Product_Op<T,1> op(x_input, y_input);
         Kernel_Base_Reducer_Helper<Vector_Dot_Product_Op<T,1> > helper(op, block_offsets_input,
                                                              number_of_blocks_input,
                                                              partition_offsets_input,
                                                              number_of_partitions_input,
                                                              x_stride_input,
                                                              y_stride_input,
                                                              block_size_input);
             return helper.Run_Parallel();
     };
       
 template <class T, class LAYOUT> inline
     double Vector_Dot_Product_Helper( const T* const x_input,const T* const y_input,
                               const LAYOUT* block_layouts_input,const int number_of_blocks_input,
                               const int* partition_offsets_input, const int number_of_partitions_input,
                               const int x_stride_input, const int y_stride_input, const int block_size_input)
     {
         Vector_Dot_Product_Op<T,1> op(x_input, y_input);
         Kernel_Mesh_Base_Reducer_Helper<Vector_Dot_Product_Op<T,1>, LAYOUT > helper(op, block_layouts_input,
                                                              number_of_blocks_input,
                                                              partition_offsets_input,
                                                              number_of_partitions_input,
                                                              x_stride_input,
                                                              y_stride_input,
                                                              block_size_input);
         return helper.Run_Parallel();
     };

 template <class T> inline
     double Vector_Dot_Product_Helper( const T* const x_input,const T* const y_input,
                                       const int number_of_entries_input,
                                       const int number_of_partitions_input)

     {
         Vector_Dot_Product_Op<T,1> op(x_input, y_input);
         Kernel_Serial_Base_Reducer_Helper<Vector_Dot_Product_Op<T,1> > helper(op, number_of_entries_input,
                                                                             number_of_partitions_input);
             return helper.Run_Parallel();
     };

#include <omp.h>

 template <class T, int stride> inline
     double Vector_Dot_Product_Helper(const T* const x_input,const T* const y_input,
                                      const int* offsets_input,
                                      const int number_of_offsets_input,
                                      const int number_of_partitions_input)
 {

     if(number_of_offsets_input == 0)
         return 0.0;

#if 0
     Vector_Dot_Product_Op<T,stride> op(x_input, y_input);
     Kernel_OpenMP_Base_Reducer_Helper<Vector_Dot_Product_Op<T,stride> >helper(op, offsets_input, 
                                                                            number_of_offsets_input,
                                                                            number_of_partitions_input);
     return helper.Run_Parallel();
#else


   
#ifdef __MIC__
     assert(stride==16);     
 
     double t_result_p[number_of_partitions_input];
     memset(t_result_p, 0, number_of_partitions_input*sizeof(double));

     int partition_offsets[number_of_partitions_input];
     int partition;
     partition_offsets[0]=0;
     for( int i=0, partition=0; i < number_of_offsets_input && partition < number_of_partitions_input;
          i+=number_of_offsets_input / number_of_partitions_input, partition++ )
         partition_offsets[partition] = i;
     

#pragma omp parallel num_threads(number_of_partitions_input) shared(partition_offsets)
         {
             const int tid=omp_get_thread_num();
             const int begin = partition_offsets[tid];
             const int end = ((tid<number_of_partitions_input-1)?partition_offsets[tid+1]:number_of_offsets_input);
             
             __m512d t_resultA_p = _mm512_setzero_pd();
             __m512d t_resultB_p = _mm512_setzero_pd();

             for( int i=begin; i<end; i++){
                 _mm_prefetch((char*)&x_input[offsets_input[i+64]],_MM_HINT_T1);
                 _mm_prefetch((char*)&x_input[offsets_input[i+8]],_MM_HINT_T0);
                 __m512 rX=_mm512_setzero_ps();
                 rX=_mm512_load_ps(&x_input[offsets_input[i]]);
                 
                 _mm_prefetch((char*)&y_input[offsets_input[i+64]],_MM_HINT_T1);
                 _mm_prefetch((char*)&y_input[offsets_input[i+8]],_MM_HINT_T0);
                 __m512 rY=_mm512_setzero_ps();
                 rY=_mm512_load_ps(&y_input[offsets_input[i]]);

                 __m512d rXd = _mm512_cvtpslo_pd(rX);
                 __m512d rYd = _mm512_cvtpslo_pd(rY);
                 t_resultA_p = _mm512_fmadd_pd(rXd, rYd, t_resultA_p);

                 rXd = _mm512_cvtpslo_pd(_mm512_permute4f128_ps(rX,_MM_PERM_DDDC));
                 rYd = _mm512_cvtpslo_pd(_mm512_permute4f128_ps(rY,_MM_PERM_DDDC));
                 t_resultB_p = _mm512_fmadd_pd(rXd, rYd, t_resultB_p);
             }
             double local_result=_mm512_reduce_add_pd(t_resultA_p);
             local_result+=_mm512_reduce_add_pd(t_resultB_p);
             t_result_p[tid]=local_result;
         }

         double result=0;
         for(int i=0; i<number_of_partitions_input; i++)
             result += t_result_p[i];    
         return result;
#else

     double t_result[stride];
     memset(t_result, 0, stride*sizeof(double));
     double t_result_p[number_of_partitions_input];
     memset(t_result_p, 0, number_of_partitions_input*sizeof(double));

     double result=0;

     int partition_offsets[number_of_partitions_input];
     int partition;
     partition_offsets[0]=0;
     for( int i=0, partition=0; i < number_of_offsets_input && partition < number_of_partitions_input;
          i+=number_of_offsets_input / number_of_partitions_input, partition++ )
         partition_offsets[partition] = i;
     

#pragma omp parallel num_threads(number_of_partitions_input) shared(partition_offsets)
         {
             const int tid=omp_get_thread_num();
             const int begin = partition_offsets[tid];
             const int end = ((tid<number_of_partitions_input-1)?partition_offsets[tid+1]:number_of_offsets_input);
             
             for( int i=begin; i<end; i++){
                 int offset = offsets_input[i];
                 for( int k=0; k<stride; k++){
                     double r = x_input[offset+k]*y_input[offset+k];
                     t_result_p[tid] += r;
                 }              
             }             
         }

         result=0;
         for(int i=0; i<number_of_partitions_input; i++)
                 result += t_result_p[i];    
     return result;


#endif 




#endif
 }
                        
}


#endif

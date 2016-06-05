//#####################################################################
// Copyright 2010, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __Vector_Times_Helper__
#define __Vector_Times_Helper__

#include "../Kernel_Base/Kernel_Base_Helper.h"
#include "../Kernel_Base/Kernel_Mesh_Base_Helper.h"
#include "../Kernel_Base/Kernel_Serial_Base_Helper.h"
#include "../Kernel_Base/Kernel_OpenMP_Base_Helper.h"

namespace MT_Streaming_Kernels{

    template <class T, int stride=1> 
class Vector_Times_Op
    {
    private:
        T* const x;
        const T a;
    public:
        explicit Vector_Times_Op(T* const x_input, T a_input) : x(x_input), a(a_input) {}
        void Execute(int index)
        {
#pragma vector always
#pragma ivdep
            for(int i=0;i<stride;i++)
                x[index+i] *= a;
        }
    };
template <class T> 
    class Vector_Times_Op<T,1>
    {
    private:
        T* const x;
        const T a;
    public:
        explicit Vector_Times_Op(T* const x_input, T a_input) : x(x_input), a(a_input) {}
        void Execute(int index)
        {
            x[index] *= a;
        }
    };


 template <class T> inline
     void Vector_Times_Helper( T* const x_input, T a_input,
                               const int* block_offsets_input,const int number_of_blocks_input,
                               const int* partition_offsets_input, const int number_of_partitions_input,
                               const int x_stride_input, const int y_stride_input, const int block_size_input)
     {
         Vector_Times_Op<T> op(x_input, a_input);
         Kernel_Base_Helper<Vector_Times_Op<T> > helper(op, block_offsets_input,  number_of_blocks_input,
                                                        partition_offsets_input,  number_of_partitions_input,
                                                        x_stride_input, y_stride_input, block_size_input);
         helper.Run_Parallel();
     };
                               
 template <class T, class LAYOUT> inline
     void Vector_Times_Helper( T* const x_input, T a_input,
                               const LAYOUT* block_layouts_input,const int number_of_blocks_input,
                               const int* partition_offsets_input, const int number_of_partitions_input,
                               const int x_stride_input, const int y_stride_input, const int block_size_input)
     {
         Vector_Times_Op<T> op(x_input, a_input);
         Kernel_Mesh_Base_Helper<Vector_Times_Op<T>,LAYOUT > helper(op, block_layouts_input,  number_of_blocks_input,
                                                        partition_offsets_input,  number_of_partitions_input,
                                                        x_stride_input, y_stride_input, block_size_input);
         helper.Run_Parallel();
     };

 template <class T> inline
     void Vector_Times_Helper( T* const x_input, T a_input,
                               const int number_of_entries_input,
                               const int number_of_partitions_input)
     {
         Vector_Times_Op<T> op(x_input, a_input);
         Kernel_Serial_Base_Helper<Vector_Times_Op<T> > helper(op, number_of_entries_input,
                                                               number_of_partitions_input);
         helper.Run_Parallel();
     };

 template <class T, int stride> inline
     void Vector_Times_Helper(T* const x_input,  T a_input,
                            const int* offsets_input,
                            const int number_of_offsets_input,
                            const int number_of_partitions_input)
 {
     if(number_of_offsets_input == 0)
         return;

#if 0
     Vector_Times_Op<T,stride> op(x_input, a_input);
     Kernel_OpenMP_Base_Helper<Vector_Times_Op<T,stride> >helper(op, offsets_input, 
                                                               number_of_offsets_input,
                                                               number_of_partitions_input);
     helper.Run_Parallel();
#else
     

#ifdef __MIC__
     enum {width=stride};
     assert(width==16);
     __m512 rC=_mm512_set1_ps(a_input);
     
#pragma omp parallel for num_threads(number_of_partitions_input)
     for(int i=0;i<number_of_offsets_input;i++){
         _mm_prefetch((char*)&x_input[offsets_input[i+64]],_MM_HINT_T1);
         _mm_prefetch((char*)&x_input[offsets_input[i+8]],_MM_HINT_T0);
         __m512 rX=_mm512_setzero_ps();
         rX=_mm512_load_ps(&x_input[offsets_input[i]]);
      
         rX=_mm512_mul_ps(rC,rX);
         _mm512_storenrngo_ps(&x_input[offsets_input[i]],rX);
     }
#else
         
#pragma omp parallel for num_threads(number_of_partitions_input)
     for(int i=0;i<number_of_offsets_input;i++){
         int offset = offsets_input[i];
         for(int j=0;j<stride;j++)
             x_input[offset+j] *= a_input;
    }

#endif 


#endif
 }

}
#endif

//#####################################################################
// Copyright 2010, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __Vector_Subtract_Helper__
#define __Vector_Subtract_Helper__

#include "../Kernel_Base/Kernel_Base_Helper.h"
#include "../Kernel_Base/Kernel_Mesh_Base_Helper.h"
#include "../Kernel_Base/Kernel_Serial_Base_Helper.h"
#include "../Kernel_Base/Kernel_OpenMP_Base_Helper.h"

namespace MT_Streaming_Kernels{


    template <class T, int stride=1> 
class Vector_Subtract_Op
    {
    private:
        T* const x;
        const T* const y;
    public:
        explicit Vector_Subtract_Op(T* const x_input, const T* const y_input) : x(x_input), y(y_input) {}
        void Execute(int index)
        {
#pragma vector always
#pragma ivdep
            for(int i=0;i<stride;i++)
                x[index+i] -= y[index+i];
        }
    };

template <class T> 
    class Vector_Subtract_Op<T,1>
    {
    private:
        T* const x;
        const T* const y;
    public:
        explicit Vector_Subtract_Op(T* const x_input, const T* const y_input) : x(x_input), y(y_input) {}
        void Execute(int index)
        {
            x[index] -= y[index];
        }
    };


 template <class T> inline
     void Vector_Subtract_Helper( T* const x_input, const T* const y_input,
                               const int* block_offsets_input,const int number_of_blocks_input,
                               const int* partition_offsets_input, const int number_of_partitions_input,
                               const int x_stride_input, const int y_stride_input, const int block_size_input)
     {
         Vector_Subtract_Op<T> op(x_input, y_input);
         Kernel_Base_Helper<Vector_Subtract_Op<T> > helper(op, block_offsets_input,  number_of_blocks_input,
                                                        partition_offsets_input,  number_of_partitions_input,
                                                        x_stride_input, y_stride_input, block_size_input);
         helper.Run_Parallel();
     };
                               
 template <class T, class LAYOUT> inline
     void Vector_Subtract_Helper( T* const x_input, const T* const y_input,
                               const LAYOUT* block_layouts_input,const int number_of_blocks_input,
                               const int* partition_offsets_input, const int number_of_partitions_input,
                               const int x_stride_input, const int y_stride_input, const int block_size_input)
     {
         Vector_Subtract_Op<T> op(x_input, y_input);
         Kernel_Mesh_Base_Helper<Vector_Subtract_Op<T>,LAYOUT > helper(op, block_layouts_input,  number_of_blocks_input,
                                                        partition_offsets_input,  number_of_partitions_input,
                                                        x_stride_input, y_stride_input, block_size_input);
         helper.Run_Parallel();
     };

 template <class T> inline
     void Vector_Subtract_Helper( T* const x_input, const T* const y_input,
                                  const int number_of_entries_input,
                                  const int number_of_partitions_input)                              
     {
         Vector_Subtract_Op<T> op(x_input, y_input);
         Kernel_Serial_Base_Helper<Vector_Subtract_Op<T> > helper(op, number_of_entries_input,
                                                                number_of_partitions_input);
         helper.Run_Parallel();
     };

 template <class T, int stride> inline
     void Vector_Subtract_Helper(T* const x_input, const T* const y_input,
                            const int* offsets_input,
                            const int number_of_offsets_input,
                            const int number_of_partitions_input)
 {
     if(number_of_offsets_input == 0)
         return;

#if 0
     Vector_Subtract_Op<T,stride> op(x_input, y_input);
     Kernel_OpenMP_Base_Helper<Vector_Subtract_Op<T,stride> >helper(op, offsets_input, 
                                                                    number_of_offsets_input,
                                                                    number_of_partitions_input);
     helper.Run_Parallel();
#else


#ifdef __MIC__
     enum {width=stride};
     assert(width==16);
     
#pragma omp parallel for num_threads(number_of_partitions_input)
     for(int i=0;i<number_of_offsets_input;i++){
         _mm_prefetch((char*)&y_input[offsets_input[i+64]],_MM_HINT_T1);
         _mm_prefetch((char*)&y_input[offsets_input[i+8]],_MM_HINT_T0);
         __m512 rY=_mm512_setzero_ps();
         rY=_mm512_load_ps(&y_input[offsets_input[i]]);
         
         _mm_prefetch((char*)&x_input[offsets_input[i+64]],_MM_HINT_T1);
         _mm_prefetch((char*)&x_input[offsets_input[i+8]],_MM_HINT_T0);
         __m512 rX=_mm512_setzero_ps();
         rX=_mm512_load_ps(&x_input[offsets_input[i]]);        
         rX=_mm512_sub_ps(rX,rY);

         _mm512_storenrngo_ps(&x_input[offsets_input[i]],rX);
     }
#else
         
#pragma omp parallel for num_threads(number_of_partitions_input)
     for(int i=0;i<number_of_offsets_input;i++){
         int offset = offsets_input[i];
         for(int j=0;j<stride;j++)
             x_input[offset+j] -= y_input[offset+j];
    }

#endif 


#endif
 }
}

#endif



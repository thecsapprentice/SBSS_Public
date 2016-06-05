//#####################################################################
// Copyright 2010, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include "Vector_Set_Helper.h"
#include <algorithm>
#include <iostream>
#include <Thread_Queueing/PTHREAD_QUEUE.h>
using namespace MT_Streaming_Kernels;
// using namespace PhysBAM;
// extern PTHREAD_QUEUE* pthread_queue;
// //#####################################################################
// // Function Run_Parallel
// //#####################################################################
// namespace{
// template<class T>
// struct Vector_Set_Thread_Helper:public PTHREAD_QUEUE::TASK
// {
//     Vector_Set_Helper<T>* const obj;
//     const int index_start,index_end;
//     Vector_Set_Thread_Helper(Vector_Set_Helper<T>* const obj_input,const int index_start_input,const int index_end_input)
//         :obj(obj_input),index_start(index_start_input),index_end(index_end_input) {}
//     void Run(){obj->Run_Index_Range(index_start,index_end);}
// };
// }



// template<class T> void Vector_Set_Helper<T>::
// Run_Parallel()
// {
//     for(int partition=0;partition<number_of_partitions;partition++){
//         int block_begin=partition_offsets[partition];
//         int block_end=((partition<number_of_partitions-1)?partition_offsets[partition+1]:number_of_blocks);
        
//         std::cout << "Partition " << partition << " -- " << block_begin << " - " << block_end << std::endl;

//         Vector_Set_Thread_Helper<T>* task=
//             new Vector_Set_Thread_Helper<T>(this,block_begin,block_end);
//         pthread_queue->Queue(task);
//     }
//     pthread_queue->Wait();
// }
// //#####################################################################
// // Function Run_Index_Range
// //#####################################################################
// template<class T> void Vector_Set_Helper<T>::
// Run_Index_Range(const int block_start,const int block_end)
// {   
//     for(int block=block_start;block<block_end;block++)
//         for(int i=0;i<block_size;i++)
//             for(int j=0;j<block_size;j++)
//                 for(int k=0;k<block_size;k++){
//                     int index = block_offsets[block] + (i)*x_stride+(j)*y_stride+(k)*1; 
//                     x[index] *= a;
//                 }
// }
// //#####################################################################
// template class Vector_Set_Helper<float>;
// template class Vector_Set_Helper<double>;

//template class Kernel_Base_Helper<Vector_Set_Op<float> >;
//template class Kernel_Base_Helper<Vector_Set_Op<double> >;

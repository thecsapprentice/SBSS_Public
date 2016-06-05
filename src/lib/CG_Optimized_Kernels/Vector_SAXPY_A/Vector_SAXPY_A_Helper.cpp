//#####################################################################
// Copyright 2010, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include "Vector_SAXPY_A_Helper.h"
#include <algorithm>
#include <iostream>
#include <Thread_Queueing/PTHREAD_QUEUE.h>
using namespace MT_Streaming_Kernels;
/*
using namespace PhysBAM;
extern PTHREAD_QUEUE* pthread_queue;
//#####################################################################
// Function Run_Parallel
//#####################################################################
namespace{
template<class T>
struct Vector_SAXPY_A_Thread_Helper:public PTHREAD_QUEUE::TASK
{
    Vector_SAXPY_A_Helper<T>* const obj;
    const int index_start,index_end;
    Vector_SAXPY_A_Thread_Helper(Vector_SAXPY_A_Helper<T>* const obj_input,const int index_start_input,const int index_end_input)
        :obj(obj_input),index_start(index_start_input),index_end(index_end_input) {}
    void Run(){obj->Run_Index_Range(index_start,index_end);}
};
}



template<class T> void Vector_SAXPY_A_Helper<T>::
Run_Parallel()
{
    for(int partition=0;partition<number_of_partitions;partition++){
        int block_begin=partition_offsets[partition];
        int block_end=((partition<number_of_partitions-1)?partition_offsets[partition+1]:number_of_blocks);
        
        Vector_SAXPY_A_Thread_Helper<T>* task=
            new Vector_SAXPY_A_Thread_Helper<T>(this,block_begin,block_end);
        pthread_queue->Queue(task);
    }
    pthread_queue->Wait();
}
//#####################################################################
// Function Run_Index_Range
//#####################################################################
template<class T> void Vector_SAXPY_A_Helper<T>::
Run_Index_Range(const int index_start,const int index_end)
{   
    for(int index=index_start;index<=index_end;index++)
        {
            x[index] = c * y[index];
        }
}
//#####################################################################
template class Vector_SAXPY_A_Helper<float>;
template class Vector_SAXPY_A_Helper<double>;

*/

//template class Kernel_Base_Helper<Vector_SAXPY_A_Op<float> >;
//template class Kernel_Base_Helper<Vector_SAXPY_A_Op<double> >;

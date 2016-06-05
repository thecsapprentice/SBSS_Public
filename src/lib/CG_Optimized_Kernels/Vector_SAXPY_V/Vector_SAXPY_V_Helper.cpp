//#####################################################################
// Copyright 2010, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include "Vector_SAXPY_V_Helper.h"
#include <algorithm>
#include <iostream>
#include <Thread_Queueing/PTHREAD_QUEUE.h>
using namespace MT_Streaming_Kernels;
using namespace PhysBAM;
extern PTHREAD_QUEUE* pthread_queue;
//#####################################################################
// Function Run_Parallel
//#####################################################################
namespace{
template<class T, bool ACCUM>
struct Vector_SAXPY_V_Thread_Helper:public PTHREAD_QUEUE::TASK
{
    Vector_SAXPY_V_Helper<T,ACCUM>* const obj;
    const int index_start,index_end;
    Vector_SAXPY_V_Thread_Helper(Vector_SAXPY_V_Helper<T,ACCUM>* const obj_input,const int index_start_input,const int index_end_input)
        :obj(obj_input),index_start(index_start_input),index_end(index_end_input) {}
    void Run(){obj->Run_Index_Range(index_start,index_end);}
};
}
template<class T, bool ACCUM> void Vector_SAXPY_V_Helper<T, ACCUM>::
Run_Parallel(const int number_of_partitions)
{
    for(int partition=0;partition<number_of_partitions;partition++){
        int first_index_of_partition=(size/number_of_partitions)*(partition)+std::min(size%number_of_partitions,partition);
        int last_index_of_partition=(size/number_of_partitions)*(partition+1)+std::min(size%number_of_partitions,partition+1)-1;
        Vector_SAXPY_V_Thread_Helper<T,ACCUM>* task=new Vector_SAXPY_V_Thread_Helper<T,ACCUM>(this,first_index_of_partition,last_index_of_partition);
        pthread_queue->Queue(task);
    }
    pthread_queue->Wait();
}
//#####################################################################
// Function Run_Index_Range
//#####################################################################
template<class T, bool ACCUM> void Vector_SAXPY_V_Helper<T,ACCUM>::
Run_Index_Range(const int index_start,const int index_end)
{   
    if( ACCUM )
        for(int index=index_start;index<=index_end;index++)
            {
                x[index] += (v[index] * y[index]);
            }
    else
        for(int index=index_start;index<=index_end;index++)
            {
                x[index] = (v[index] * y[index]);
            }   
}
//#####################################################################
template class Vector_SAXPY_V_Helper<float,true>;
template class Vector_SAXPY_V_Helper<double,true>;

template class Vector_SAXPY_V_Helper<float,false>;
template class Vector_SAXPY_V_Helper<double,false>;


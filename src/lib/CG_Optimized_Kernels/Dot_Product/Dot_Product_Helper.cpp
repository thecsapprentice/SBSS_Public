//#####################################################################
// Copyright 2010, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include "Dot_Product_Helper.h"
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
struct Dot_Product_Thread_Helper:public PTHREAD_QUEUE::TASK
{
    Dot_Product_Helper<T>* const obj;
    const int index_start,index_end;
    double &result;
    Dot_Product_Thread_Helper(Dot_Product_Helper<T>* const obj_input,const int index_start_input,const int index_end_input,double &result_input)
        :obj(obj_input),index_start(index_start_input),index_end(index_end_input),result(result_input) {}
    void Run(){result=obj->Run_Index_Range(index_start,index_end);}
};
}
template<class T> double Dot_Product_Helper<T>::
Run_Parallel(const int number_of_partitions)
{
    double partial_results[number_of_partitions];
    for(int partition=0;partition<number_of_partitions;partition++){
        int first_index_of_partition=(size/number_of_partitions)*(partition)+std::min(size%number_of_partitions,partition);
        int last_index_of_partition=(size/number_of_partitions)*(partition+1)+std::min(size%number_of_partitions,partition+1)-1;
        Dot_Product_Thread_Helper<T>* task=new Dot_Product_Thread_Helper<T>(this,first_index_of_partition,last_index_of_partition,partial_results[partition]);
        pthread_queue->Queue(task);
    }
    pthread_queue->Wait();

    double result=0.;
    for(int partition=0;partition<number_of_partitions;partition++) result+=partial_results[partition];
    return result;
}
//#####################################################################
// Function Run_Index_Range
//#####################################################################
template<class T> double Dot_Product_Helper<T>::
Run_Index_Range(const int index_start,const int index_end)
{   
    double result=0.;
    for(int index=index_start;index<=index_end;index++)
        {
            result+=(double)x[index]*(double)y[index];
        }
    return result;
}
//#####################################################################
template class Dot_Product_Helper<float>;
template class Dot_Product_Helper<double>;

*/


//template class Kernel_Base_Reducer_Helper<Vector_Dot_Product_Op<float> >;
//template class Kernel_Base_Reducer_Helper<Vector_Dot_Product_Op<double> >;

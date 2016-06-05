//###########5##########################################################
// Copyright 2010, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include "Project_Helper.h"
#include <algorithm>
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
struct Project_Thread_Helper:public PTHREAD_QUEUE::TASK
{
    Project_Helper<T>* const obj;
    const int index_start,index_end;
    Project_Thread_Helper(Project_Helper<T>* const obj_input,const int index_start_input,const int index_end_input)
        :obj(obj_input),index_start(index_start_input),index_end(index_end_input){}
    void Run(){obj->Run_Index_Range(index_start,index_end);}
};
}
template<class T> void Project_Helper<T>::
Run_Parallel(const int number_of_partitions)
{
    for(int partition=0;partition<number_of_partitions;partition++){
	int first_index_of_partition=(size/number_of_partitions)*(partition)+std::min(size%number_of_partitions,partition);
	int last_index_of_partition=(size/number_of_partitions)*(partition+1)+std::min(size%number_of_partitions,partition+1)-1;
	Project_Thread_Helper<T>* task=new Project_Thread_Helper<T>(this,first_index_of_partition,last_index_of_partition);
	pthread_queue->Queue(task);}
    pthread_queue->Wait();
}
//#####################################################################
// Function Run_Index_Range
//#####################################################################
template<class T> void Project_Helper<T>::
Run_Index_Range(const int index_start,const int index_end)
{   
    for(int index=index_start;index<=index_end;index++)
        if(d[index]) x[index] = T();
}
//#####################################################################
template class Project_Helper<float>;
template class Project_Helper<double>;

*/

//template class Kernel_Base_Helper<Vector_Project_Op<float> >;
//template class Kernel_Base_Helper<Vector_Project_Op<double> >;

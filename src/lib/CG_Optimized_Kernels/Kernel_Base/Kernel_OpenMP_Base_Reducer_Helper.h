//#####################################################################
// Copyright 2010, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __Kernel_OpenMP_Base_Reducer_Helper__
#define __Kernel_OpenMP_Base_Reducer_Helper__

#include <Thread_Queueing/PTHREAD_QUEUE.h>
using namespace PhysBAM;
extern PTHREAD_QUEUE* pthread_queue;

namespace MT_Streaming_Kernels{

    template<class OP> class Kernel_OpenMP_Base_Reducer_Helper;

#if !defined(_OPENMP) 
    namespace{
        template<class OP>
            struct Kernel_OpenMP_Base_Reducer_Thread_Helper:public PTHREAD_QUEUE::TASK
            {
                Kernel_OpenMP_Base_Reducer_Helper<OP>* const obj;
                const int offset_start,offset_end;
                double &result;
            Kernel_OpenMP_Base_Reducer_Thread_Helper(Kernel_OpenMP_Base_Reducer_Helper<OP>* const obj_input,const int offset_start_input,const int offset_end_input, double &result_input)
                :obj(obj_input),offset_start(offset_start_input),offset_end(offset_end_input),result(result_input) {}
                void Run(){result=obj->Run_Index_Range(offset_start,offset_end);}
            };
    }
#endif    


    template<class OP>
        class Kernel_OpenMP_Base_Reducer_Helper
        {
            OP& op;
            const int *offsets;
            const int number_of_offsets;
            const int number_of_partitions;

        public:
            explicit Kernel_OpenMP_Base_Reducer_Helper(OP& op_input,
                                                       const int* offsets_input,
                                                       const int number_of_offsets_input,
                                                       const int number_of_partitions_input
                                        )
                :op(op_input),
                offsets(offsets_input),
                number_of_offsets(number_of_offsets_input),
                number_of_partitions(number_of_partitions_input)
                {}
    
            double Run()
            {
                return Run_Index_Range(0,number_of_offsets);
            }
  
//#####################################################################

            double Run_Parallel()
            { 
                if(!number_of_offsets)
                    return 0;
#if defined(_OPENMP)                               
                double partial_results[number_of_partitions];
                for(int i=0;i<number_of_partitions;i++)
                    partial_results[i] = 0.0;

                double result=0;
                int partition;
                const int chunk_size = number_of_offsets/number_of_partitions;
                
#pragma omp parallel for schedule(static,chunk_size) shared(partial_results) private(partition) num_threads(number_of_partitions)
                for(int i=0; i<number_of_offsets; i++){
                    partition = i/chunk_size;                    
                    op.Execute(offsets[i], partial_results[partition]);
                }
                
                for(int i=0;i<number_of_partitions;i++){
                    op.Join(result, partial_results[i]);
                }
                         
                return result;
#else
                int partition_offsets[number_of_partitions];
                int partition;
                partition_offsets[0]=0;
                for( int i=0, partition=0; i < number_of_offsets && partition < number_of_partitions;
                     i+=number_of_offsets / number_of_partitions, partition++ )
                    partition_offsets[partition] = i;

                double partial_results[number_of_partitions];

                for(partition=0;partition<number_of_partitions;partition++){
                    int offset_begin=partition_offsets[partition];
                    int offset_end=((partition<number_of_partitions-1)?partition_offsets[partition+1]:number_of_offsets);
            
                    Kernel_OpenMP_Base_Reducer_Thread_Helper<OP>* task=
                        new Kernel_OpenMP_Base_Reducer_Thread_Helper<OP>(this,offset_begin,offset_end, partial_results[partition]);
                    pthread_queue->Queue(task);
                }
                pthread_queue->Wait();
                
                double result=0.;
                for(partition=0;partition<number_of_partitions;partition++)
                    op.Join(result,partial_results[partition]);

                return result;
#endif
            }

//#####################################################################


            double Run_Index_Range(const int offset_start,const int offset_end)
            {
                double result=0.;
                for(int offset=offset_start;offset<offset_end;offset++)
                    op.Execute( offset, result ); 
                return result;
            }


//#####################################################################
        };
}
#endif

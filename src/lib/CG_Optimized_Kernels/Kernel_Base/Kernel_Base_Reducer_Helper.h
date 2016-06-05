//#####################################################################
// Copyright 2010, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __Kernel_Base_Reducer_Helper__
#define __Kernel_Base_Reducer_Helper__

#include <Thread_Queueing/PTHREAD_QUEUE.h>
using namespace PhysBAM;
extern PTHREAD_QUEUE* pthread_queue;

namespace MT_Streaming_Kernels{

    template<class OP> class Kernel_Base_Reducer_Helper;

    namespace{
        template<class OP>
            struct Kernel_Base_Reducer_Thread_Helper:public PTHREAD_QUEUE::TASK
            {
                Kernel_Base_Reducer_Helper<OP>* const obj;
                const int block_start,block_end;
                double &result;
            Kernel_Base_Reducer_Thread_Helper(Kernel_Base_Reducer_Helper<OP>* const obj_input,const int block_start_input,const int block_end_input, double &result_input)
                :obj(obj_input),block_start(block_start_input),block_end(block_end_input),result(result_input) {}
                void Run(){result=obj->Run_Index_Range(block_start,block_end);}
            };
    }
    

    template<class OP>
        class Kernel_Base_Reducer_Helper
        {
            OP& op;
            const int* block_offsets;
            const int number_of_blocks;
            const int* partition_offsets;
            const int number_of_partitions;
            const int x_stride;
            const int y_stride;
            const int block_size;


        public:
            explicit Kernel_Base_Reducer_Helper(OP& op_input,
                                        const int* block_offsets_input,const int number_of_blocks_input,
                                        const int* partition_offsets_input, const int number_of_partitions_input,
                                        const int x_stride_input, const int y_stride_input, const int block_size_input)
                :op(op_input),
                block_offsets(block_offsets_input),
                number_of_blocks(number_of_blocks_input),
                partition_offsets(partition_offsets_input),
                number_of_partitions(number_of_partitions_input),
                x_stride(x_stride_input),
                y_stride(y_stride_input),
                block_size(block_size_input)
                {}
    
            double Run()
            {
                return Run_Index_Range(0,number_of_blocks);
            }
  
//#####################################################################

            double Run_Parallel()
            { 
                double partial_results[number_of_partitions];

                for(int partition=0;partition<number_of_partitions;partition++){
                    int block_begin=partition_offsets[partition];
                    int block_end=((partition<number_of_partitions-1)?partition_offsets[partition+1]:number_of_blocks);
            
                    Kernel_Base_Reducer_Thread_Helper<OP>* task=
                        new Kernel_Base_Reducer_Thread_Helper<OP>(this,block_begin,block_end, partial_results[partition]);
                    pthread_queue->Queue(task);
                }
                pthread_queue->Wait();
                
                double result=0.;
                for(int partition=0;partition<number_of_partitions;partition++)
                    op.Join(result,partial_results[partition]);
                return result;
            }

//#####################################################################


            double Run_Index_Range(const int block_start,const int block_end)
            {
                double result=0.;
                for(int block=block_start;block<block_end;block++)
                    for(int i=0;i<block_size;i++)
                        for(int j=0;j<block_size;j++)
                            for(int k=0;k<block_size;k++){
                                int index = block_offsets[block] + (i)*x_stride+(j)*y_stride+(k)*1; 
                                op.Execute(index, result);
                            }
                return result;
            }


//#####################################################################
        };
}
#endif

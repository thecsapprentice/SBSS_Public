//#####################################################################
// Copyright 2010, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __Kernel_Mesh_Base_Helper__
#define __Kernel_Mesh_Base_Helper__

#include <Thread_Queueing/PTHREAD_QUEUE.h>

using namespace PhysBAM;
extern PTHREAD_QUEUE* pthread_queue;

namespace MT_Streaming_Kernels{

    template<class OP, class LAYOUT> class Kernel_Mesh_Base_Helper;

    namespace{
        template<class OP, class LAYOUT>
            struct Kernel_Mesh_Base_Thread_Helper:public PTHREAD_QUEUE::TASK
            {
                Kernel_Mesh_Base_Helper<OP, LAYOUT>* const obj;
                const int block_start,block_end;
            Kernel_Mesh_Base_Thread_Helper(Kernel_Mesh_Base_Helper<OP, LAYOUT>* const obj_input,const int block_start_input,const int block_end_input)
        :obj(obj_input),block_start(block_start_input),block_end(block_end_input) {}
                void Run(){obj->Run_Index_Range(block_start,block_end);}
            };
    }
    

    template<class OP, class LAYOUT>
        class Kernel_Mesh_Base_Helper
        {
            OP& op;
            const LAYOUT* block_layouts;
            const int number_of_blocks;
            const int* partition_offsets;
            const int number_of_partitions;
            const int x_stride;
            const int y_stride;
            const int block_size;


        public:
            explicit Kernel_Mesh_Base_Helper(OP& op_input,
                                        const LAYOUT* block_layouts_input,const int number_of_blocks_input,
                                        const int* partition_offsets_input, const int number_of_partitions_input,
                                        const int x_stride_input, const int y_stride_input, const int block_size_input)
                :op(op_input),
                block_layouts(block_layouts_input),
                number_of_blocks(number_of_blocks_input),
                partition_offsets(partition_offsets_input),
                number_of_partitions(number_of_partitions_input),
                x_stride(x_stride_input),
                y_stride(y_stride_input),
                block_size(block_size_input)
                {}
    
            void Run()
            {
                Run_Index_Range(0,number_of_blocks);
            }
  
//#####################################################################

            void Run_Parallel()
            { 
                for(int partition=0;partition<number_of_partitions;partition++){
                    int block_begin=partition_offsets[partition];
                    int block_end=((partition<number_of_partitions-1)?partition_offsets[partition+1]:number_of_blocks);
            
                    Kernel_Mesh_Base_Thread_Helper<OP, LAYOUT>* task=
                        new Kernel_Mesh_Base_Thread_Helper<OP, LAYOUT>(this,block_begin,block_end);
                    pthread_queue->Queue(task);
                }
                pthread_queue->Wait();
            }

//#####################################################################


            void Run_Index_Range(const int block_start,const int block_end)
            {
                for(int block=block_start;block<block_end;block++){
                    if( block_layouts[block].type == LAYOUT::GRID ){
                        assert( 0 == 1 );
                    }
                    else{
                        const typename LAYOUT::LAYOUT_MIX& block_mix = block_layouts[block].block_layout;

                        for( int index=0; index<block_mix.SIZE; index++){
                            if( block_mix[index] ) // Skip all grid nodes
                                continue;
                            int mesh_index = block_mix.cell[index].mesh_index - 1 ; // This is the index into the mesh element
                            op.Execute( mesh_index ); 
                        }
                    }
                }
            }


//#####################################################################
        };
}
#endif

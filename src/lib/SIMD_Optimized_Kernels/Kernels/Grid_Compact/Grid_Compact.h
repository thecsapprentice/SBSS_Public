//#####################################################################
//  Copyright (c) 2011-2013 Nathan Mitchell, Eftychios Sifakis.
//  This file is covered by the FreeBSD license. Please refer to the 
//  license.txt file for more information.
//#####################################################################



template<class T>
class Grid_Compact_Task : public PhysBAM::PTHREAD_QUEUE::TASK
{
 private:
    T (*_u_compact)[3][27];
    T (*_p_compact)[8];

    const T* _uX_ptr;
    const T* _uY_ptr;
    const T* _uZ_ptr;
    const T* _p_ptr;

    const int  _first_block;
    const int* _block_offsets;
    const int  _number_of_blocks;
    const int  _x_stride_cell;
    const int  _y_stride_cell;
    const int  _x_stride_node;
    const int  _y_stride_node;

    const bool _compact_p;

 public:
    Grid_Compact_Task(T (*u_compact)[3][27], T (*p_compact)[8], const T * uX_ptr,
                      const T* uY_ptr, const T* uZ_ptr, const T* p_ptr,
                      const int first_block, const int* block_offsets,const int number_of_blocks,
                      const int x_stride_cell, const int y_stride_cell,
                      const int x_stride_node, const int y_stride_node,
                      const bool compact_p) : 
    _u_compact(u_compact), _p_compact(p_compact), 
        _uX_ptr(uX_ptr), _uY_ptr(uY_ptr), _uZ_ptr(uZ_ptr), _p_ptr(p_ptr),
        _first_block(first_block),
        _block_offsets(block_offsets),
        _number_of_blocks(number_of_blocks),
        _x_stride_cell(x_stride_cell),
        _y_stride_cell(y_stride_cell),
        _x_stride_node(x_stride_node),
        _y_stride_node(y_stride_node),
        _compact_p(compact_p)
        {
        }

    void Run()
    {
        for(int block=0;block<_number_of_blocks;block++)
            {
                int linear_index=0;
                for(int i=0;i<3;i++)
                    for(int j=0;j<3;j++)
                        for(int k=0;k<3;k++){
                            
                            int index = _block_offsets[block] + (i)*_x_stride_node+(j)*_y_stride_node+(k)*1; 

                            _u_compact[_first_block+block][0][linear_index]=*(_uX_ptr+index);
                            _u_compact[_first_block+block][1][linear_index]=*(_uY_ptr+index);
                            _u_compact[_first_block+block][2][linear_index]=*(_uZ_ptr+index);

                            linear_index++;
                        }

                if(_compact_p)
                    for(int i=0;i<2;i++)
                        for(int j=0;j<2;j++)
                            for(int k=0;k<2;k++)
                                _p_compact[_first_block+block][4*i+2*j+k]=(_p_ptr+(_block_offsets[block]))[i*_x_stride_cell+j*_y_stride_cell+k];
                
            }
    }
    
};

template<class T>
void Grid_Compact( T (*u_compact)[3][27], T (*p_compact)[8], const T * uX_ptr,
                   const T* uY_ptr, const T* uZ_ptr, const T* p_ptr,
                   const int* block_offsets,const int number_of_blocks,
                   const int* partition_offsets, const int number_of_partitions,
                   const int x_stride_cell, const int y_stride_cell,
                   const int x_stride_node, const int y_stride_node,
                   const bool compact_p = true)
{
	for(int partition=0;partition<number_of_partitions;partition++){
        int block_begin=partition_offsets[partition];
        int block_end=((partition<number_of_partitions-1)?partition_offsets[partition+1]:number_of_blocks);
        //std::cout << "Compacting data from block "<< block_begin << " to block " << block_end << std::endl;
        //std::cout << "      Number of blocks "<< block_end-block_begin << std::endl;

        Grid_Compact_Task<T>* task = new Grid_Compact_Task<T>(u_compact,
                                                              p_compact,
                                                              uX_ptr,
                                                              uY_ptr,
                                                              uZ_ptr,
                                                              p_ptr,
                                                              block_begin,
                                                              block_offsets+block_begin,
                                                              block_end-block_begin,
                                                              x_stride_cell,
                                                              y_stride_cell,
                                                              x_stride_node,
                                                              y_stride_node,
                                                              compact_p
                                                              );
        pthread_queue->Queue(task);
	}
	pthread_queue->Wait();
}

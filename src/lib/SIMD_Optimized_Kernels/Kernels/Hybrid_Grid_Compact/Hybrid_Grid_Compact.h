
#ifdef ENABLE_MIC_INSTRUCTION_SET
#include <immintrin.h>
#include <zmmintrin.h>
#endif

template<class T>
class Hybrid_Grid_Compact_Task : public PhysBAM::PTHREAD_QUEUE::TASK
{
 private:
    T (*_u_compact)[3][27];
    T (*_p_compact)[8];

    const T* _uX_ptr;
    const T* _uY_ptr;
    const T* _uZ_ptr;
    const T* _p_ptr;

    const T* _uX_mesh_ptr;
    const T* _uY_mesh_ptr;
    const T* _uZ_mesh_ptr;
    const T* _p_mesh_ptr;

    const int  _first_block;
    const int* _block_offsets;
    const int* _block_cells;
    const int* _block_nodes;
    const int  _number_of_blocks;
    const int  _x_stride_cell;
    const int  _y_stride_cell;
    const int  _x_stride_node;
    const int  _y_stride_node;

    const bool _compact_p;

 public:
    Hybrid_Grid_Compact_Task(T (*u_compact)[3][27], T (*p_compact)[8],
                    const T * uX_ptr, const T* uY_ptr, const T* uZ_ptr, const T* p_ptr,
                    const T * uX_mesh_ptr, const T* uY_mesh_ptr, const T* uZ_mesh_ptr, const T* p_mesh_ptr,
                    const int first_block, const int* block_offsets, const int* block_cells,
                    const int* block_nodes, const int number_of_blocks,
                    const int x_stride_cell, const int y_stride_cell,
                    const int x_stride_node, const int y_stride_node,
                    const bool compact_p         
                             ) : 
    _u_compact(u_compact), _p_compact(p_compact), 
        _uX_ptr(uX_ptr), _uY_ptr(uY_ptr), _uZ_ptr(uZ_ptr), _p_ptr(p_ptr),
        _uX_mesh_ptr(uX_mesh_ptr), _uY_mesh_ptr(uY_mesh_ptr), _uZ_mesh_ptr(uZ_mesh_ptr), _p_mesh_ptr(p_mesh_ptr),
        _first_block(first_block),
        _block_offsets(block_offsets),
        _block_cells(block_cells),
        _block_nodes(block_nodes),
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
                            
                            int base_index = _block_offsets[block*3+0] * _x_stride_node + 
                                _block_offsets[block*3+1] * _y_stride_node + 
                                _block_offsets[block*3+2]; 
                            int node_index = base_index + (i)*_x_stride_node+(j)*_y_stride_node+(k)*1; 
                            
                            switch( _block_nodes[block*27 + linear_index] ){
                            case -1:
                                // Do nothing for this location.
                                _u_compact[_first_block+block][0][linear_index]=0;
                                _u_compact[_first_block+block][1][linear_index]=0;
                                _u_compact[_first_block+block][2][linear_index]=0;                        
                                break;
                            case 0:
                                // Copy from grid for this location
                                _u_compact[_first_block+block][0][linear_index]=*(_uX_ptr+node_index);
                                _u_compact[_first_block+block][1][linear_index]=*(_uY_ptr+node_index);
                                _u_compact[_first_block+block][2][linear_index]=*(_uZ_ptr+node_index);
                                break;
                            default:
                                // Any other value must be from the mesh
                                int mesh_index = _block_nodes[block*27 + linear_index]-1;
                                _u_compact[_first_block+block][0][linear_index]=*(_uX_mesh_ptr+mesh_index);
                                _u_compact[_first_block+block][1][linear_index]=*(_uY_mesh_ptr+mesh_index);
                                _u_compact[_first_block+block][2][linear_index]=*(_uZ_mesh_ptr+mesh_index);
                                break;
                            }
                            linear_index++;
                        }
                
                if(_compact_p){
                linear_index=0;
                for(int i=0;i<2;i++)
                    for(int j=0;j<2;j++)
                        for(int k=0;k<2;k++){
                            int base_index = _block_offsets[block*3+0] * _x_stride_cell + 
                                _block_offsets[block*3+1] * _y_stride_cell + 
                                _block_offsets[block*3+2]; 
                            int cell_index = base_index + (i)*_x_stride_cell+(j)*_y_stride_cell+(k)*1; 
                            
                            switch( _block_cells[block*8 + linear_index] ){
                            case -1:
                                // Do nothing for this location.
                                _p_compact[_first_block+block][linear_index]=0;
                                break;
                            case 0:
                                // Copy from grid for this location
                                _p_compact[_first_block+block][linear_index]=*(_p_ptr+cell_index);
                                break;
                            default:
                                // Any other value must be from the mesh
                                int mesh_index = _block_cells[block*8 + linear_index]-1;
                                _p_compact[_first_block+block][linear_index]=*(_p_mesh_ptr+mesh_index);;
                                break;
                            }
                            linear_index++;
                        }
                }
            }
    }
};

template<class T>
void Hybrid_Grid_Compact( T (*u_compact)[3][27], T (*p_compact)[8],
                          const T * uX_ptr, const T* uY_ptr, const T* uZ_ptr, const T* p_ptr,
                          const T * uX_mesh_ptr, const T* uY_mesh_ptr, const T* uZ_mesh_ptr, const T* p_mesh_ptr,
                          const int number_of_blocks, 
                          const int* block_offsets,
                          const int* block_cells,
                          const int* block_nodes,
                          const int number_of_partitions,
                          const int x_stride_cell, const int y_stride_cell,
                          const int x_stride_node, const int y_stride_node,
                          const bool compact_p = true)
{
    Hybrid_Grid_Compact_Task<T>** tasks = new Hybrid_Grid_Compact_Task<T>*[number_of_partitions];

	for(int partition=0;partition<number_of_partitions;partition++){
        int block_begin=(partition)*(number_of_blocks/number_of_partitions)+min(partition,number_of_blocks%number_of_partitions);
        int block_end=(partition+1)*(number_of_blocks/number_of_partitions)+min((partition+1),number_of_blocks%number_of_partitions);
        //std::cout << "Compacting data from block "<< block_begin << " to block " << block_end << std::endl;
        //std::cout << "      Number of blocks "<< block_end-block_begin << std::endl;

        tasks[partition] = new Hybrid_Grid_Compact_Task<T>(u_compact,
                                                              p_compact,
                                                              uX_ptr,
                                                              uY_ptr,
                                                              uZ_ptr,
                                                              p_ptr,
                                                              uX_mesh_ptr,
                                                              uY_mesh_ptr,
                                                              uZ_mesh_ptr,
                                                              p_mesh_ptr,
                                                              block_begin,
                                                              block_offsets+(block_begin*3),
                                                              block_cells+(block_begin*8),
                                                              block_nodes+(block_begin*27),
                                                              block_end-block_begin,
                                                              x_stride_cell,
                                                              y_stride_cell,
                                                              x_stride_node,
                                                              y_stride_node,
                                                              compact_p
                                                              );
	}

#pragma omp parallel for num_threads(number_of_partitions)
    for(int partition=0;partition<number_of_partitions;partition++){
        tasks[partition]->Run();
        delete tasks[partition];
    }
    delete[] tasks;
}

template<class T>
void Hybrid_Grid_Compact( T (*u_compact)[3][27], T (*p_compact)[8],
                          const T * uX_ptr, const T* uY_ptr, const T* uZ_ptr, const T* p_ptr,
                          const T * uX_mesh_ptr, const T* uY_mesh_ptr, const T* uZ_mesh_ptr, const T* p_mesh_ptr,
                          const int number_of_blocks, 
                          const int* mic_block_offsets_low,
                          const int* mic_block_offsets_high,
                          const int* mic_block_offsets_grid_mask,
                          const int* mic_block_offsets_mesh_mask,
                          const int* block_offsets,
                          const int* block_cells,
                          const int* block_nodes,
                          const int number_of_partitions,
                          const int x_stride_cell, const int y_stride_cell,
                          const int x_stride_node, const int y_stride_node,
                          const bool compact_p = true)
{
#ifdef ENABLE_MIC_INSTRUCTION_SET
#pragma omp parallel for
    for( int block=0; block < number_of_blocks; block++){
        int high_bit_mask = 0xFFFF0000;
        int low_bit_mask = 0x0000FFFF;
        __mmask16 low_grid_mask=_mm512_int2mask( (mic_block_offsets_grid_mask[block] & low_bit_mask));
        __mmask16 high_grid_mask=_mm512_int2mask( ((mic_block_offsets_grid_mask[block] & high_bit_mask) >> 16) );
        __mmask16 low_mesh_mask=_mm512_int2mask( (mic_block_offsets_mesh_mask[block] & low_bit_mask));
        __mmask16 high_mesh_mask=_mm512_int2mask( ((mic_block_offsets_mesh_mask[block] & high_bit_mask) >> 16) );
        __mmask16 low_mask=_mm512_int2mask(  0x0000FFFF );
        __mmask16 high_mask=_mm512_int2mask( 0x000007FF );
        __m512i offsets_low;
        __m512i offsets_high;
        __m512 data_low;
        __m512 data_high;

        offsets_low = _mm512_load_epi32(&mic_block_offsets_low[block*16]);
        offsets_high = _mm512_load_epi32(&mic_block_offsets_high[block*16]);       

        const T* grid_ptrs[3]; grid_ptrs[0] = uX_ptr; grid_ptrs[1] = uY_ptr; grid_ptrs[2] = uZ_ptr;
        const T* mesh_ptrs[3]; mesh_ptrs[0] = uX_mesh_ptr; mesh_ptrs[1] = uY_mesh_ptr; mesh_ptrs[2] = uZ_mesh_ptr;
        
        for( int v = 0; v < 3; v++)
        {             
            data_low = _mm512_setzero_ps();
            data_high = _mm512_setzero_ps();

            data_low =  _mm512_mask_i32gather_ps(data_low, low_grid_mask, offsets_low, grid_ptrs[v], 4);
            data_high =  _mm512_mask_i32gather_ps(data_high, high_grid_mask, offsets_high, grid_ptrs[v], 4);
            data_low =  _mm512_mask_i32gather_ps(data_low, low_mesh_mask, offsets_low, mesh_ptrs[v], 4);
            data_high =  _mm512_mask_i32gather_ps(data_high, high_mesh_mask, offsets_high, mesh_ptrs[v], 4);

            _mm512_mask_packstorelo_ps(&(u_compact[block][v][0]),low_mask,data_low);
            _mm512_mask_packstorehi_ps(((void*)&(u_compact[block][v][0]))+64,low_mask,data_low);
            _mm512_mask_packstorelo_ps(&(u_compact[block][v][16]),high_mask,data_high);
            _mm512_mask_packstorehi_ps(((void*)&(u_compact[block][v][16]))+64,high_mask,data_high);
        }
        
        if(compact_p){
            int linear_index=0;
            for(int i=0;i<2;i++)
                for(int j=0;j<2;j++)
                    for(int k=0;k<2;k++){
                        int base_index = block_offsets[block*3+0] * x_stride_cell + 
                            block_offsets[block*3+1] * y_stride_cell + 
                            block_offsets[block*3+2]; 
                        int cell_index = base_index + (i)*x_stride_cell+(j)*y_stride_cell+(k)*1; 
                            
                        switch( block_cells[block*8 + linear_index] ){
                        case -1:
                            // Do nothing for this location.
                            p_compact[block][linear_index]=0;
                            break;
                        case 0:
                            // Copy from grid for this location
                            p_compact[block][linear_index]=*(p_ptr+cell_index);
                            break;
                        default:
                            // Any other value must be from the mesh
                            int mesh_index = block_cells[block*8 + linear_index]-1;
                            p_compact[block][linear_index]=*(p_mesh_ptr+mesh_index);;
                            break;
                        }
                        linear_index++;
                    }
        }
    }    
#endif
}

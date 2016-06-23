//#####################################################################
// Copyright 2011, Taylor Patterson, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class NONLINEAR_ELASTICITY
//#####################################################################
#include <PhysBAM_Tools/Data_Structures/HASHTABLE_ITERATOR.h>

#include "SKINNING_NONLINEAR_ELASTICITY.h"
//#include "NONLINEAR_ELASTICITY.h"
#include <Common_Tools/Math_Tools/RANGE_ITERATOR.h>
#include "SPECIALIZED_KERNELS.h"

#include <Thread_Queueing/PTHREAD_QUEUE.h>
using namespace PhysBAM;
extern PTHREAD_QUEUE* pthread_queue;

#include <SIMD_Optimized_Kernels/Kernels/Grid_Compact/Grid_Compact.h>
#include <SIMD_Optimized_Kernels/Kernels/Grid_UnCompact/Grid_UnCompact.h>

#ifdef LOG_DETAILED_PERFORMANCE
#define LOG_DETAILED_PERFORMANCE_NE
#endif
template<class T,int d>
VECTOR<typename NONLINEAR_ELASTICITY<T,d>::T_SCALAR_VARIABLE_VIEW_CONST,d> NONLINEAR_ELASTICITY<T,d>::
View_Convert(const VECTOR<T_SCALAR_VARIABLE_VIEW,d>& in)
{
    return VECTOR<T_SCALAR_VARIABLE_VIEW_CONST,d>(in.x,in.y,in.z);
};

//#####################################################################
// Function Update_Position_Based_State
//#####################################################################
template<class T,int d> template<bool enable_constraints,bool enable_muscles> void NONLINEAR_ELASTICITY<T,d>::
Update_Position_Based_State_Specialized(T_VECTOR_VARIABLE_VIEW_CONST u,T_SCALAR_VARIABLE_VIEW_CONST p, T_VECTOR_VARIABLE_VIEW diag)
{
#ifdef LOG_DETAILED_PERFORMANCE_NE
    LOG::SCOPE scope("NONLINEAR_ELASTICITY::Update_Position_Based_State_Specialized");
#endif
    
#ifdef USE_SPECIALIZED_KERNELS

    CompactData_Specialized(u,p);
    typedef T (*u_compact_type)[3][27];
    typedef T (*p_compact_type)[8];
    u_compact_type u_compact=reinterpret_cast<u_compact_type>(specialized_data.u_compact_array.Get_Array_Pointer());
    p_compact_type p_compact=reinterpret_cast<p_compact_type>(specialized_data.p_compact_array.Get_Array_Pointer());
    u_compact_type d_compact=reinterpret_cast<u_compact_type>(specialized_data.d_compact_array.Get_Array_Pointer());

    Update_Position_Based_State_With_Specialized_Kernels<T,enable_constraints,enable_muscles>(
        u_compact,p_compact,d_compact,
        cell_block_base_offsets.m,
        cell_block_partition_offsets.m,
        muscle_fiber_max_stresses,muscle_activations,
        specialized_data);

    UnCompactData_Specialized(diag);

    // TODO: Figure out something productive to do with this
#if 0
    for(int block=1;block<=cell_block_base_indices.m;block++){
        const T_INDEX& base_index=cell_block_base_indices(block);

        int cell=0;
        for(RANGE_ITERATOR<d> cell_iterator(RANGE<T_INDEX>(base_index,base_index+1));cell_iterator.Valid();cell_iterator.Next(),cell++){
            const T_INDEX& cell_index=cell_iterator.Index();
            U_flat(block).Get(U(cell_index),cell);
            V_flat(block).Get(V(cell_index),cell);
            Sigma_flat(block).Get(Sigma(cell_index),cell);
            Q_hat_flat(block).Get(Q_hat(cell_index),cell);
            dPdF_flat(block).Get(dP_dF(cell_index),cell);}

        int muscle_block_offset = muscle_base_offsets(block);

        for( int layer=0; layer < muscle_base_lengths(block); layer++)
            {
                int muscle_block = muscle_block_offset + layer;
                cell=0;
                for(RANGE_ITERATOR<d> cell_iterator(RANGE<T_INDEX>(base_index,base_index+1));cell_iterator.Valid();cell_iterator.Next(),cell++){
                    const T_INDEX& cell_index=cell_iterator.Index();
                          
                    if(cell_F_fibers(cell_index).Size() > layer )
                        {
                            muscle_F_fiber(muscle_block).Get(cell_F_fibers(cell_index)(layer+1),cell);
                            muscle_c1(muscle_block).Get(cell_c1(cell_index)(layer+1),cell);
                            muscle_c2(muscle_block).Get(cell_c2(cell_index)(layer+1),cell);
                        }
                }
            }
    }
#endif
#endif
}
//#####################################################################
// Function CompactData_Specialized
//#####################################################################
template<class T,int d> void NONLINEAR_ELASTICITY<T,d>::
CompactData_Specialized(T_VECTOR_VARIABLE_VIEW_CONST x, T_SCALAR_VARIABLE_VIEW_CONST p) const 
{
#ifdef USE_SPECIALIZED_KERNELS

    VECTOR<int,3> cell_grid_dims = padded_cell_domain.max_corner - padded_cell_domain.min_corner +1;
    VECTOR<int,3> node_grid_dims = padded_node_domain.max_corner - padded_node_domain.min_corner +1 ;

    int X_StrideCell = cell_grid_dims(2)*cell_grid_dims(3);
    int Y_StrideCell =                   cell_grid_dims(3);

    int X_StrideNode = node_grid_dims(2)*node_grid_dims(3);
    int Y_StrideNode =                   node_grid_dims(3);
    
    typedef T (*u_compact_type)[3][27];
    typedef T (*p_compact_type)[8];
    u_compact_type u_or_f_compact=reinterpret_cast<u_compact_type>(specialized_data.u_compact_array.Get_Array_Pointer());
    p_compact_type p_or_q_compact=reinterpret_cast<p_compact_type>(specialized_data.p_compact_array.Get_Array_Pointer());

    Grid_Compact(u_or_f_compact,
                 p_or_q_compact,
                 &x(1)(0,0,0),
                 &x(2)(0,0,0),
                 &x(3)(0,0,0),
                 &p(0,0,0),
                 cell_block_base_offsets.Get_Array_Pointer(),cell_block_base_offsets.m,
                 cell_block_partition_offsets.Get_Array_Pointer(),cell_block_partition_offsets.m,
                 X_StrideCell, Y_StrideCell,
                 X_StrideNode, Y_StrideNode);
#endif
}
//#####################################################################
// Function UnCompactData_Specialized
//#####################################################################
template<class T,int d> void NONLINEAR_ELASTICITY<T,d>::
UnCompactData_Specialized(T_VECTOR_VARIABLE_VIEW x,T_SCALAR_VARIABLE_VIEW p) const 
{
#ifdef USE_SPECIALIZED_KERNELS

    VECTOR<int,3> cell_grid_dims = padded_cell_domain.max_corner - padded_cell_domain.min_corner +1;
    VECTOR<int,3> node_grid_dims = padded_node_domain.max_corner - padded_node_domain.min_corner +1 ;

    int X_StrideCell = cell_grid_dims(2)*cell_grid_dims(3);
    int Y_StrideCell =                   cell_grid_dims(3);

    int X_StrideNode = node_grid_dims(2)*node_grid_dims(3);
    int Y_StrideNode =                   node_grid_dims(3);
    
    typedef T (*u_compact_type)[3][27];
    typedef T (*p_compact_type)[8];
    u_compact_type u_or_f_compact=reinterpret_cast<u_compact_type>(specialized_data.u_compact_array.Get_Array_Pointer());
    p_compact_type p_or_q_compact=reinterpret_cast<p_compact_type>(specialized_data.p_compact_array.Get_Array_Pointer());

    Grid_UnCompact(u_or_f_compact,
                   p_or_q_compact,
                   &x(1)(0,0,0),
                   &x(2)(0,0,0),
                   &x(3)(0,0,0),
                   &p(0,0,0),
                   cell_block_base_offsets.Get_Array_Pointer(),cell_block_base_offsets.m,
                   cell_block_partition_offsets.Get_Array_Pointer(),cell_block_partition_offsets.m,
                   X_StrideCell, Y_StrideCell,
                   X_StrideNode, Y_StrideNode);
#endif
}
//#####################################################################
// Function CompactData_Specialized
//#####################################################################
template<class T,int d> void NONLINEAR_ELASTICITY<T,d>::
CompactData_Specialized(T_VECTOR_VARIABLE_VIEW_CONST x) const 
{
#ifdef USE_SPECIALIZED_KERNELS

    VECTOR<int,3> cell_grid_dims = padded_cell_domain.max_corner - padded_cell_domain.min_corner +1;
    VECTOR<int,3> node_grid_dims = padded_node_domain.max_corner - padded_node_domain.min_corner +1 ;

    int X_StrideCell = cell_grid_dims(2)*cell_grid_dims(3);
    int Y_StrideCell =                   cell_grid_dims(3);

    int X_StrideNode = node_grid_dims(2)*node_grid_dims(3);
    int Y_StrideNode =                   node_grid_dims(3);
    
    typedef T (*u_compact_type)[3][27];
    typedef T (*p_compact_type)[8];
    u_compact_type d_compact=reinterpret_cast<u_compact_type>(specialized_data.d_compact_array.Get_Array_Pointer());
    p_compact_type null_compact=reinterpret_cast<p_compact_type>(NULL);

    Grid_Compact(d_compact,
                 null_compact,
                 &x(1)(0,0,0),
                 &x(2)(0,0,0),
                 &x(3)(0,0,0),
                 (T*)(NULL),
                 cell_block_base_offsets.Get_Array_Pointer(),cell_block_base_offsets.m,
                 cell_block_partition_offsets.Get_Array_Pointer(),cell_block_partition_offsets.m,
                 X_StrideCell, Y_StrideCell,
                 X_StrideNode, Y_StrideNode,
                 false);
#endif
}
//#####################################################################
// Function UnCompactData_Specialized
//#####################################################################
template<class T,int d> void NONLINEAR_ELASTICITY<T,d>::
UnCompactData_Specialized(T_VECTOR_VARIABLE_VIEW x) const 
{
#ifdef USE_SPECIALIZED_KERNELS

    VECTOR<int,3> cell_grid_dims = padded_cell_domain.max_corner - padded_cell_domain.min_corner +1;
    VECTOR<int,3> node_grid_dims = padded_node_domain.max_corner - padded_node_domain.min_corner +1 ;

    int X_StrideCell = cell_grid_dims(2)*cell_grid_dims(3);
    int Y_StrideCell =                   cell_grid_dims(3);

    int X_StrideNode = node_grid_dims(2)*node_grid_dims(3);
    int Y_StrideNode =                   node_grid_dims(3);
    
    typedef T (*u_compact_type)[3][27];
    typedef T (*p_compact_type)[8];
    u_compact_type d_compact=reinterpret_cast<u_compact_type>(specialized_data.d_compact_array.Get_Array_Pointer());
    p_compact_type null_compact=reinterpret_cast<p_compact_type>(NULL);

    Grid_UnCompact(d_compact,
                   null_compact,
                   &x(1)(0,0,0),
                   &x(2)(0,0,0),
                   &x(3)(0,0,0),
                   (T*)(NULL),
                   cell_block_base_offsets.Get_Array_Pointer(),cell_block_base_offsets.m,
                   cell_block_partition_offsets.Get_Array_Pointer(),cell_block_partition_offsets.m,
                   X_StrideCell, Y_StrideCell,
                   X_StrideNode, Y_StrideNode,
                   false);
#endif
}
//#####################################################################
// Function Update_Position_Based_State
//#####################################################################
template<class T,int d>  template<bool enable_constraints,bool enable_muscles> void NONLINEAR_ELASTICITY<T,d>::
Add_Force_First_Order_Elasticity_Specialized(T_VECTOR_VARIABLE_VIEW_CONST u,T_SCALAR_VARIABLE_VIEW_CONST p, T_VECTOR_VARIABLE_VIEW f,T_SCALAR_VARIABLE_VIEW q) 
{
#ifdef LOG_DETAILED_PERFORMANCE_NE
    LOG::SCOPE scope("NONLINEAR_ELASTICITY::Add_Force_First_Order_Elasticity_Specialized");
#endif
#ifdef USE_SPECIALIZED_KERNELS

    typedef T (*u_compact_type)[3][27];
    typedef T (*p_compact_type)[8];
    u_compact_type u_or_f_compact=reinterpret_cast<u_compact_type>(specialized_data.u_compact_array.Get_Array_Pointer());
    p_compact_type p_or_q_compact=reinterpret_cast<p_compact_type>(specialized_data.p_compact_array.Get_Array_Pointer());

    CompactData_Specialized(u,p);

    Add_Force_First_Order_With_Specialized_Kernels<T,enable_constraints,enable_muscles>(
        u_or_f_compact,p_or_q_compact,
        cell_block_base_offsets.m,
        cell_block_partition_offsets.m,
        specialized_data);

    UnCompactData_Specialized(f,q);
#endif
}
//#####################################################################
// Function Update_Position_Based_State
//#####################################################################
template<class T,int d>  template<bool enable_constraints,bool enable_muscles> void NONLINEAR_ELASTICITY<T,d>::
Add_Force_Differential_Elasticity_Specialized(T_VECTOR_VARIABLE_VIEW_CONST du,T_SCALAR_VARIABLE_VIEW_CONST dp,T_VECTOR_VARIABLE_VIEW df,T_SCALAR_VARIABLE_VIEW dq) const
{
#ifdef LOG_DETAILED_PERFORMANCE_NE
    //LOG::SCOPE scope("NONLINEAR_ELASTICITY::Add_Force_Differential_Elasticity_Specialized");
#endif
#ifdef USE_SPECIALIZED_KERNELS

    typedef T (*u_compact_type)[3][27];
    typedef T (*p_compact_type)[8];

    u_compact_type du_or_df_compact=reinterpret_cast<u_compact_type>(specialized_data.u_compact_array.Get_Array_Pointer());
    p_compact_type dp_or_dq_compact=reinterpret_cast<p_compact_type>(specialized_data.p_compact_array.Get_Array_Pointer());

    CompactData_Specialized(du,dp);

    Add_Force_Differential_With_Specialized_Kernels<T,enable_constraints,enable_muscles>(
        du_or_df_compact,dp_or_dq_compact,
        cell_block_base_offsets.m,
        cell_block_partition_offsets.m,
        specialized_data );

    UnCompactData_Specialized(df,dq);
#endif

}
//#####################################################################
// Function Update_Position_Based_State
//#####################################################################
template<class T,int d>  template<bool enable_constraints,bool enable_muscles> void NONLINEAR_ELASTICITY<T,d>::
Add_Force_Differential_Elasticity_Specialized(T_VECTOR_VARIABLE_VIEW_CONST du,T_SCALAR_VARIABLE_VIEW_CONST dp,T_VECTOR_VARIABLE_VIEW df,T_SCALAR_VARIABLE_VIEW dq, const int subdomain) const
{
#ifdef LOG_DETAILED_PERFORMANCE_NE
    //LOG::SCOPE scope("NONLINEAR_ELASTICITY::Add_Force_Differential_Elasticity_Specialized");
#endif
#ifdef USE_SPECIALIZED_KERNELS

    typedef T (*u_compact_type)[3][27];
    typedef T (*p_compact_type)[8];

    u_compact_type du_or_df_compact=reinterpret_cast<u_compact_type>(specialized_data.u_compact_array.Get_Array_Pointer());
    p_compact_type dp_or_dq_compact=reinterpret_cast<p_compact_type>(specialized_data.p_compact_array.Get_Array_Pointer());

    CompactData_Specialized(du,dp);

    Add_Force_Differential_With_Specialized_Kernels<T,enable_constraints,enable_muscles>(
        du_or_df_compact,dp_or_dq_compact,
        cell_block_base_offsets.m,
        cell_block_partition_offsets.m,
        specialized_data );

    UnCompactData_Specialized(df,dq);
#endif

}
//#####################################################################
// Function Initialize_Blocks
//#####################################################################
template<class T,int d> void NONLINEAR_ELASTICITY<T,d>::
Initialize_Blocks_Specialized(const int number_of_partitions)
{
#ifdef USE_SPECIALIZED_KERNELS

//    static const int d=3;
    LOG::SCOPE scope("Initializing blocks");

    constant_partitions = number_of_partitions;

    T_FLAG block_is_active(padded_cell_domain);
    T_FLAG nodeblock_is_active(padded_node_domain);









    int number_of_active_cells=0;
    int number_of_active_blocks=0;
    int number_of_inactive_cells=0;
    int number_of_inactive_blocks=0;

    // Find out how many active cells/blocks we have

    for(RANGE_ITERATOR<d,2> block_iterator(unpadded_cell_domain);block_iterator.Valid();block_iterator.Next()){
        const T_INDEX& base_index=block_iterator.Index();
        for(RANGE_ITERATOR<d> cell_iterator(RANGE<T_INDEX>(base_index,base_index+1));cell_iterator.Valid();cell_iterator.Next()){
            const T_INDEX& cell_index=cell_iterator.Index();
            if(padded_cell_domain.Lazy_Outside(cell_index)) continue;
            if(cell_type(cell_index)==INTERIOR_CELL_TYPE || cell_type(cell_index)==BOUNDARY_CELL_TYPE){
                block_is_active(base_index)=true;
                number_of_active_cells++;
            }
            else
                number_of_inactive_cells++;
        }
        if(block_is_active(base_index)) number_of_active_blocks++; else number_of_inactive_blocks++;
    }

    LOG::cout<<"Total active cells = "<<number_of_active_cells<<std::endl;
    LOG::cout<<"  Active cell ratio: " << T(number_of_active_cells) / ( number_of_active_cells + number_of_inactive_cells) << std::endl;
    LOG::cout<<"Total active blocks = "<<number_of_active_blocks<<std::endl;
    LOG::cout<<"  Active block ratio: " << T(number_of_active_blocks) / ( number_of_active_blocks + number_of_inactive_blocks) << std::endl;



    LOG::cout << std::endl << std::endl;


    int number_of_active_nodes=0;
    int number_of_active_nodeblocks=0;
    int number_of_inactive_nodes=0;
    int number_of_inactive_nodeblocks=0;
    // Find out how many active nodes/blocks we have

    for(RANGE_ITERATOR<d,3> nodeblock_iterator(unpadded_node_domain);
        nodeblock_iterator.Valid();
        nodeblock_iterator.Next()){
        const T_INDEX& base_index=nodeblock_iterator.Index();
        for(RANGE_ITERATOR<d> node_iterator(RANGE<T_INDEX>(base_index,base_index+2));node_iterator.Valid();node_iterator.Next()){
            const T_INDEX& node_index=node_iterator.Index();
            if(padded_node_domain.Lazy_Outside(node_index)) continue;
            if(node_is_active(node_index)){
                nodeblock_is_active(base_index)=true;
                number_of_active_nodes++;
            }
            else
                number_of_inactive_nodes++;
        }
        if(nodeblock_is_active(base_index)) number_of_active_nodeblocks++; else number_of_inactive_nodeblocks++;
    }

    LOG::cout<<"Total active nodes = "<<number_of_active_nodes<<std::endl;
    LOG::cout<<"  Active node ratio: " << T(number_of_active_nodes) / ( number_of_active_nodes + number_of_inactive_nodes) << std::endl;
    LOG::cout<<"Total active nodeblocks = "<<number_of_active_nodeblocks<<std::endl;
    LOG::cout<<"  Active nodeblock ratio: " << T(number_of_active_nodeblocks) / ( number_of_active_nodeblocks + number_of_inactive_nodeblocks) << std::endl;


    LOG::cout << std::endl << std::endl;
















    // Generate partitions

    ARRAY<int,T_INDEX> block_partition_id(padded_cell_domain);
    HASHTABLE<T_INDEX> current_partition, next_partition;
    RANGE_ITERATOR<d,2> block_iterator(unpadded_cell_domain);
    int blocks_assigned=0;

    for(int partition=1;partition<=number_of_partitions;partition++){

        LOG::cout<<"Partition #"<<partition<<" :"<<std::endl;
        cell_block_partition_offsets.Append(blocks_assigned);

        // Compute optimal final block count

        int optimal_block_begin=(partition-1)*(number_of_active_blocks/number_of_partitions)+min(partition-1,number_of_active_blocks%number_of_partitions)+1;
        int optimal_block_end=partition*(number_of_active_blocks/number_of_partitions)+min(partition,number_of_active_blocks%number_of_partitions);

        LOG::cout<<"Ideal beginning block = "<<optimal_block_begin<<std::endl;
        LOG::cout<<"Actual beginning block = "<<blocks_assigned+1<<std::endl;
        LOG::cout<<"Ideal ending block = "<<optimal_block_end<<std::endl;

        // Process all the mandatory assignments

        HASHTABLE<T_INDEX>::Exchange_Hashtables(current_partition,next_partition);
        next_partition.Remove_All();

        for(HASHTABLE_ITERATOR<T_INDEX> hashtable_iterator(current_partition);hashtable_iterator.Valid();hashtable_iterator.Next()){
            const T_INDEX& index=hashtable_iterator.Key();
            if(!block_partition_id(index)){
                blocks_assigned++;
                block_partition_id(index)=partition;
                for(RANGE_ITERATOR<d,2> neighbor_iterator(RANGE<T_INDEX>(index-2,index+2));neighbor_iterator.Valid();neighbor_iterator.Next()){
                    const T_INDEX& neighbor_index=neighbor_iterator.Index();
                    if(unpadded_cell_domain.Lazy_Inside(neighbor_index) && block_is_active(neighbor_index) && !block_partition_id(neighbor_index))
                        next_partition.Set(neighbor_index);}}

        }

        LOG::cout<<"Ending block after mandatory insertions = "<<blocks_assigned<<std::endl;

        // Insert enough blocks from lexicographical order to meet quota

        while(blocks_assigned<optimal_block_end){

            PHYSBAM_ASSERT(block_iterator.Valid());
            const T_INDEX index=block_iterator.Index();
            block_iterator.Next();

            if(block_is_active(index) && !block_partition_id(index)){
                blocks_assigned++;
                block_partition_id(index)=partition;
                for(RANGE_ITERATOR<d,2> neighbor_iterator(RANGE<T_INDEX>(index-2,index+2));neighbor_iterator.Valid();neighbor_iterator.Next()){
                    const T_INDEX& neighbor_index=neighbor_iterator.Index();
                    if(unpadded_cell_domain.Lazy_Inside(neighbor_index) && block_is_active(neighbor_index) && !block_partition_id(neighbor_index))
                        next_partition.Set(neighbor_index);}}
           

        }

        LOG::cout<<"Ending block after mandatory and lexicographical insertions = "<<blocks_assigned<<std::endl;
    }







    LOG::cout << std::endl << std::endl;






    // Generate Node partitions

    ARRAY<int,T_INDEX> nodeblock_partition_id(padded_node_domain);
    HASHTABLE<T_INDEX> current_nodepartition, next_nodepartition;
    RANGE_ITERATOR<d,3> nodeblock_iterator(unpadded_node_domain);
    int nodeblocks_assigned=0;

    for(int partition=1;partition<=number_of_partitions;partition++){

        LOG::cout<<"Partition #"<<partition<<" :"<<std::endl;
        node_block_partition_offsets.Append(nodeblocks_assigned);

        // Compute optimal final nodeblock count

        int optimal_nodeblock_begin=(partition-1)*(number_of_active_nodeblocks/number_of_partitions)+min(partition-1,number_of_active_nodeblocks%number_of_partitions)+1;
        int optimal_nodeblock_end=partition*(number_of_active_nodeblocks/number_of_partitions)+min(partition,number_of_active_nodeblocks%number_of_partitions);

        LOG::cout<<"Ideal beginning nodeblock = "<<optimal_nodeblock_begin<<std::endl;
        LOG::cout<<"Actual beginning nodeblock = "<<nodeblocks_assigned+1<<std::endl;
        LOG::cout<<"Ideal ending nodeblock = "<<optimal_nodeblock_end<<std::endl;

        // Process all the mandatory assignments

        HASHTABLE<T_INDEX>::Exchange_Hashtables(current_nodepartition,next_nodepartition);
        next_nodepartition.Remove_All();

        for(HASHTABLE_ITERATOR<T_INDEX> hashtable_iterator(current_nodepartition);hashtable_iterator.Valid();hashtable_iterator.Next()){
            const T_INDEX& index=hashtable_iterator.Key();
            if(!nodeblock_partition_id(index)){
                nodeblocks_assigned++;
                nodeblock_partition_id(index)=partition;
                for(RANGE_ITERATOR<d,2> neighbor_iterator(RANGE<T_INDEX>(index-2,index+2));neighbor_iterator.Valid();neighbor_iterator.Next()){
                    const T_INDEX& neighbor_index=neighbor_iterator.Index();
                    if(unpadded_node_domain.Lazy_Inside(neighbor_index) && nodeblock_is_active(neighbor_index) && !nodeblock_partition_id(neighbor_index))
                        next_nodepartition.Set(neighbor_index);}}

        }

        LOG::cout<<"Ending nodeblock after mandatory insertions = "<<nodeblocks_assigned<<std::endl;

        // Insert enough nodeblocks from lexicographical order to meet quota

        while(nodeblocks_assigned<optimal_nodeblock_end){

            PHYSBAM_ASSERT(nodeblock_iterator.Valid());
            const T_INDEX index=nodeblock_iterator.Index();
            nodeblock_iterator.Next();

            if(nodeblock_is_active(index) && !nodeblock_partition_id(index)){
                nodeblocks_assigned++;
                nodeblock_partition_id(index)=partition;
                for(RANGE_ITERATOR<d,2> neighbor_iterator(RANGE<T_INDEX>(index-2,index+2));neighbor_iterator.Valid();neighbor_iterator.Next()){
                    const T_INDEX& neighbor_index=neighbor_iterator.Index();
                    if(unpadded_node_domain.Lazy_Inside(neighbor_index) && nodeblock_is_active(neighbor_index) && !nodeblock_partition_id(neighbor_index))
                        next_nodepartition.Set(neighbor_index);}}
           

        }

        LOG::cout<<"Ending nodeblock after mandatory and lexicographical insertions = "<<nodeblocks_assigned<<std::endl;
    }







    LOG::cout << std::endl << std::endl;















    LOG::cout<<"Total assigned blocks = "<<blocks_assigned<<std::endl;

    int missing_blocks = 0;
    for(;block_iterator.Valid();block_iterator.Next()){
        //   PHYSBAM_ASSERT(!block_is_active(block_iterator.Index()));
        if( block_is_active(block_iterator.Index()) )
            missing_blocks++;
    }

    LOG::cout << "Missing blocks: " << missing_blocks << std::endl;

    PHYSBAM_ASSERT(blocks_assigned==number_of_active_blocks);









    LOG::cout<<"Total assigned nodeblocks = "<<nodeblocks_assigned<<std::endl;

    int missing_nodeblocks = 0;
    for(;block_iterator.Valid();block_iterator.Next()){
        //   PHYSBAM_ASSERT(!block_is_active(block_iterator.Index()));
        if( block_is_active(block_iterator.Index()) )
            missing_nodeblocks++;
    }

    LOG::cout << "Missing nodeblocks: " << missing_nodeblocks << std::endl;

    PHYSBAM_ASSERT(nodeblocks_assigned==number_of_active_nodeblocks);











    // Order blocks within each partition

    ARRAY<int> blocks_in_partition(number_of_partitions);
    cell_block_base_indices.Resize(number_of_active_blocks);
    cell_block_base_offsets.Resize(number_of_active_blocks);

    VECTOR<int,3> grid_dims = padded_node_domain.max_corner - padded_node_domain.min_corner +1;
    int X_Stride = grid_dims(2)*grid_dims(3);
    int Y_Stride =              grid_dims(3);


    for(RANGE_ITERATOR<d,2> iterator(unpadded_cell_domain);iterator.Valid();iterator.Next()){
        const T_INDEX& index=iterator.Index();
        if(block_is_active(index)){
            const int partition=block_partition_id(index);
            cell_block_base_indices(++blocks_in_partition(partition)+cell_block_partition_offsets(partition))=index;
            cell_block_base_offsets(blocks_in_partition(partition)+cell_block_partition_offsets(partition))=index(1)*X_Stride + index(2)*Y_Stride + index(3);

      
        }}






    // Order nodeblocks within each partition

    ARRAY<int> nodeblocks_in_partition(number_of_partitions);
    node_block_base_indices.Resize(number_of_active_nodeblocks);
    node_block_base_offsets.Resize(number_of_active_nodeblocks);

    VECTOR<int,3> nodegrid_dims = padded_node_domain.max_corner - padded_node_domain.min_corner +1;
    int node_X_Stride = nodegrid_dims(2)*nodegrid_dims(3);
    int node_Y_Stride =              nodegrid_dims(3);


    for(RANGE_ITERATOR<d,3> iterator(unpadded_node_domain);iterator.Valid();iterator.Next()){
        const T_INDEX& index=iterator.Index();
        if(nodeblock_is_active(index)){
            const int partition=nodeblock_partition_id(index);
            node_block_base_indices(++nodeblocks_in_partition(partition)+node_block_partition_offsets(partition))=index;
            node_block_base_offsets(nodeblocks_in_partition(partition)+node_block_partition_offsets(partition))=index(1)*node_X_Stride + index(2)*node_Y_Stride + index(3);

      
        }}









    
    // Sanity check

    for(int partition=1;partition<=number_of_partitions;partition++){
        int block_begin=cell_block_partition_offsets(partition)+1;
        int block_end=(partition<number_of_partitions)?cell_block_partition_offsets(partition+1):number_of_active_blocks;
        for(int block=block_begin;block<=block_end;block++)
            PHYSBAM_ASSERT(block_partition_id(cell_block_base_indices(block))==partition);}


    for(int partition=1;partition<=number_of_partitions;partition++){
        int block_begin=node_block_partition_offsets(partition)+1;
        int block_end=(partition<number_of_partitions)?node_block_partition_offsets(partition+1):number_of_active_nodeblocks;
        for(int block=block_begin;block<=block_end;block++)
            PHYSBAM_ASSERT(nodeblock_partition_id(node_block_base_indices(block))==partition);}


    // Allocate flattened arrays
    ARRAY<BLOCKED_TYPE<T,8> > mu_flat;
    ARRAY<BLOCKED_TYPE<T,8> > mu_stab_flat;
    ARRAY<BLOCKED_TYPE<T,8> > kappa_flat;
    ARRAY<BLOCKED_TYPE<T,8> > alpha_flat;
    ARRAY<BLOCKED_TYPE<T,8> > alpha_sqr_over_kappa_flat;

    mu_flat.Resize(number_of_active_blocks);        
    mu_stab_flat.Resize(number_of_active_blocks);        
    kappa_flat.Resize(number_of_active_blocks);  
    alpha_sqr_over_kappa_flat.Resize(number_of_active_blocks);   
    alpha_flat.Resize(number_of_active_blocks);     

    // Initialize spatially varying material parameters

    for(int block=1;block<=cell_block_base_indices.m;block++){
        const T_INDEX& base_index=cell_block_base_indices(block);
        int cell=0;
        for(RANGE_ITERATOR<d> cell_iterator(RANGE<T_INDEX>(base_index,base_index+1));cell_iterator.Valid();cell_iterator.Next(),cell++){
            const T_INDEX& cell_index=cell_iterator.Index();
            if(cell_type(cell_index)==INTERIOR_CELL_TYPE || cell_type(cell_index)==BOUNDARY_CELL_TYPE){
                mu_flat(block).Set(constant_mu,cell);
                mu_stab_flat(block).Set(constant_mu_stab,cell);
                kappa_flat(block).Set(constant_kappa,cell);
                alpha_flat(block).Set(constant_alpha,cell);
                alpha_sqr_over_kappa_flat(block).Set((T)constant_alpha * constant_alpha / constant_kappa,cell);
            }
            else{
                mu_flat(block).Set(T(),cell);
                mu_stab_flat(block).Set(T(),cell);
                kappa_flat(block).Set(T(),cell);
                alpha_flat(block).Set(T(),cell);
                alpha_sqr_over_kappa_flat(block).Set(T(),cell);
            }}}

    // Initialize global material parameters


    specialized_data.InitializeMaterialParameters(mu_flat, mu_stab_flat,kappa_flat,alpha_flat,
                                                  alpha_sqr_over_kappa_flat,cutoff_value,
                                                  (T)h, (T)1.f/h, pow<d>(h));
    
#endif

}
//#####################################################################
template class NONLINEAR_ELASTICITY<float,3>;
template void NONLINEAR_ELASTICITY<float,3>::Update_Position_Based_State_Specialized<true,true>(T_VECTOR_VARIABLE_VIEW_CONST u,T_SCALAR_VARIABLE_VIEW_CONST p,T_VECTOR_VARIABLE_VIEW diag);
template void NONLINEAR_ELASTICITY<float,3>::Update_Position_Based_State_Specialized<true,false>(T_VECTOR_VARIABLE_VIEW_CONST u,T_SCALAR_VARIABLE_VIEW_CONST p,T_VECTOR_VARIABLE_VIEW diag);
template void NONLINEAR_ELASTICITY<float,3>::Update_Position_Based_State_Specialized<false,true>(T_VECTOR_VARIABLE_VIEW_CONST u,T_SCALAR_VARIABLE_VIEW_CONST p,T_VECTOR_VARIABLE_VIEW diag);
template void NONLINEAR_ELASTICITY<float,3>::Update_Position_Based_State_Specialized<false,false>(T_VECTOR_VARIABLE_VIEW_CONST u,T_SCALAR_VARIABLE_VIEW_CONST p,T_VECTOR_VARIABLE_VIEW diag);

template void NONLINEAR_ELASTICITY<float,3>::Add_Force_First_Order_Elasticity_Specialized<true,true>(T_VECTOR_VARIABLE_VIEW_CONST u,T_SCALAR_VARIABLE_VIEW_CONST p, T_VECTOR_VARIABLE_VIEW f,T_SCALAR_VARIABLE_VIEW q) ;
template void NONLINEAR_ELASTICITY<float,3>::Add_Force_First_Order_Elasticity_Specialized<true,false>(T_VECTOR_VARIABLE_VIEW_CONST u,T_SCALAR_VARIABLE_VIEW_CONST p, T_VECTOR_VARIABLE_VIEW f,T_SCALAR_VARIABLE_VIEW q) ;
template void NONLINEAR_ELASTICITY<float,3>::Add_Force_First_Order_Elasticity_Specialized<false,true>(T_VECTOR_VARIABLE_VIEW_CONST u,T_SCALAR_VARIABLE_VIEW_CONST p, T_VECTOR_VARIABLE_VIEW f,T_SCALAR_VARIABLE_VIEW q) ;
template void NONLINEAR_ELASTICITY<float,3>::Add_Force_First_Order_Elasticity_Specialized<false,false>(T_VECTOR_VARIABLE_VIEW_CONST u,T_SCALAR_VARIABLE_VIEW_CONST p, T_VECTOR_VARIABLE_VIEW f,T_SCALAR_VARIABLE_VIEW q) ;

template void NONLINEAR_ELASTICITY<float,3>::Add_Force_Differential_Elasticity_Specialized<true,true>(T_VECTOR_VARIABLE_VIEW_CONST du,T_SCALAR_VARIABLE_VIEW_CONST dp,T_VECTOR_VARIABLE_VIEW df,T_SCALAR_VARIABLE_VIEW dq) const;
template void NONLINEAR_ELASTICITY<float,3>::Add_Force_Differential_Elasticity_Specialized<true,false>(T_VECTOR_VARIABLE_VIEW_CONST du,T_SCALAR_VARIABLE_VIEW_CONST dp,T_VECTOR_VARIABLE_VIEW df,T_SCALAR_VARIABLE_VIEW dq) const;
template void NONLINEAR_ELASTICITY<float,3>::Add_Force_Differential_Elasticity_Specialized<false,true>(T_VECTOR_VARIABLE_VIEW_CONST du,T_SCALAR_VARIABLE_VIEW_CONST dp,T_VECTOR_VARIABLE_VIEW df,T_SCALAR_VARIABLE_VIEW dq) const;
template void NONLINEAR_ELASTICITY<float,3>::Add_Force_Differential_Elasticity_Specialized<false,false>(T_VECTOR_VARIABLE_VIEW_CONST du,T_SCALAR_VARIABLE_VIEW_CONST dp,T_VECTOR_VARIABLE_VIEW df,T_SCALAR_VARIABLE_VIEW dq) const;

template void NONLINEAR_ELASTICITY<float,3>::Add_Force_Differential_Elasticity_Specialized<true,true>(T_VECTOR_VARIABLE_VIEW_CONST du,T_SCALAR_VARIABLE_VIEW_CONST dp,T_VECTOR_VARIABLE_VIEW df,T_SCALAR_VARIABLE_VIEW dq, const int subdomain) const;
template void NONLINEAR_ELASTICITY<float,3>::Add_Force_Differential_Elasticity_Specialized<true,false>(T_VECTOR_VARIABLE_VIEW_CONST du,T_SCALAR_VARIABLE_VIEW_CONST dp,T_VECTOR_VARIABLE_VIEW df,T_SCALAR_VARIABLE_VIEW dq, const int subdomain) const;
template void NONLINEAR_ELASTICITY<float,3>::Add_Force_Differential_Elasticity_Specialized<false,true>(T_VECTOR_VARIABLE_VIEW_CONST du,T_SCALAR_VARIABLE_VIEW_CONST dp,T_VECTOR_VARIABLE_VIEW df,T_SCALAR_VARIABLE_VIEW dq, const int subdomain) const;
template void NONLINEAR_ELASTICITY<float,3>::Add_Force_Differential_Elasticity_Specialized<false,false>(T_VECTOR_VARIABLE_VIEW_CONST du,T_SCALAR_VARIABLE_VIEW_CONST dp,T_VECTOR_VARIABLE_VIEW df,T_SCALAR_VARIABLE_VIEW dq, const int subdomain) const;

#ifndef USE_SPECIALIZED_KERNELS
template class NONLINEAR_ELASTICITY<double,3>;
template void NONLINEAR_ELASTICITY<double,3>::Update_Position_Based_State_Specialized<true,true>(T_VECTOR_VARIABLE_VIEW_CONST u,T_SCALAR_VARIABLE_VIEW_CONST p,T_VECTOR_VARIABLE_VIEW diag);
template void NONLINEAR_ELASTICITY<double,3>::Update_Position_Based_State_Specialized<true,false>(T_VECTOR_VARIABLE_VIEW_CONST u,T_SCALAR_VARIABLE_VIEW_CONST p,T_VECTOR_VARIABLE_VIEW diag);
template void NONLINEAR_ELASTICITY<double,3>::Update_Position_Based_State_Specialized<false,true>(T_VECTOR_VARIABLE_VIEW_CONST u,T_SCALAR_VARIABLE_VIEW_CONST p,T_VECTOR_VARIABLE_VIEW diag);
template void NONLINEAR_ELASTICITY<double,3>::Update_Position_Based_State_Specialized<false,false>(T_VECTOR_VARIABLE_VIEW_CONST u,T_SCALAR_VARIABLE_VIEW_CONST p,T_VECTOR_VARIABLE_VIEW diag);

template void NONLINEAR_ELASTICITY<double,3>::Add_Force_First_Order_Elasticity_Specialized<true,true>(T_VECTOR_VARIABLE_VIEW_CONST u,T_SCALAR_VARIABLE_VIEW_CONST p, T_VECTOR_VARIABLE_VIEW f,T_SCALAR_VARIABLE_VIEW q) ;
template void NONLINEAR_ELASTICITY<double,3>::Add_Force_First_Order_Elasticity_Specialized<true,false>(T_VECTOR_VARIABLE_VIEW_CONST u,T_SCALAR_VARIABLE_VIEW_CONST p, T_VECTOR_VARIABLE_VIEW f,T_SCALAR_VARIABLE_VIEW q) ;
template void NONLINEAR_ELASTICITY<double,3>::Add_Force_First_Order_Elasticity_Specialized<false,true>(T_VECTOR_VARIABLE_VIEW_CONST u,T_SCALAR_VARIABLE_VIEW_CONST p, T_VECTOR_VARIABLE_VIEW f,T_SCALAR_VARIABLE_VIEW q) ;
template void NONLINEAR_ELASTICITY<double,3>::Add_Force_First_Order_Elasticity_Specialized<false,false>(T_VECTOR_VARIABLE_VIEW_CONST u,T_SCALAR_VARIABLE_VIEW_CONST p, T_VECTOR_VARIABLE_VIEW f,T_SCALAR_VARIABLE_VIEW q) ;

template void NONLINEAR_ELASTICITY<double,3>::Add_Force_Differential_Elasticity_Specialized<true,true>(T_VECTOR_VARIABLE_VIEW_CONST du,T_SCALAR_VARIABLE_VIEW_CONST dp,T_VECTOR_VARIABLE_VIEW df,T_SCALAR_VARIABLE_VIEW dq) const;
template void NONLINEAR_ELASTICITY<double,3>::Add_Force_Differential_Elasticity_Specialized<true,false>(T_VECTOR_VARIABLE_VIEW_CONST du,T_SCALAR_VARIABLE_VIEW_CONST dp,T_VECTOR_VARIABLE_VIEW df,T_SCALAR_VARIABLE_VIEW dq) const;
template void NONLINEAR_ELASTICITY<double,3>::Add_Force_Differential_Elasticity_Specialized<false,true>(T_VECTOR_VARIABLE_VIEW_CONST du,T_SCALAR_VARIABLE_VIEW_CONST dp,T_VECTOR_VARIABLE_VIEW df,T_SCALAR_VARIABLE_VIEW dq) const;
template void NONLINEAR_ELASTICITY<double,3>::Add_Force_Differential_Elasticity_Specialized<false,false>(T_VECTOR_VARIABLE_VIEW_CONST du,T_SCALAR_VARIABLE_VIEW_CONST dp,T_VECTOR_VARIABLE_VIEW df,T_SCALAR_VARIABLE_VIEW dq) const;

template void NONLINEAR_ELASTICITY<double,3>::Add_Force_Differential_Elasticity_Specialized<true,true>(T_VECTOR_VARIABLE_VIEW_CONST du,T_SCALAR_VARIABLE_VIEW_CONST dp,T_VECTOR_VARIABLE_VIEW df,T_SCALAR_VARIABLE_VIEW dq, const int subdomain) const;
template void NONLINEAR_ELASTICITY<double,3>::Add_Force_Differential_Elasticity_Specialized<true,false>(T_VECTOR_VARIABLE_VIEW_CONST du,T_SCALAR_VARIABLE_VIEW_CONST dp,T_VECTOR_VARIABLE_VIEW df,T_SCALAR_VARIABLE_VIEW dq, const int subdomain) const;
template void NONLINEAR_ELASTICITY<double,3>::Add_Force_Differential_Elasticity_Specialized<false,true>(T_VECTOR_VARIABLE_VIEW_CONST du,T_SCALAR_VARIABLE_VIEW_CONST dp,T_VECTOR_VARIABLE_VIEW df,T_SCALAR_VARIABLE_VIEW dq, const int subdomain) const;
template void NONLINEAR_ELASTICITY<double,3>::Add_Force_Differential_Elasticity_Specialized<false,false>(T_VECTOR_VARIABLE_VIEW_CONST du,T_SCALAR_VARIABLE_VIEW_CONST dp,T_VECTOR_VARIABLE_VIEW df,T_SCALAR_VARIABLE_VIEW dq, const int subdomain) const;
#endif

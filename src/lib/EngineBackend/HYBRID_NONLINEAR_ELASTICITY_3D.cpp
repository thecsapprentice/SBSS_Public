#include "HYBRID_NONLINEAR_ELASTICITY.h"
#include <Common/RANGE_ITERATOR.h>
#include "SPECIALIZED_KERNELS.h"
#include <PhysBAM_Tools/Data_Structures/HASHTABLE.h>
#include <PhysBAM_Tools/Data_Structures/HASHTABLE_ITERATOR.h>
#include <PhysBAM_Tools/Data_Structures/UNION_FIND.h>
#include <PhysBAM_Tools/Data_Structures/OPERATION_HASH.h>
#include <PhysBAM_Tools/Data_Structures/TRIPLE.h>
#include <PhysBAM_Tools/Arrays_Computations/SORT.h>


#include <Thread_Queueing/PTHREAD_QUEUE.h>
using namespace PhysBAM;
extern PTHREAD_QUEUE* pthread_queue;


#ifdef ENABLE_MIC
#define ENABLE_MIC_INSTRUCTION_SET
#endif
#include <SIMD_Optimized_Kernels/Kernels/Hybrid_Grid_Compact/Hybrid_Grid_Compact.h>
#include <SIMD_Optimized_Kernels/Kernels/Hybrid_Grid_UnCompact/Hybrid_Grid_UnCompact.h>
#ifdef ENABLE_MIC_INSTRUCTION_SET
#undef ENABLE_MIC_INSTRUCTION_SET
#endif



namespace BLOCK_HELPERS{
    const int COMPACT_CELL_MAP[] = { 0, 1, 3, 4, 9, 10, 12, 13 };
};

#ifdef LOG_DETAILED_PERFORMANCE
#define LOG_DETAILED_PERFORMANCE_NE
#endif

template<class T_NONLINEAR_ELASTICITY>
HYBRID_NONLINEAR_ELASTICITY_STATE<T_NONLINEAR_ELASTICITY>::HYBRID_NONLINEAR_ELASTICITY_STATE() :
    x_mesh(T_SCALAR_VARIABLE_MESH_VIEW(0,NULL),T_SCALAR_VARIABLE_MESH_VIEW(0,NULL),T_SCALAR_VARIABLE_MESH_VIEW(0,NULL)), p_mesh(0,NULL)
{}

template<class T_NONLINEAR_ELASTICITY>
typename HYBRID_NONLINEAR_ELASTICITY<T_NONLINEAR_ELASTICITY>::T_VECTOR_VARIABLE_MESH_VIEW_CONST 
HYBRID_NONLINEAR_ELASTICITY<T_NONLINEAR_ELASTICITY>::
View_Convert(const VECTOR<T_SCALAR_VARIABLE_MESH_VIEW,d>& in)
{
    return VECTOR<T_SCALAR_VARIABLE_MESH_VIEW_CONST,d>(in.x,in.y,in.z);
};

//#####################################################################
// Function Initialize_Blocks_Elasticity
//#####################################################################

void Block_Merge( VECTOR<int,8>& Merge_Tox, VECTOR<int,27>& Merge_Toy, const VECTOR<int,8>& Merge_Fromx, const VECTOR<int,27>& Merge_Fromy )
{
    for( int i = 1; i <=8; i++){
        PHYSBAM_ASSERT( Merge_Tox(i) == Merge_Fromx(i) || Merge_Tox(i) == -1 || Merge_Fromx(i) == -1);
        if(Merge_Fromx(i) != -1) Merge_Tox(i) = Merge_Fromx(i);
    }

    for( int i = 1; i <=27; i++){
        PHYSBAM_ASSERT( Merge_Toy(i) == Merge_Fromy(i) || Merge_Toy(i) == -1 || Merge_Fromy(i) == -1);
        if(Merge_Fromy(i) != -1) Merge_Toy(i) = Merge_Fromy(i);
    }
}

bool Block_Connects( const VECTOR<int,8>& Ax, const VECTOR<int,27>& Ay, const VECTOR<int,8>& Bx, const VECTOR<int,27>& By )
{
    for( int i = 1; i <=8; i++)
        if( Ax(i) != Bx(i) && Ax(i) != -1 && Bx(i) != -1 )
            return false;
    for( int i = 1; i <=27; i++)
        if( Ay(i) != By(i) && Ay(i) != -1 && By(i) != -1 )
            return false;
    return true;
}


void Collapse_Blocks(HASHTABLE<int,ARRAY<int> >& hash, HASHTABLE<int,ARRAY<int> >& clusters)
{
    UNION_FIND<int> union_find;
    clusters.Remove_All();
    union_find.Initialize(hash.Size());
    OPERATION_HASH<int> operation_hash(hash.Size());
    int initial_elements = hash.Size();

    //LOG::cout << "Pre-Iteration.\n" << std::endl;
    //hash.Print_Table(LOG::cout);

    int passes = 0;
    while(1){
        //LOG::cout << "\n---------------------------------------\nIteration " << ++passes << std::endl; 
        VECTOR<int,2> best_pair;
        int best_mismatches=hash.Size();
    
        for(HASHTABLE_ITERATOR<int,ARRAY<int> > it1(hash);it1.Valid();it1.Next())
            for(HASHTABLE_ITERATOR<int,ARRAY<int> > it2(hash);it2.Valid();it2.Next()){
                const int& key1=it1.Key(),key2=it2.Key();
                if(key1<key2){
                    const ARRAY<int> &array1=it1.Data(),&array2=it2.Data();
                    operation_hash.Next_Operation();
                    for(int m=1;m<=array1.m;m++)
                        operation_hash.Mark(array1(m));
                    int matches=0;
                    int is_valid=0;
                    for(int m=1;m<=array2.m;m++)
                        if(operation_hash.Is_Marked_Current(array2(m))){
                            matches++;
                            if(array2(m)==key1 || array2(m)==key2)
                                is_valid++;
                        }
                    PHYSBAM_ASSERT( is_valid <=2 && is_valid >=0 );
                    int mismatches=array1.m+array2.m-2*matches;
                    if(mismatches<best_mismatches && matches>0 && is_valid==2){
                        best_mismatches=mismatches;
                        best_pair=VECTOR<int,2>(key1,key2);
                    }
                }
            }
        //LOG::cout << "Best pair in this iteration: " << best_pair << std::endl;
        //LOG::cout << "Best mismatches in this iteration: " << best_mismatches  << std::endl;
        if(best_pair==VECTOR<int,2>()) break; // no off-diagonals

        // Collapse
        int i,j;best_pair.Get(i,j);
        union_find.Union(i,j);

        for(HASHTABLE_ITERATOR<int,ARRAY<int> > it(hash);it.Valid();it.Next()){
            ARRAY<int>& array=it.Data();
            int ij_found=0;
            for(int k=array.m;k>=1;k--)
                if(array(k)==i || array(k)==j){
                    array.Remove_Index_Lazy(k);
                    ij_found++;
                }
            if(ij_found==2)
                array.Append(i);
        }
        hash.Delete(j);
        //Clean up step to preserve symmetry
        ARRAY<int> array = hash.Get(i);
        ARRAY<int> clean_array;
        for(int m=1;m <= array.m; m++){
            const ARRAY<int>& array2 = hash.Get(array(m));
            int i_found=0;
            for(int k=1;k<=array2.m;k++)
                if(array2(k) == i ) i_found++;
            if(i_found)
                clean_array.Append(array(m));
        }
        hash.Set(i, clean_array);

        //LOG::cout << std::endl;
        //hash.Print_Table(LOG::cout);
    }

    // LOG::cout << "\n\nAfter Last Iteration" << std::endl;
    // hash.Print_Table(LOG::cout);

    HASHTABLE<int,int> cluster_roots;
    int cluster_count=0;
    for(int m = 1; m <= initial_elements; m++){
        int cluster_root = union_find.Find(m);
        if( !cluster_roots.Contains(cluster_root) ){
            cluster_roots.Insert(cluster_root,++cluster_count);
            clusters.Insert(cluster_count,ARRAY<int>());
        }
        clusters.Get(cluster_roots.Get(cluster_root)).Append(m);
    }      

    //LOG::cout << "\n\nFinal Clusters" << std::endl;
    //clusters.Print_Table(LOG::cout);

}

struct HYBRID_BLOCK_SORTER{
    typedef VECTOR<int,3> T_INDEX;
    typedef T_INDEX HYBRID_BLOCK_BASE;
    typedef ARRAY<VECTOR<int,8> > HYBRID_BLOCK_CELLS;
    typedef ARRAY<VECTOR<int,27> > HYBRID_BLOCK_NODES;
    typedef TRIPLE<HYBRID_BLOCK_BASE, HYBRID_BLOCK_CELLS, HYBRID_BLOCK_NODES > HYBRID_BLOCK_STACK; 

    const ARRAY<int,T_INDEX> * _block_partition_id;
    HYBRID_BLOCK_SORTER(const ARRAY<int,T_INDEX>* bpi) : _block_partition_id( bpi ) {};

    bool operator() (const HYBRID_BLOCK_STACK& A, const HYBRID_BLOCK_STACK& B){
        if( (*_block_partition_id)(A.x) < (*_block_partition_id)(B.x) )
            return true;
        else if( (*_block_partition_id)(A.x) == (*_block_partition_id)(B.x) && LEXICOGRAPHIC_COMPARE()(A.x,B.x) )
            return true;
        
        return false;
    }
};

template<class T_NONLINEAR_ELASTICITY> void HYBRID_NONLINEAR_ELASTICITY<T_NONLINEAR_ELASTICITY>::
Initialize_Blocks_Elasticity(int number_of_partitions)
{

#ifdef USE_SPECIALIZED_KERNELS

    LOG::SCOPE scope("HYBRID_NONLINEAR_ELASTICITY::Initialize_Blocks_Elasticity()");
    
    constant_partitions = number_of_partitions;

    ARRAY< ARRAY<int>, T_INDEX> cell_bins;
    HASHTABLE< T_INDEX, int > block_index;
    ARRAY< HYBRID_BLOCK_STACK > collapsed_blocks;
    number_of_blocks = 0;
    cell_bins.Resize( padded_cell_domain );
    
    for(RANGE_ITERATOR<d> cell_iterator(padded_cell_domain); cell_iterator.Valid(); cell_iterator.Next()){
        const T_INDEX& cell_index = cell_iterator.Index();
        if(cell_type(cell_index) == INTERIOR_CELL_TYPE || cell_type(cell_index) == BOUNDARY_CELL_TYPE || cell_type(cell_index) == DIRICHLET_CELL_TYPE)
            cell_bins(cell_index).Append( 0 );
    }
    
    for( int m = 1; m <= Number_Of_Mesh_Cells(); m++){
        const T_INDEX& cell_index = cell_indices_mesh(m);
        if(cell_type_mesh(m) == INTERIOR_CELL_TYPE || cell_type_mesh(m) == BOUNDARY_CELL_TYPE || cell_type_mesh(m) == DIRICHLET_CELL_TYPE)
            cell_bins(cell_index).Append( m );
    }
    
    for(RANGE_ITERATOR<d,2> intrablock_iterator(unpadded_cell_domain); intrablock_iterator.Valid(); intrablock_iterator.Next()){
        const T_INDEX& block_base = intrablock_iterator.Index();
        //LOG::SCOPE scope("Examining block base");
        //LOG::cout << block_base << std::endl;
        
        collapsed_blocks.Append(  HYBRID_BLOCK_STACK() );
        block_index.Insert( block_base, collapsed_blocks.m );
        
        HYBRID_BLOCK_CELLS hbc;
        HYBRID_BLOCK_NODES hbn;
          
        int cell_in_block = 1;
        for(RANGE_ITERATOR<d> interblock_iterator(RANGE<T_INDEX>(block_base,block_base+1)); interblock_iterator.Valid(); interblock_iterator.Next(), cell_in_block++){
            const T_INDEX& cell_index = interblock_iterator.Index();
            //LOG::SCOPE scope("Examining cell");
            //LOG::cout << cell_index << std::endl;

            //LOG::cout << "This cell has " << cell_bins(cell_index).m << " duplicates." << std::endl;
            for( int c = 1; c <= cell_bins(cell_index).m; c++){
                hbc.Append(VECTOR<int,8>::All_Ones_Vector()*-1);
                hbn.Append(VECTOR<int,27>::All_Ones_Vector()*-1);              
                
                hbc.Last()(cell_in_block) = cell_bins(cell_index)(c);
                int node_in_cell = 1;
                for(RANGE_ITERATOR<d> node_iterator(RANGE<T_INDEX>(T_INDEX(),T_INDEX()+1)); node_iterator.Valid(); node_iterator.Next(), node_in_cell++){
                    if( cell_bins(cell_index)(c) == 0 ){ // A Grid Cell
                        hbn.Last()( BLOCK_HELPERS::COMPACT_CELL_MAP[cell_in_block-1] + BLOCK_HELPERS::COMPACT_CELL_MAP[node_in_cell-1] + 1 ) = 0;
                    }
                    else{ // A Mesh Cell
                        hbn.Last()( BLOCK_HELPERS::COMPACT_CELL_MAP[cell_in_block-1] + BLOCK_HELPERS::COMPACT_CELL_MAP[node_in_cell-1] + 1 ) = cells_mesh(cell_bins(cell_index)(c))(node_in_cell);
                    }
                }
            }
        }
        
        HASHTABLE<int,ARRAY<int> > connectivity;
        HASHTABLE<int,ARRAY<int> > clusters;
        
        PHYSBAM_ASSERT( hbc.m == hbn.m );
        
        //LOG::cout << "\n-------------------- DETERMINING CONNECTIVITY ---------------------------\n" << std::endl;
        for( int m = 1; m <= hbc.m; m++)
            for( int n = 1; n <= hbc.m; n++ ){
                //LOG::cout << "Determining Connectivity for pair " << m << " and " << n;
                if( Block_Connects( hbc(m), hbn(m), hbc(n), hbn(n))){
                    if( !connectivity.Contains(m) )
                        connectivity.Insert(m,ARRAY<int>());
                    connectivity.Get(m).Append(n);
                    //LOG::cout << "   Connected!";
                }
                //LOG::cout << std::endl;
            }
        
        //LOG::cout << "\n-------------------- INITIAL CONNECTIVITY ---------------------------\n" << std::endl;
        //connectivity.Print_Table(LOG::cout);

        //LOG::cout << "\n-------------------- COLLAPSING BLOCKS ---------------------------\n" << std::endl;
        Collapse_Blocks(connectivity, clusters);
        
        HYBRID_BLOCK_BASE&  c_hbb = collapsed_blocks(block_index.Get(block_base)).x;
        HYBRID_BLOCK_CELLS& c_hbc = collapsed_blocks(block_index.Get(block_base)).y;
        HYBRID_BLOCK_NODES& c_hbn = collapsed_blocks(block_index.Get(block_base)).z;

        c_hbb = block_base;
        for(HASHTABLE_ITERATOR<int,ARRAY<int> > it(clusters);it.Valid();it.Next()){
            c_hbc.Append(VECTOR<int,8>::All_Ones_Vector()*-1);
            c_hbn.Append(VECTOR<int,27>::All_Ones_Vector()*-1);  
            for( int m = 1; m <= it.Data().m; m++)
                Block_Merge( c_hbc.Last(), c_hbn.Last(), hbc(it.Data()(m)), hbn(it.Data()(m)));
        }
        number_of_blocks+= c_hbc.m;
    }

    // Generate Partitions

    ARRAY<int,T_INDEX> block_partition_id(padded_cell_domain);
    HASHTABLE<T_INDEX> current_partition, next_partition;
    RANGE_ITERATOR<d,2> block_iterator(unpadded_cell_domain);
    int blocks_assigned=0;

    for(int partition=1;partition<=number_of_partitions;partition++){

        LOG::cout<<"\n\n\nPartition #"<<partition<<" :"<<std::endl;
        cell_block_partition_offsets.Append(blocks_assigned);

        //PHYSBAM_ASSERT( blocks_assigned <= number_of_blocks );

        int optimal_block_begin=(partition-1)*(number_of_blocks/number_of_partitions)+min(partition-1,number_of_blocks%number_of_partitions)+1;
        int optimal_block_end=partition*(number_of_blocks/number_of_partitions)+min(partition,number_of_blocks%number_of_partitions);

        LOG::cout<<"Ideal beginning block = "<<optimal_block_begin<<std::endl;
        LOG::cout<<"Actual beginning block = "<<blocks_assigned+1<<std::endl;
        LOG::cout<<"Ideal ending block = "<<optimal_block_end<<std::endl;

        // Process all the mandatory assignments

        HASHTABLE<T_INDEX>::Exchange_Hashtables(current_partition,next_partition);
        next_partition.Remove_All();

        for(HASHTABLE_ITERATOR<T_INDEX> hashtable_iterator(current_partition);hashtable_iterator.Valid();hashtable_iterator.Next()){
            const T_INDEX& index=hashtable_iterator.Key();
            const int bx = block_index.Get_Or_Insert(index,0);
            if(!block_partition_id(index)){
                blocks_assigned+=collapsed_blocks(bx).y.m; // Increase block count by number of blocks that share this base index               
                block_partition_id(index)=partition;
                for(RANGE_ITERATOR<d,2> neighbor_iterator(RANGE<T_INDEX>(index-2,index+2));neighbor_iterator.Valid();neighbor_iterator.Next()){
                    const T_INDEX& neighbor_index=neighbor_iterator.Index();
                    const int nx = block_index.Get_Or_Insert(neighbor_index,0);
                    if(unpadded_cell_domain.Lazy_Inside(neighbor_index) && nx && !block_partition_id(neighbor_index))
                        next_partition.Set(neighbor_index);                    
                }}            
        }

        LOG::cout<<"Ending block after mandatory insertions = "<<blocks_assigned<<std::endl;

        // Insert enough blocks from lexicographical order to meet quota

        while(blocks_assigned<optimal_block_end){
            //PHYSBAM_ASSERT(block_iterator.Valid());
            if(!block_iterator.Valid())
                break;
            const T_INDEX index=block_iterator.Index();
            block_iterator.Next();
            const int bx = block_index.Get_Or_Insert(index,0);
            if(  bx != 0 && !block_partition_id(index)){                    
                blocks_assigned+=collapsed_blocks(bx).y.m; // Increase block count by number of blocks that share this base index
                block_partition_id(index)=partition;
                for(RANGE_ITERATOR<d,2> neighbor_iterator(RANGE<T_INDEX>(index-2,index+2));neighbor_iterator.Valid();neighbor_iterator.Next()){
                    const T_INDEX& neighbor_index=neighbor_iterator.Index();
                    const int nx = block_index.Get_Or_Insert(neighbor_index,0);
                    if(unpadded_cell_domain.Lazy_Inside(neighbor_index) && nx && !block_partition_id(neighbor_index))
                        next_partition.Set(neighbor_index);}}
           

        }

        LOG::cout<<"Ending block after mandatory and lexicographical insertions = "<<blocks_assigned<<std::endl;
    }


    // Sort collapsed_blocks to make them line up with the ordering in partitions.

    Sort( collapsed_blocks, HYBRID_BLOCK_SORTER(&block_partition_id) );

#if 0
    LOG::cout << "Partition spread: " << std::endl;
    for( int p = 1; p < number_of_partitions; p++)
        LOG::cout << "Partition " << p << "(" << cell_block_partition_offsets(p) << "): " << cell_block_partition_offsets(p+1) - cell_block_partition_offsets(p) << std::endl;
    LOG::cout << "Partition " << number_of_partitions << "(" << cell_block_partition_offsets(number_of_partitions) << "): " << number_of_blocks - cell_block_partition_offsets(number_of_partitions) << std::endl;

    LOG::cout << "Partition correctness " << std::endl;
    for(RANGE_ITERATOR<d,2> intrablock_iterator(unpadded_cell_domain); intrablock_iterator.Valid(); intrablock_iterator.Next()){
        const T_INDEX& index = intrablock_iterator.Index();
        for(RANGE_ITERATOR<d,2> neighbor_iterator(RANGE<T_INDEX>(index-2,index+2));neighbor_iterator.Valid();neighbor_iterator.Next()){
            const T_INDEX& neighbor_index = neighbor_iterator.Index();
            const int bx = block_index.Get_Or_Insert(index,0);
            const int nx = block_index.Get_Or_Insert(neighbor_index,0);
            if(unpadded_cell_domain.Lazy_Inside(neighbor_index))
                try{
                    PHYSBAM_ASSERT( block_partition_id(index) == block_partition_id(neighbor_index) ||
                                    abs(block_partition_id(index) - block_partition_id(neighbor_index))==1 ||
                                    block_partition_id(neighbor_index) == 0);
                }
                catch(...){
                    LOG::cout << index << std::endl;
                    LOG::cout << neighbor_index << std::endl;
                    
                    LOG::cout << block_partition_id(index) << std::endl;
                    LOG::cout << block_partition_id(neighbor_index) << std::endl;
                    exit(1);
                }
        }
    }
#endif

    // Assign compacted blocks to master arrays

    for( int cb = 1; cb <= collapsed_blocks.m; cb++ ){
        for( int m = 1; m <= collapsed_blocks(cb).y.m; m++ ){
            block_bases.Append( collapsed_blocks(cb).x );
            hybrid_block_cells.Append( collapsed_blocks(cb).y(m) );
            hybrid_block_nodes.Append( collapsed_blocks(cb).z(m) );
        }
    }

#if 1
    // Build Block Reverse Maps
    hybrid_block_grid_reverse_map.Resize(padded_cell_domain);
    hybrid_block_mesh_reverse_map.Resize(number_of_mesh_cells);
    hybrid_block_grid_reverse_map.Fill( VECTOR<int,2>(-1,-1) );
    hybrid_block_mesh_reverse_map.Fill( VECTOR<int,2>(-1,-1) );
    for(int block=1;block<=number_of_blocks;block++){
        const T_INDEX& base_index=block_bases(block);
        const VECTOR<int,8>& block_cell=hybrid_block_cells(block); 
        int cell=1;
        for(RANGE_ITERATOR<d> cell_iterator(RANGE<T_INDEX>(base_index,base_index+1));cell_iterator.Valid();cell_iterator.Next(),cell++){
            if( block_cell(cell) == 0 )
                hybrid_block_grid_reverse_map( cell_iterator.Index() ) = VECTOR<int,2>( block, cell );
            if( block_cell(cell) > 0 )
                hybrid_block_mesh_reverse_map( block_cell(cell) ) = VECTOR<int,2>( block, cell );           
        }
    }
#endif


#ifdef ENABLE_MIC
    // Fill out the special mic offset/mask arrays for fast access
    {        
        VECTOR<int,3> cell_grid_dims = padded_cell_domain.max_corner - padded_cell_domain.min_corner +1;
        VECTOR<int,3> node_grid_dims = padded_node_domain.max_corner - padded_node_domain.min_corner +1 ;
        
        int X_StrideCell = cell_grid_dims(2)*cell_grid_dims(3);
        int Y_StrideCell =                   cell_grid_dims(3);
        
        int X_StrideNode = node_grid_dims(2)*node_grid_dims(3);
        int Y_StrideNode =                   node_grid_dims(3);

        mic_block_offsets_low = (int*)(_mm_malloc (number_of_blocks*16*sizeof(int), 64));
        mic_block_offsets_high = (int*)(_mm_malloc (number_of_blocks*16*sizeof(int), 64));
        mic_block_offsets_grid_mask = (int*)(_mm_malloc (number_of_blocks*sizeof(int), 64)); 
        mic_block_offsets_mesh_mask = (int*)(_mm_malloc (number_of_blocks*sizeof(int), 64));       
        for( int block = 0 ; block < number_of_blocks; block++){
            const T_INDEX& block_base = block_bases(block+1);
            const VECTOR<int,27> block_nodes = hybrid_block_nodes(block+1);
            int base_index = block_base.x * X_StrideNode + 
                             block_base.y * Y_StrideNode + 
                             block_base.z; 
            int linear_index=0;
            mic_block_offsets_grid_mask[block] = 0;
            mic_block_offsets_mesh_mask[block] = 0;
            for(int i=0;i<3;i++)
                for(int j=0;j<3;j++)
                    for(int k=0;k<3;k++){
                        int node_index = base_index + (i)*X_StrideNode+(j)*Y_StrideNode+(k)*1;
                        int mesh_index = block_nodes(linear_index+1)-1; 
                        int toggle_bit;
                        switch( block_nodes(linear_index+1) ){
                        case -1:
                            if( linear_index < 16 )
                                mic_block_offsets_low[ (block*16)+linear_index ] = 0;
                            else
                                mic_block_offsets_high[ (block*16)+(linear_index-16) ] = 0;
                            break;
                        case 0:
                            if( linear_index < 16 )
                                mic_block_offsets_low[ (block*16)+linear_index ] = node_index;
                            else
                                mic_block_offsets_high[ (block*16)+(linear_index-16) ] = node_index;
                            toggle_bit = 1 << linear_index;
                            mic_block_offsets_grid_mask[block] |= toggle_bit;
                            break;
                        default:                            
                            if( linear_index < 16 )
                                mic_block_offsets_low[ (block*16)+linear_index ] = mesh_index;
                            else
                                mic_block_offsets_high[ (block*16)+(linear_index-16) ] = mesh_index;
                            toggle_bit = 1 << linear_index;
                            mic_block_offsets_mesh_mask[block] |= toggle_bit;
                            break;
                        }

                        linear_index++;
                    }            
        }
    }
#endif

    PHYSBAM_ASSERT( block_bases.m == hybrid_block_cells.m );
    PHYSBAM_ASSERT( block_bases.m == hybrid_block_nodes.m );
    PHYSBAM_ASSERT( number_of_blocks == block_bases.m );
    LOG::cout << "Generated " << number_of_blocks << " total active blocks." << std::endl;


    // Allocate flattened arrays
    ARRAY<BLOCKED_TYPE<T,8> > mu_flat;
    ARRAY<BLOCKED_TYPE<T,8> > mu_stab_flat;
    ARRAY<BLOCKED_TYPE<T,8> > kappa_flat;
    ARRAY<BLOCKED_TYPE<T,8> > alpha_flat;
    ARRAY<BLOCKED_TYPE<T,8> > alpha_sqr_over_kappa_flat;

    mu_flat.Resize(number_of_blocks);        
    mu_stab_flat.Resize(number_of_blocks);        
    kappa_flat.Resize(number_of_blocks);  
    alpha_sqr_over_kappa_flat.Resize(number_of_blocks);   
    alpha_flat.Resize(number_of_blocks);     

    // Initialize spatially varying material parameters

    for(int block=1;block<=number_of_blocks;block++){
        const T_INDEX& base_index=block_bases(block);
        const VECTOR<int,8>& block_cell=hybrid_block_cells(block);

        int cell=0;
        for(RANGE_ITERATOR<d> cell_iterator(RANGE<T_INDEX>(base_index,base_index+1));cell_iterator.Valid();cell_iterator.Next(),cell++){
            const T_INDEX& cell_index=cell_iterator.Index();
            switch(block_cell(cell+1)){
            case -1:
                mu_flat(block).Set(T(),cell);
                mu_stab_flat(block).Set(T(),cell);
                kappa_flat(block).Set(T(),cell);
                alpha_flat(block).Set(T(),cell);
                alpha_sqr_over_kappa_flat(block).Set(T(),cell);
                break;
            case 0:
                if(cell_type(cell_index)==INTERIOR_CELL_TYPE || cell_type(cell_index)==BOUNDARY_CELL_TYPE){
                    mu_flat(block).Set(BASE::constant_mu,cell);
                    mu_stab_flat(block).Set(BASE::constant_mu_stab,cell);
                    kappa_flat(block).Set(BASE::constant_kappa,cell);
                    alpha_flat(block).Set(BASE::constant_alpha,cell);
                    alpha_sqr_over_kappa_flat(block).Set((T)BASE::constant_alpha * BASE::constant_alpha / BASE::constant_kappa,cell);
                }
                else{
                    mu_flat(block).Set(T(),cell);
                    mu_stab_flat(block).Set(T(),cell);
                    kappa_flat(block).Set(T(),cell);
                    alpha_flat(block).Set(T(),cell);
                    alpha_sqr_over_kappa_flat(block).Set(T(),cell);
                }
                break;
            default:
                int mesh_index=block_cell(cell+1);
                if(cell_type_mesh(mesh_index)==INTERIOR_CELL_TYPE || cell_type_mesh(mesh_index)==BOUNDARY_CELL_TYPE){
                    mu_flat(block).Set(BASE::constant_mu,cell);
                    mu_stab_flat(block).Set(BASE::constant_mu_stab,cell);
                    kappa_flat(block).Set(BASE::constant_kappa,cell);
                    alpha_flat(block).Set(BASE::constant_alpha,cell);
                    alpha_sqr_over_kappa_flat(block).Set((T)BASE::constant_alpha * BASE::constant_alpha / BASE::constant_kappa,cell);
                }
                else{
                    mu_flat(block).Set(T(),cell);
                    mu_stab_flat(block).Set(T(),cell);
                    kappa_flat(block).Set(T(),cell);
                    alpha_flat(block).Set(T(),cell);
                    alpha_sqr_over_kappa_flat(block).Set(T(),cell);
                }
                break;
            }
        }
    }

    BASE::specialized_data.InitializeMaterialParameters(mu_flat, mu_stab_flat,kappa_flat,alpha_flat,
                                                        alpha_sqr_over_kappa_flat,BASE::cutoff_value,
                                                        (T)h,(T)1.f/h, pow<d>(h));


#endif

    // Generate CG blocks (Only grid. We pass mesh stuff through as partitioned serial lists.)

    VECTOR<int,3> nodegrid_dims = padded_node_domain.max_corner - padded_node_domain.min_corner +1;
    int node_X_Stride = nodegrid_dims(2)*nodegrid_dims(3);
    int node_Y_Stride =              nodegrid_dims(3);
    VECTOR<int,3> grid_dims = padded_node_domain.max_corner - padded_node_domain.min_corner +1;
    int X_Stride = grid_dims(2)*grid_dims(3);
    int Y_Stride =              grid_dims(3);

    // This needs to be true. If not, then padding is wrong and needs to be fixed.
    PHYSBAM_ASSERT( node_X_Stride == X_Stride );
    PHYSBAM_ASSERT( node_Y_Stride == Y_Stride );


 // This was originally written with an internal u variable. But we can fake it just fine
#if 0
    VECTOR<T*,d+1> cg_base;
    VECTOR<int,d+1> cg_max_index;
    for(int w=1;w<=d;w++){
        cg_base(w) = &BASE::internal_state.x(w)(0,0,0);
        PHYSBAM_ASSERT( cg_base(w) == BASE::internal_state.x(w).array.Get_Array_Pointer());
        cg_max_index(w) = nodegrid_dims(1) * nodegrid_dims(2) * nodegrid_dims(3);
    }
    cg_base(d+1) = &BASE::internal_state.p(0,0,0);
    PHYSBAM_ASSERT( cg_base(d+1) == BASE::internal_state.p.array.Get_Array_Pointer());
    cg_max_index(d+1) = grid_dims(1) * grid_dims(2) * grid_dims(3);
    PHYSBAM_ASSERT(cg_max_index(d+1) == cg_max_index(d));
#else
    VECTOR<T*,d+1> cg_base;
    VECTOR<int,d+1> cg_max_index;
    for(int w=1;w<=d;w++){
        cg_base(w) = 0x0;
        cg_max_index(w) = nodegrid_dims(1) * nodegrid_dims(2) * nodegrid_dims(3);
    }
    cg_base(d+1) = 0x0;
    cg_max_index(d+1) = grid_dims(1) * grid_dims(2) * grid_dims(3);
#endif


    const int STRIDE = 16;
    {
        LOG::SCOPE scope( "Generating Cache Line Offsets - Grid" );

        LOG::cout << "There are " << cg_max_index(1)/STRIDE << " cache lines to examine." << std::endl;
        
        if(cg_max_index(1) % STRIDE != 0){
            LOG::cout << "   There are uneven cache lines! This could be a problem." << std::endl;
            LOG::cout << "   " << cg_max_index(1) % STRIDE << " extra floats..." << std::endl;
        }
        
        int i;
        int max_index = cg_max_index(1) % STRIDE ? cg_max_index(1)-(cg_max_index(1)%STRIDE) : cg_max_index(1);
        for(i=0; i<max_index; i+=STRIDE){
            bool line_active=false;
            for( int k=i; k<i+STRIDE; k++){
                T_INDEX ND_index;
                PHYSBAM_ASSERT( (k) < cg_max_index(1) );
                ND_index(1) = k / node_X_Stride;
                ND_index(2) = (k % node_X_Stride) / node_Y_Stride;
                ND_index(3) = ((k % node_X_Stride) % node_Y_Stride);
                if(node_is_active(ND_index) || node_is_dirichlet(ND_index) )
                    line_active = true;                
            }
            if(line_active)
                cg_offsets.Append(i);
        }
        LOG::cout << "There are " << cg_offsets.m << " active cache lines." << std::endl;
    }
    {
        LOG::SCOPE scope( "Generating Cache Line Offsets - Mesh:Cell" );
        int i;
        int max_index = Number_Of_Mesh_Cells() % STRIDE ? Number_Of_Mesh_Cells()+(STRIDE - (Number_Of_Mesh_Cells()%STRIDE)) : Number_Of_Mesh_Cells(); // We preadjust these elsewhere...

        for(i=0; i<max_index; i+=STRIDE){
            bool line_active=false;
            for( int k=i; k<i+STRIDE; k++){
                if( (k) < Number_Of_Mesh_Cells() )
                    if(/*cell_type_mesh(k+1) != EXTERIOR_CELL_TYPE*/true)
                        line_active = true;                
            }
            if(line_active)
                cg_offsets_cell_mesh.Append(i);
        }
        LOG::cout << "There are " << cg_offsets_cell_mesh.m << " active cache lines." << std::endl;
    }
    {
        LOG::SCOPE scope( "Generating Cache Line Offsets - Mesh:Node" );
        int i;
        int max_index = Number_Of_Mesh_Nodes() % STRIDE ? Number_Of_Mesh_Nodes()+(STRIDE - (Number_Of_Mesh_Nodes()%STRIDE)) : Number_Of_Mesh_Nodes(); // We preadjust these elsewhere...

        for(i=0; i<max_index; i+=STRIDE){
            bool line_active=false;
            for( int k=i; k<i+STRIDE; k++){
                if( (k) < Number_Of_Mesh_Nodes())
                    if(/*node_is_active_mesh(k+1)*/true)
                        line_active = true;                
            }
            if(line_active)
                cg_offsets_node_mesh.Append(i);
        }
        LOG::cout << "There are " << cg_offsets_node_mesh.m << " active cache lines." << std::endl;
    }

    ARRAY<bool,T_INDEX> cgblock_is_active;
    cgblock_base_offsets.Resize(0);
    cgblock_is_active.Resize(unpadded_node_domain);


    const int CG_BLOCK_SIZE = 3;  

    for(RANGE_ITERATOR<d,CG_BLOCK_SIZE> cgblock_iterator(unpadded_node_domain);
        cgblock_iterator.Valid();
        cgblock_iterator.Next()){
        const T_INDEX& base_index=cgblock_iterator.Index();
        for(RANGE_ITERATOR<d> cg_iterator(RANGE<T_INDEX>(base_index,base_index+CG_BLOCK_SIZE-1));cg_iterator.Valid();cg_iterator.Next()){
            const T_INDEX& cg_index=cg_iterator.Index();
            if(padded_node_domain.Lazy_Outside(cg_index)) continue;
            if(node_is_active(cg_index))
                cgblock_is_active(base_index) = true;               
        }
        
        if(cgblock_is_active(base_index))
            cgblock_base_offsets.Append( base_index.x * X_Stride + base_index.y * Y_Stride + base_index.z );
    }
    LOG::cout << "There are " << cgblock_is_active.Number_True() << " active blocks for CG operations." << std::endl;


}


//#####################################################################
// Function CompactData_Specialized
//#####################################################################
template<class T_NONLINEAR_ELASTICITY> void HYBRID_NONLINEAR_ELASTICITY<T_NONLINEAR_ELASTICITY>::
CompactData_Specialized(T_VECTOR_VARIABLE_VIEW_CONST x, T_SCALAR_VARIABLE_VIEW_CONST p, T_VECTOR_VARIABLE_MESH_VIEW_CONST x_mesh, T_SCALAR_VARIABLE_MESH_VIEW_CONST p_mesh) const 
{
#ifdef LOG_DETAILED_PERFORMANCE_NE
    LOG::SCOPE scope("HYBRID_NONLINEAR_ELASTICITY::CompactData_Specialized");
#endif
#ifdef USE_SPECIALIZED_KERNELS

    VECTOR<int,3> cell_grid_dims = padded_cell_domain.max_corner - padded_cell_domain.min_corner +1;
    VECTOR<int,3> node_grid_dims = padded_node_domain.max_corner - padded_node_domain.min_corner +1 ;

    int X_StrideCell = cell_grid_dims(2)*cell_grid_dims(3);
    int Y_StrideCell =                   cell_grid_dims(3);

    int X_StrideNode = node_grid_dims(2)*node_grid_dims(3);
    int Y_StrideNode =                   node_grid_dims(3);
    
    typedef T (*u_compact_type)[3][27];
    typedef T (*p_compact_type)[8];
    u_compact_type u_or_f_compact=reinterpret_cast<u_compact_type>(BASE::specialized_data.u_compact_array.Get_Array_Pointer());
    p_compact_type p_or_q_compact=reinterpret_cast<p_compact_type>(BASE::specialized_data.p_compact_array.Get_Array_Pointer());


    Hybrid_Grid_Compact(u_or_f_compact,
                        p_or_q_compact,
                        &x(1)(0,0,0),
                        &x(2)(0,0,0),
                        &x(3)(0,0,0),
                        &p(0,0,0),
                        x_mesh(1).Get_Array_Pointer(),
                        x_mesh(2).Get_Array_Pointer(),
                        x_mesh(3).Get_Array_Pointer(),
                        p_mesh.Get_Array_Pointer(),
                        number_of_blocks,
#ifdef ENABLE_MIC
                        mic_block_offsets_low,
                        mic_block_offsets_high,
                        mic_block_offsets_grid_mask,
                        mic_block_offsets_mesh_mask,
#endif
                        reinterpret_cast<const int*>(block_bases.Get_Array_Pointer()),
                        reinterpret_cast<const int*>(hybrid_block_cells.Get_Array_Pointer()),
                        reinterpret_cast<const int*>(hybrid_block_nodes.Get_Array_Pointer()),
                        constant_partitions,
                        X_StrideCell, Y_StrideCell,
                        X_StrideNode, Y_StrideNode,
                        false);

#endif
}
//#####################################################################
// Function UnCompactData_Specialized
//#####################################################################
template<class T_NONLINEAR_ELASTICITY> void HYBRID_NONLINEAR_ELASTICITY<T_NONLINEAR_ELASTICITY>::
UnCompactData_Specialized(T_VECTOR_VARIABLE_VIEW x,T_SCALAR_VARIABLE_VIEW p, T_VECTOR_VARIABLE_MESH_VIEW x_mesh,T_SCALAR_VARIABLE_MESH_VIEW p_mesh) const 
{
#ifdef LOG_DETAILED_PERFORMANCE_NE
    LOG::SCOPE scope("HYBRID_NONLINEAR_ELASTICITY::UnCompactData_Specialized");
#endif
#ifdef USE_SPECIALIZED_KERNELS

    VECTOR<int,3> cell_grid_dims = padded_cell_domain.max_corner - padded_cell_domain.min_corner +1;
    VECTOR<int,3> node_grid_dims = padded_node_domain.max_corner - padded_node_domain.min_corner +1 ;

    int X_StrideCell = cell_grid_dims(2)*cell_grid_dims(3);
    int Y_StrideCell =                   cell_grid_dims(3);

    int X_StrideNode = node_grid_dims(2)*node_grid_dims(3);
    int Y_StrideNode =                   node_grid_dims(3);
    
    typedef T (*u_compact_type)[3][27];
    typedef T (*p_compact_type)[8];
    u_compact_type u_or_f_compact=reinterpret_cast<u_compact_type>(BASE::specialized_data.u_compact_array.Get_Array_Pointer());
    p_compact_type p_or_q_compact=reinterpret_cast<p_compact_type>(BASE::specialized_data.p_compact_array.Get_Array_Pointer());


    Hybrid_Grid_UnCompact(u_or_f_compact,
                          p_or_q_compact,
                          &x(1)(0,0,0),
                          &x(2)(0,0,0),
                          &x(3)(0,0,0),
                          &p(0,0,0),
                          x_mesh(1).Get_Array_Pointer(),
                          x_mesh(2).Get_Array_Pointer(),
                          x_mesh(3).Get_Array_Pointer(),
                          p_mesh.Get_Array_Pointer(),
                          number_of_blocks,
#ifdef ENABLE_MIC
                          mic_block_offsets_low,
                          mic_block_offsets_high,
                          mic_block_offsets_grid_mask,
                          mic_block_offsets_mesh_mask,
#endif
                          reinterpret_cast<const int*>(block_bases.Get_Array_Pointer()),
                          reinterpret_cast<const int*>(hybrid_block_cells.Get_Array_Pointer()),
                          reinterpret_cast<const int*>(hybrid_block_nodes.Get_Array_Pointer()),
                          cell_block_partition_offsets.Get_Array_Pointer(),cell_block_partition_offsets.m,
                          X_StrideCell, Y_StrideCell,
                          X_StrideNode, Y_StrideNode,
                          false);

#endif
}
//#####################################################################
// Function CompactData_Specialized
//#####################################################################
template<class T_NONLINEAR_ELASTICITY> void HYBRID_NONLINEAR_ELASTICITY<T_NONLINEAR_ELASTICITY>::
CompactData_Specialized(T_VECTOR_VARIABLE_VIEW_CONST x, T_VECTOR_VARIABLE_MESH_VIEW_CONST x_mesh) const 
{
#ifdef LOG_DETAILED_PERFORMANCE_NE
    LOG::SCOPE scope("HYBRID_NONLINEAR_ELASTICITY::CompactData_Specialized");
#endif
#ifdef USE_SPECIALIZED_KERNELS

    VECTOR<int,3> cell_grid_dims = padded_cell_domain.max_corner - padded_cell_domain.min_corner +1;
    VECTOR<int,3> node_grid_dims = padded_node_domain.max_corner - padded_node_domain.min_corner +1 ;

    int X_StrideCell = cell_grid_dims(2)*cell_grid_dims(3);
    int Y_StrideCell =                   cell_grid_dims(3);

    int X_StrideNode = node_grid_dims(2)*node_grid_dims(3);
    int Y_StrideNode =                   node_grid_dims(3);
    
    typedef T (*u_compact_type)[3][27];
    typedef T (*p_compact_type)[8];
    u_compact_type d_compact=reinterpret_cast<u_compact_type>(BASE::specialized_data.d_compact_array.Get_Array_Pointer());
    p_compact_type null_compact=reinterpret_cast<p_compact_type>(NULL);

    Hybrid_Grid_Compact(d_compact,
                        null_compact,
                        &x(1)(0,0,0),
                        &x(2)(0,0,0),
                        &x(3)(0,0,0),
                        (T*)(NULL),
                        x_mesh(1).Get_Array_Pointer(),
                        x_mesh(2).Get_Array_Pointer(),
                        x_mesh(3).Get_Array_Pointer(),
                        (T*)(NULL),
                        number_of_blocks,
#ifdef ENABLE_MIC
                        mic_block_offsets_low,
                        mic_block_offsets_high,
                        mic_block_offsets_grid_mask,
                        mic_block_offsets_mesh_mask,
#endif
                        reinterpret_cast<const int*>(block_bases.Get_Array_Pointer()),
                        reinterpret_cast<const int*>(hybrid_block_cells.Get_Array_Pointer()),
                        reinterpret_cast<const int*>(hybrid_block_nodes.Get_Array_Pointer()),
                        cell_block_partition_offsets.m,
                        X_StrideCell, Y_StrideCell,
                        X_StrideNode, Y_StrideNode,
                        false);
#endif
}
//#####################################################################
// Function UnCompactData_Specialized
//#####################################################################
template<class T_NONLINEAR_ELASTICITY> void HYBRID_NONLINEAR_ELASTICITY<T_NONLINEAR_ELASTICITY>::
UnCompactData_Specialized(T_VECTOR_VARIABLE_VIEW x, T_VECTOR_VARIABLE_MESH_VIEW x_mesh) const 
{
#ifdef LOG_DETAILED_PERFORMANCE_NE
    LOG::SCOPE scope("HYBRID_NONLINEAR_ELASTICITY::UnCompactData_Specialized");
#endif
#ifdef USE_SPECIALIZED_KERNELS

    VECTOR<int,3> cell_grid_dims = padded_cell_domain.max_corner - padded_cell_domain.min_corner +1;
    VECTOR<int,3> node_grid_dims = padded_node_domain.max_corner - padded_node_domain.min_corner +1 ;

    int X_StrideCell = cell_grid_dims(2)*cell_grid_dims(3);
    int Y_StrideCell =                   cell_grid_dims(3);

    int X_StrideNode = node_grid_dims(2)*node_grid_dims(3);
    int Y_StrideNode =                   node_grid_dims(3);
    
    typedef T (*u_compact_type)[3][27];
    typedef T (*p_compact_type)[8];
    u_compact_type d_compact=reinterpret_cast<u_compact_type>(BASE::specialized_data.d_compact_array.Get_Array_Pointer());
    p_compact_type null_compact=reinterpret_cast<p_compact_type>(NULL);

    Hybrid_Grid_UnCompact(d_compact,
                          null_compact,
                          &x(1)(0,0,0),
                          &x(2)(0,0,0),
                          &x(3)(0,0,0),
                          (T*)(NULL),
                          x_mesh(1).Get_Array_Pointer(),
                          x_mesh(2).Get_Array_Pointer(),
                          x_mesh(3).Get_Array_Pointer(),
                          (T*)(NULL),
                          number_of_blocks,
#ifdef ENABLE_MIC
                          mic_block_offsets_low,
                          mic_block_offsets_high,
                          mic_block_offsets_grid_mask,
                          mic_block_offsets_mesh_mask,
#endif
                          reinterpret_cast<const int*>(block_bases.Get_Array_Pointer()),
                          reinterpret_cast<const int*>(hybrid_block_cells.Get_Array_Pointer()),
                          reinterpret_cast<const int*>(hybrid_block_nodes.Get_Array_Pointer()),
                          cell_block_partition_offsets.Get_Array_Pointer(),cell_block_partition_offsets.m,
                          X_StrideCell, Y_StrideCell,
                          X_StrideNode, Y_StrideNode,
                          false);
#endif
}
//#####################################################################
// Function Update_Position_Based_State
//#####################################################################
template<class T_NONLINEAR_ELASTICITY> void HYBRID_NONLINEAR_ELASTICITY<T_NONLINEAR_ELASTICITY>::
Update_Position_Based_State_Specialized(T_VECTOR_VARIABLE_VIEW_CONST u,T_SCALAR_VARIABLE_VIEW_CONST p,T_VECTOR_VARIABLE_VIEW diag,T_VECTOR_VARIABLE_MESH_VIEW_CONST u_mesh,T_SCALAR_VARIABLE_MESH_VIEW_CONST p_mesh,T_VECTOR_VARIABLE_MESH_VIEW diag_mesh)
{
#ifdef LOG_DETAILED_PERFORMANCE_NE
    LOG::SCOPE scope("HYBRID_NONLINEAR_ELASTICITY::Update_Position_Based_State_Specialized");
#endif
    
#ifdef USE_SPECIALIZED_KERNELS
    CompactData_Specialized(u,p,u_mesh,p_mesh);
    typedef T (*u_compact_type)[3][27];
    typedef T (*p_compact_type)[8];
    u_compact_type u_compact=reinterpret_cast<u_compact_type>(BASE::specialized_data.u_compact_array.Get_Array_Pointer());
    p_compact_type p_compact=reinterpret_cast<p_compact_type>(BASE::specialized_data.p_compact_array.Get_Array_Pointer());
    u_compact_type d_compact=reinterpret_cast<u_compact_type>(BASE::specialized_data.d_compact_array.Get_Array_Pointer());

#ifdef LOG_DETAILED_PERFORMANCE_NE
    {
        LOG::SCOPE scope("HYBRID_NONLINEAR_ELASTICITY::Update_Position_Based_State_Specialized_KERNEL");
#endif
    Update_Position_Based_State_With_Specialized_Kernels<T,supports_constraints,supports_muscles>(
        u_compact,p_compact,d_compact,
        number_of_blocks,
        constant_partitions,
        muscle_fiber_max_stresses,muscle_activations,
        BASE::specialized_data);
#ifdef LOG_DETAILED_PERFORMANCE_NE
    }
#endif

    //UnCompactData_Specialized(diag,diag_mesh);
#endif
}
//#####################################################################
// Function Update_Position_Based_State
//#####################################################################
template<class T_NONLINEAR_ELASTICITY> void HYBRID_NONLINEAR_ELASTICITY<T_NONLINEAR_ELASTICITY>::
Add_Force_First_Order_Elasticity_Specialized(T_VECTOR_VARIABLE_VIEW_CONST u,T_SCALAR_VARIABLE_VIEW_CONST p,T_VECTOR_VARIABLE_MESH_VIEW_CONST u_mesh,T_SCALAR_VARIABLE_MESH_VIEW_CONST p_mesh, T_VECTOR_VARIABLE_VIEW f,T_SCALAR_VARIABLE_VIEW q,T_VECTOR_VARIABLE_MESH_VIEW f_mesh,T_SCALAR_VARIABLE_MESH_VIEW q_mesh) 
{
#ifdef LOG_DETAILED_PERFORMANCE_NE
    LOG::SCOPE scope("HYBRID_NONLINEAR_ELASTICITY::Add_Force_First_Order_Elasticity_Specialized");
#endif
#ifdef USE_SPECIALIZED_KERNELS

    typedef T (*u_compact_type)[3][27];
    typedef T (*p_compact_type)[8];
    u_compact_type u_or_f_compact=reinterpret_cast<u_compact_type>(BASE::specialized_data.u_compact_array.Get_Array_Pointer());
    p_compact_type p_or_q_compact=reinterpret_cast<p_compact_type>(BASE::specialized_data.p_compact_array.Get_Array_Pointer());

    CompactData_Specialized(u,p,u_mesh,p_mesh);

#ifdef LOG_DETAILED_PERFORMANCE_NE
    {
        LOG::SCOPE scope("HYBRID_NONLINEAR_ELASTICITY::Add_Force_First_Order_Specialized_KERNEL");
#endif
    Add_Force_First_Order_With_Specialized_Kernels<T,supports_constraints,supports_muscles>(
        u_or_f_compact,p_or_q_compact,
        number_of_blocks,
        constant_partitions,
        BASE::specialized_data);
#ifdef LOG_DETAILED_PERFORMANCE_NE
    }
#endif

    UnCompactData_Specialized(f,q,f_mesh,q_mesh);
#endif
}
//#####################################################################
// Function Update_Position_Based_State
//#####################################################################
template<class T_NONLINEAR_ELASTICITY> void HYBRID_NONLINEAR_ELASTICITY<T_NONLINEAR_ELASTICITY>::
Add_Force_Differential_Elasticity_Specialized(T_VECTOR_VARIABLE_VIEW_CONST du,T_SCALAR_VARIABLE_VIEW_CONST dp,T_VECTOR_VARIABLE_MESH_VIEW_CONST du_mesh,T_SCALAR_VARIABLE_MESH_VIEW_CONST dp_mesh,
                                              T_VECTOR_VARIABLE_VIEW df,T_SCALAR_VARIABLE_VIEW dq, T_VECTOR_VARIABLE_MESH_VIEW df_mesh,T_SCALAR_VARIABLE_MESH_VIEW dq_mesh) const
{
#ifdef LOG_DETAILED_PERFORMANCE_NE
    LOG::SCOPE scope("HYBRID_NONLINEAR_ELASTICITY::Add_Force_Differential_Elasticity_Specialized");
#endif
#ifdef USE_SPECIALIZED_KERNELS

    typedef T (*u_compact_type)[3][27];
    typedef T (*p_compact_type)[8];

    u_compact_type du_or_df_compact=reinterpret_cast<u_compact_type>(BASE::specialized_data.u_compact_array.Get_Array_Pointer());
    p_compact_type dp_or_dq_compact=reinterpret_cast<p_compact_type>(BASE::specialized_data.p_compact_array.Get_Array_Pointer());

    CompactData_Specialized(du,dp,du_mesh,dp_mesh);


#ifdef LOG_DETAILED_PERFORMANCE_NE
    {
        LOG::SCOPE scope("HYBRID_NONLINEAR_ELASTICITY::Add_Force_Differential_Elasticity_Specialized_KERNEL");
#endif
        Add_Force_Differential_With_Specialized_Kernels<T,supports_constraints,supports_muscles>(
                      du_or_df_compact,dp_or_dq_compact,
                      number_of_blocks,
                      constant_partitions,
                      BASE::specialized_data );
#ifdef LOG_DETAILED_PERFORMANCE_NE
    }
#endif

    UnCompactData_Specialized(df,dq,df_mesh,dq_mesh);
#endif

}
//#####################################################################
// Function Update_Position_Based_State
//#####################################################################
template<class T_NONLINEAR_ELASTICITY> void HYBRID_NONLINEAR_ELASTICITY<T_NONLINEAR_ELASTICITY>::
Add_Force_Differential_Elasticity_Specialized(T_VECTOR_VARIABLE_VIEW_CONST du,T_SCALAR_VARIABLE_VIEW_CONST dp,T_VECTOR_VARIABLE_MESH_VIEW_CONST du_mesh,T_SCALAR_VARIABLE_MESH_VIEW_CONST dp_mesh,
                                              T_VECTOR_VARIABLE_VIEW df,T_SCALAR_VARIABLE_VIEW dq, T_VECTOR_VARIABLE_MESH_VIEW df_mesh,T_SCALAR_VARIABLE_MESH_VIEW dq_mesh, const int subdomain) const
{
#ifdef LOG_DETAILED_PERFORMANCE_NE
    LOG::SCOPE scope("HYBRID_NONLINEAR_ELASTICITY::Add_Force_Differential_Elasticity_Specialized");
#endif
	PHYSBAM_NOT_IMPLEMENTED();
#if 0
#ifdef USE_SPECIALIZED_KERNELS

    typedef T (*u_compact_type)[3][27];
    typedef T (*p_compact_type)[8];

    u_compact_type du_or_df_compact=reinterpret_cast<u_compact_type>(BASE::specialized_data.u_compact_array.Get_Array_Pointer());
    p_compact_type dp_or_dq_compact=reinterpret_cast<p_compact_type>(BASE::specialized_data.p_compact_array.Get_Array_Pointer());

    ARRAY<int> list_of_subdomain_blocks;
    int new_partition_count;

    for(int block=1; block<=number_of_blocks; block++){
        const T_INDEX& base_index = block_bases(block);
        int color = cell_color(base_index);
        if( color == subdomain )
            list_of_subdomain_blocks.Append( block );
    }
    // Change partition count such that the new partition sizes are equivilent to the old ones.
    new_partition_count = list_of_subdomain_blocks.m / (number_of_blocks / cell_block_partition_offsets.m);


    CompactData_Specialized(du,dp,du_mesh,dp_mesh);


#ifdef LOG_DETAILED_PERFORMANCE_NE
    {
        LOG::SCOPE scope("HYBRID_NONLINEAR_ELASTICITY::Add_Force_Differential_Elasticity_Specialized_KERNEL");
#endif
        Add_Force_Differential_With_Specialized_Kernels_Domain<T,supports_constraints,supports_muscles>(
                      du_or_df_compact,dp_or_dq_compact,
                      list_of_subdomain_blocks,
                      new_partition_count,
                      BASE::specialized_data );
#ifdef LOG_DETAILED_PERFORMANCE_NE
    }
#endif

    UnCompactData_Specialized(df,dq,df_mesh,dq_mesh);
#endif
#endif
}
//#####################################################################
// Function Stress_Mesh
//#####################################################################
template<class T_NONLINEAR_ELASTICITY>
HYBRID_NONLINEAR_ELASTICITY<T_NONLINEAR_ELASTICITY>::TV
HYBRID_NONLINEAR_ELASTICITY<T_NONLINEAR_ELASTICITY>::
Stress_Mesh( const int& mesh_element, const TV& multilinear_coordinates ) {
    int block = hybrid_block_mesh_reverse_map( mesh_element )(1);
    int cell =  hybrid_block_mesh_reverse_map( mesh_element )(2);
    int bundle = ((block-1) / BASE::specialized_data.VECTOR_MULT) + 1;
    int bundle_offset = (block-1) % BASE::specialized_data.VECTOR_MULT;
    int bundle_cell = (cell-1) + 8 * bundle_offset;

    DIAGONAL_MATRIX<T,d> result;
    BASE::specialized_data.P_hat_bundled(bundle).Get( result, bundle_cell );
    
    return TV( result.x11, result.x22, result.x33 );
}
//#####################################################################
// Function Stress_Grid
//#####################################################################
template<class T_NONLINEAR_ELASTICITY>
HYBRID_NONLINEAR_ELASTICITY<T_NONLINEAR_ELASTICITY>::TV
HYBRID_NONLINEAR_ELASTICITY<T_NONLINEAR_ELASTICITY>::
Stress_Grid( const T_INDEX& cell_index, const TV& multilinear_coordinates ) {
    int block = hybrid_block_grid_reverse_map( cell_index )(1);
    int cell =  hybrid_block_grid_reverse_map( cell_index )(2);
    int bundle = ((block-1) / BASE::specialized_data.VECTOR_MULT) + 1;
    int bundle_offset = (block-1) % BASE::specialized_data.VECTOR_MULT;
    int bundle_cell = (cell-1) + 8 * bundle_offset;

    DIAGONAL_MATRIX<T,d> result;
    BASE::specialized_data.P_hat_bundled(bundle).Get( result, bundle_cell );
    
    return TV( result.x11, result.x22, result.x33 );
}
//#####################################################################
// Function Strain_Mesh
//#####################################################################
template<class T_NONLINEAR_ELASTICITY>
HYBRID_NONLINEAR_ELASTICITY<T_NONLINEAR_ELASTICITY>::TV
HYBRID_NONLINEAR_ELASTICITY<T_NONLINEAR_ELASTICITY>::
Strain_Mesh( const int& mesh_element, const TV& multilinear_coordinates ) {

    int block = hybrid_block_mesh_reverse_map( mesh_element )(1);
    int cell =  hybrid_block_mesh_reverse_map( mesh_element )(2);
    int bundle = ((block-1) / BASE::specialized_data.VECTOR_MULT) + 1;
    int bundle_offset = (block-1) % BASE::specialized_data.VECTOR_MULT;
    int bundle_cell = (cell-1) + 8 * bundle_offset;

    DIAGONAL_MATRIX<T,d> result;
    BASE::specialized_data.Sigma_bundled(bundle).Get( result, bundle_cell );
    
    return TV( result.x11-1, result.x22-1, result.x33-1 );

}
//#####################################################################
// Function Strain_Grid
//#####################################################################
template<class T_NONLINEAR_ELASTICITY>
HYBRID_NONLINEAR_ELASTICITY<T_NONLINEAR_ELASTICITY>::TV
HYBRID_NONLINEAR_ELASTICITY<T_NONLINEAR_ELASTICITY>::
Strain_Grid( const T_INDEX& cell_index, const TV& multilinear_coordinates ) {

    int block = hybrid_block_grid_reverse_map( cell_index )(1);
    int cell =  hybrid_block_grid_reverse_map( cell_index )(2);
    int bundle = ((block-1) / BASE::specialized_data.VECTOR_MULT) + 1;
    int bundle_offset = (block-1) % BASE::specialized_data.VECTOR_MULT;
    int bundle_cell = (cell-1) + 8 * bundle_offset;

    DIAGONAL_MATRIX<T,d> result;
    BASE::specialized_data.Sigma_bundled(bundle).Get( result, bundle_cell );
    
    return TV( result.x11-1, result.x22-1, result.x33-1 );
}


template struct HYBRID_NONLINEAR_ELASTICITY_STATE<NONLINEAR_ELASTICITY<float,3> >;
template struct HYBRID_NONLINEAR_ELASTICITY_STATE<SKINNING_NONLINEAR_ELASTICITY<float,3,true,true> >;
template struct HYBRID_NONLINEAR_ELASTICITY_STATE<SKINNING_NONLINEAR_ELASTICITY<float,3,true,false> >;
template class HYBRID_NONLINEAR_ELASTICITY<NONLINEAR_ELASTICITY<float,3> >;
template class HYBRID_NONLINEAR_ELASTICITY<SKINNING_NONLINEAR_ELASTICITY<float,3,true,true> >;
template class HYBRID_NONLINEAR_ELASTICITY<SKINNING_NONLINEAR_ELASTICITY<float,3,true,false> >;
#ifndef USE_SPECIALIZED_KERNELS
template struct HYBRID_NONLINEAR_ELASTICITY_STATE<NONLINEAR_ELASTICITY<double,3> >;
template struct HYBRID_NONLINEAR_ELASTICITY_STATE<SKINNING_NONLINEAR_ELASTICITY<double,3,true,true> >;
template struct HYBRID_NONLINEAR_ELASTICITY_STATE<SKINNING_NONLINEAR_ELASTICITY<double,3,true,false> >;
template class HYBRID_NONLINEAR_ELASTICITY<NONLINEAR_ELASTICITY<double,3> >;
template class HYBRID_NONLINEAR_ELASTICITY<SKINNING_NONLINEAR_ELASTICITY<double,3,true,true> >;
template class HYBRID_NONLINEAR_ELASTICITY<SKINNING_NONLINEAR_ELASTICITY<double,3,true,false> >;
#endif

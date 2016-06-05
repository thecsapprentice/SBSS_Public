#include "OVERRIDES.h"


#include <PhysBAM_Tools/Arrays_Computations/ARRAY_COPY.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Data_Structures/QUEUE.h>
#include <PhysBAM_Tools/Data_Structures/HASHTABLE_ITERATOR.h>

#include "REGION_GENERATOR.h"
#include "RANGE_ITERATOR.h"
#include <Thread_Queueing/PTHREAD_QUEUE.h>

using namespace PhysBAM;
extern PTHREAD_QUEUE* pthread_queue;

namespace{
    template<class T, int d> struct EquivalenceClass_Thread_Helper:public PTHREAD_QUEUE::TASK
    {
        REGION_GENERATOR<T,d>* const obj;
        UNION_FIND<int>* rc;
        UNION_FIND<int>* ec;
        const int cell_start,cell_end;
        EquivalenceClass_Thread_Helper(REGION_GENERATOR<T,d>* const obj_input,const int cell_start_input,const int cell_end_input, UNION_FIND<int>* rc_input, UNION_FIND<int>* ec_input )
            :obj(obj_input),cell_start(cell_start_input),cell_end(cell_end_input), rc(rc_input), ec(ec_input) {}
                void Run(){
                    obj->Generate_EquivalenceClass_Range(cell_start,cell_end, *ec, *rc);
                }
    };
}


std::ostream& PhysBAM::operator<<(std::ostream& out, const T_CELLTYPE<float, 2>& cell)
{
    out << "| " << cell.indices[1] << " " << cell.indices[3] << " |"<<std::endl;
    out << "| " << cell.indices[0] << " " << cell.indices[2] << " |"<<std::endl;
    return out;
}
std::ostream& PhysBAM::operator<<(std::ostream& out, const T_CELLTYPE<double, 2>& cell)
{
    out << "| " << cell.indices[1] << " " << cell.indices[3] << " |"<<std::endl;
    out << "| " << cell.indices[0] << " " << cell.indices[2] << " |"<<std::endl;
    return out;
}
std::ostream& PhysBAM::operator<<(std::ostream& out, const T_CELLTYPE<float,3>& cell)
{
    out << "| " << cell.indices[2] << " " << cell.indices[3] << " |"<<std::endl;
    out << "| " << cell.indices[6] << " " << cell.indices[7] << " |"<<std::endl;
    out << std::endl;
    out << "| " << cell.indices[0] << " " << cell.indices[1] << " |"<<std::endl;
    out << "| " << cell.indices[4] << " " << cell.indices[5] << " |"<<std::endl;
    return out;
}
std::ostream& PhysBAM::operator<<(std::ostream& out, const T_CELLTYPE<double,3>& cell)
{
    out << "| " << cell.indices[2] << " " << cell.indices[3] << " |"<<std::endl;
    out << "| " << cell.indices[6] << " " << cell.indices[7] << " |"<<std::endl;
    out << std::endl;
    out << "| " << cell.indices[0] << " " << cell.indices[1] << " |"<<std::endl;
    out << "| " << cell.indices[4] << " " << cell.indices[5] << " |"<<std::endl;
    return out;
}


template<class T, int d>
REGION_GENERATOR<T,d>::REGION_GENERATOR(const T_GRID& grid_input) : grid(grid_input)
{
    
}

template<class T, int d> void
REGION_GENERATOR<T,d>::Generate()
{
    LOG::SCOPE scope("Generating Domain Cuts");


/* 
   PUESDO-CODE

   A) Generate individual root cells from grid. Assign unqiue vertices
*/

    
    CreateRootCells();


/*
   B) For each root cell, call ResolveCuts to break it into X sub-cells, each with
      a side of a cut that passes through the cell. Each sub-cell will have unique
      vertices.
*/
    rootnode_to_subnode.Resize( root_grid_node_mapping.m );
    ResolveCuts();


/*
   C) For each pair of root-cell and its neighbors, compare all sub-cells between
      them with IsMaterialContinous. If two sub-cells are connected, via this
      function, add their shared vertices to a growing Union.
*/
    equivalence_sets.Initialize( sub_cells.m * T_CELL::vertices_per_cell );
    region_sets.Initialize( sub_cells.m );

    RANGE<T_INDEX> grid_domain( grid.Cell_Indices());   
    
    {
        LOG::SCOPE scope("Generating Equivalence Classes for SubCells.");
#if 1
    const int partitions = 32;
    int partition_offsets[partitions];

    UNION_FIND<int> partitioned_ec[partitions];
    UNION_FIND<int> partitioned_rc[partitions];

    int partition;
    for( int i=1, partition=0; i <= root_cells.m && partition < partitions;
         i+=root_cells.m / partitions, partition++ )
        partition_offsets[partition] = i;

    for(partition=0;partition<partitions;partition++){
        int cell_begin=partition_offsets[partition];
        int cell_end=((partition<partitions-1)?partition_offsets[partition+1]:root_cells.m);


        EquivalenceClass_Thread_Helper<T,d>* task = 
            new EquivalenceClass_Thread_Helper<T,d>(this, cell_begin, cell_end, partitioned_rc+partition,
                                                    partitioned_ec+partition);

        pthread_queue->Queue(task);
    }

    pthread_queue->Wait();

    for(partition=0;partition<partitions;partition++){
        region_sets.Merge( partitioned_rc[partition] );
        equivalence_sets.Merge( partitioned_ec[partition] );
    }

#else
    int last_progress = 0;

    for( int root_cell = 1; root_cell <= root_cells.m; root_cell++ )
        {
            if( (int)(((T)(root_cell) / root_cells.m) * 100) % 10 == 0  &&
                (int)(((T)(root_cell) / root_cells.m) * 100) != last_progress){
                last_progress = (int)(((T)(root_cell) / root_cells.m) * 100);
                LOG::cout << "Progress: "<< (int)(((T)(root_cell) / root_cells.m) * 100) << "%" << std::endl;
            }

            const T_INDEX& index = root_grid_cell_mapping( root_cell );

            for( int i = 1; i <= d; i++)
                {
                    T_INDEX offset;
                    offset(i) = 1;
                    T_INDEX neighbor_index = index + offset;
                    
                    if( grid_domain.Lazy_Outside(neighbor_index))
                        continue;

                    int root_cell_neighbor = root_index_to_root_id.Get( neighbor_index );

                    int last_sub_cell = 0;
                    for( int t=1; t <= root_cell_to_subcell(root_cell).m; t++)
                        {
                            last_sub_cell=root_cell_to_subcell(root_cell)(t);
                            int last_sub_cell_neighbor = 0;
                            for( int r=1; r <= root_cell_to_subcell(root_cell_neighbor).m ; r++)
                                {
                                    last_sub_cell_neighbor=root_cell_to_subcell(root_cell_neighbor)(r);
                                    const T_CELL& cellA = sub_cells( last_sub_cell );
                                    const T_CELL& cellB = sub_cells( last_sub_cell_neighbor );

                                       
                                    AXIS faceA_axis;
                                    AXIS faceB_axis;

                                    switch(i){
                                    case 1:
                                        faceA_axis = A_X;
                                        faceB_axis = A_mX;
                                        break;
                                    case 2:
                                        faceA_axis = A_Y;
                                        faceB_axis = A_mY;
                                        break;
                                    case 3:
                                        faceA_axis = A_Z;
                                        faceB_axis = A_mZ;
                                        break; 
                                    }
                                    
                                    bool isContinous = IsMaterialContinous(root_cell,
                                                                           last_sub_cell,
                                                                           faceA_axis,
                                                                           root_cell_neighbor,
                                                                           last_sub_cell_neighbor,
                                                                           faceB_axis );
                                    

                                    if( isContinous ){
                                        T_FACE faceA;
                                        T_FACE faceB;
                                        cellA.GetAxisFace(faceA, faceA_axis);
                                        cellB.GetAxisFace(faceB, faceB_axis);
                                        for( int j=0; j < T_FACE::vertices_per_face; j++){
                                            equivalence_sets.Union( faceA.indices[j], faceB.indices[j] );
                                        }
                                        region_sets.Union( last_sub_cell,  last_sub_cell_neighbor );
                                    }
                                }            
                        }

                }

        }
    
#endif
    }


/*
   D) Relabel sub-cells using the results of the Union-Find to create new cells
      with a minimal set of connected vertices.
*/
    regions.regions.Resize( 0 );
    regions.grid_regions.Resize( 0 );
    regions.max_index = 0;
    regions.vertices.Resize( 0 );

    int max_index=0;

    {
        LOG::SCOPE scope("Relabeling cells to form connected meshes.");


    for( int i = 1; i <= sub_cells.m; i++)
        {
            T_CELL& cell = sub_cells(i);
            int root_cell = root_sub_mapping( i );
            const T_CELL& root = root_cells( root_cell );

            for( int j = 0; j < T_CELL::vertices_per_cell; j++){
                cell.indices[j] = equivalence_sets.Find( cell.indices[j] );
                max_index = max( max_index, cell.indices[j] );

                if( node_to_subcell.m < cell.indices[j])
                    node_to_subcell.Resize( cell.indices[j], true, true, ARRAY<int>() ); 

                if( regions.vertices.m < cell.indices[j])
                    regions.vertices.Resize( cell.indices[j], true, true, TV() );

                //LOG::cout << root_grid_node_mapping(root.indices[j]) << std::endl;
                //LOG::cout << i << std::endl;
                //LOG::cout << cell << std::endl;
                //LOG::cout << std::endl;

                node_to_subcell( cell.indices[j] ).Append( i );
                node_to_subcell( cell.indices[j] ).Prune_Duplicates();
                regions.vertices(cell.indices[j]) = grid.Node(root_grid_node_mapping(root.indices[j]));
            }          
        }
    }
    regions.max_index = max_index;

    {
        LOG::SCOPE scope("Building Incident Lists.");
        BuildIncidentLists();
    }

/*
   E) Organize cells into sets of regions and return the region array.
*/

    // Don't have a good thought on how to do this efficently yet. 
    // For now, simply group all sub_cells together and let thier
    // shared, or not, vertices define the regions.


    ARRAY<int> region_parents;
    {
        LOG::SCOPE scope("Forming all Regions.");

    for( int s = 1; s <= sub_cells.m; s++ )
        {
            const T_INDEX& grid_cell = root_grid_cell_mapping(root_sub_mapping( s ));
            int parent = region_sets.Find( s );
            int insert_region;

            if( region_parents.Find( parent, insert_region) ){
                regions.regions(insert_region).Append( sub_cells(s) );
                regions.grid_regions(insert_region).Append(grid_cell);
            }
            else{
                region_parents.Append(parent);
                regions.regions.Append( typename REGIONS<T,d>::T_REGION() );
                regions.grid_regions.Append( typename REGIONS<T,d>::T_GRID_REGION() );
                regions.regions.Last().Append( sub_cells(s) );
                regions.grid_regions.Last().Append(grid_cell);
            }
        }

    }
    LOG::cout << regions.regions.m << " Regions have been found." << std::endl;

}



template<class T, int d> void
REGION_GENERATOR<T,d>::Generate_EquivalenceClass_Range(int cell_start, int cell_end, 
                                                       UNION_FIND<int>& equivalence_sets,
                                                       UNION_FIND<int>& region_sets)
{

    equivalence_sets.Initialize( sub_cells.m * T_CELL::vertices_per_cell );
    region_sets.Initialize( sub_cells.m );
    
    RANGE<T_INDEX> grid_domain( grid.Cell_Indices());   
    
    {                
        for( int root_cell = cell_start; root_cell <= cell_end; root_cell++ )
            {               
                const T_INDEX& index = root_grid_cell_mapping( root_cell );
                
                for( int i = 1; i <= d; i++)
                    {
                        T_INDEX offset;
                        offset(i) = 1;
                        T_INDEX neighbor_index = index + offset;
                        
                        if( grid_domain.Lazy_Outside(neighbor_index))
                            continue;
                        
                        int root_cell_neighbor = root_index_to_root_id.Get( neighbor_index );
                        
                        int last_sub_cell = 0;
                        for( int t=1; t <= root_cell_to_subcell(root_cell).m; t++)
                            {
                                last_sub_cell=root_cell_to_subcell(root_cell)(t);
                                int last_sub_cell_neighbor = 0;
                                for( int r=1; r <= root_cell_to_subcell(root_cell_neighbor).m ; r++)
                                    {
                                        last_sub_cell_neighbor=root_cell_to_subcell(root_cell_neighbor)(r);
                                        const T_CELL& cellA = sub_cells( last_sub_cell );
                                        const T_CELL& cellB = sub_cells( last_sub_cell_neighbor );
                                        
                                        
                                        AXIS faceA_axis;
                                        AXIS faceB_axis;
                                        
                                        switch(i){
                                        case 1:
                                            faceA_axis = A_X;
                                            faceB_axis = A_mX;
                                            break;
                                        case 2:
                                            faceA_axis = A_Y;
                                            faceB_axis = A_mY;
                                            break;
                                        case 3:
                                            faceA_axis = A_Z;
                                            faceB_axis = A_mZ;
                                            break; 
                                        }
                                        
                                        bool isContinous = IsMaterialContinous(root_cell,
                                                                               last_sub_cell,
                                                                               faceA_axis,
                                                                               root_cell_neighbor,
                                                                               last_sub_cell_neighbor,
                                                                               faceB_axis );
                                        
                                        
                                        if( isContinous ){
                                            T_FACE faceA;
                                            T_FACE faceB;
                                            cellA.GetAxisFace(faceA, faceA_axis);
                                            cellB.GetAxisFace(faceB, faceB_axis);
                                            for( int j=0; j < T_FACE::vertices_per_face; j++){
                                                equivalence_sets.Union( faceA.indices[j], faceB.indices[j] );
                                            }
                                            region_sets.Union( last_sub_cell,  last_sub_cell_neighbor );
                                        }
                                    }            
                            }
                        
                    }
                
            }
    }
}


template<class T, int d> const REGIONS<T,d>* 
REGION_GENERATOR<T,d>::GetRegionData() const
{
    return &regions;
}

template<class T, int d> int
REGION_GENERATOR<T,d>::DuplicatesAtCoarseIndex(const T_INDEX& cell_index) const
{
    int root_cell_index = root_index_to_root_id.Get( cell_index );
    return root_cell_to_subcell( root_cell_index ).m;
}

template<class T,int d> bool
REGION_GENERATOR<T,d>::IsIncident( const T_INDEX& cell_index, const int subcell, const T_INDEX& test_index )
{
    int root_cell_index = root_index_to_root_id.Get( cell_index );
    int sub_cell_index = root_cell_to_subcell( root_cell_index )( subcell );

    int otherroot_cell_index = root_index_to_root_id.Get( test_index );
    for( int test_subcells = 1; test_subcells<=root_cell_to_subcell( otherroot_cell_index ).m; test_subcells++){
        int test_sub_cell_index = root_cell_to_subcell( root_cell_index )( test_subcells );
        if( subcell_incident_list( sub_cell_index ).Count_Matches( test_sub_cell_index ) )
            return true;
    }
    return false;    
}

template<class T, int d> bool
REGION_GENERATOR<T,d>::IsMeshMappable( const T_INDEX& cell_index) const
{
    int root_cell_index = root_index_to_root_id.Get( cell_index );
 
    for( RANGE_ITERATOR<d> iter(RANGE<T_INDEX>(cell_index,cell_index+1)); iter.Valid(); iter.Next()){
        int last_group=-1;
        for( int r=1; r<=node_root_grid_mapping(iter.Index()).m; r++){
            int node_index = node_root_grid_mapping(iter.Index())(r);
            for( int j=1; j <= rootnode_to_subnode( node_index ).m; j++){
                int pre_mapped_node = rootnode_to_subnode( node_index )(j);
                int mapped_node = equivalence_sets.Find( pre_mapped_node );
                if( last_group == -1 || last_group == mapped_node )
                    last_group = mapped_node;
                else
                    return true;
            }
        }
    }
    return false;
}

template<class T, int d> VECTOR<int,T_CELLTYPE<T,d>::vertices_per_cell> 
REGION_GENERATOR<T,d>::GetCellVertices( const T_INDEX& cell_index, int duplicate) const
{
    int root_cell_index = root_index_to_root_id.Get( cell_index );
    PHYSBAM_ASSERT( root_cell_to_subcell( root_cell_index ).m >= duplicate );
    int sub_cell_index = root_cell_to_subcell( root_cell_index )( duplicate );
    
    VECTOR<int,T_CELL::vertices_per_cell> vertices;
    for( int i=1; i<=T_CELL::vertices_per_cell; i++)
        vertices(i) = sub_cells( sub_cell_index ).indices[i-1];

    return vertices;

}

template<class T, int d> void
REGION_GENERATOR<T,d>::CreateRootCells( )
{
    LOG::SCOPE scope("Create Root Cells");

    root_grid_cell_mapping.Resize(0);
    root_grid_node_mapping.Resize(0);
    node_root_grid_mapping.Resize( grid.Node_Indices() );

    int flat_index = 1;

    for(RANGE_ITERATOR<d> iterator(grid.Cell_Indices());iterator.Valid();iterator.Next())
        {
            const T_INDEX& cell_index = iterator.Index();
            
            T_CELL cell;
            int vertex_number = 0;
            for(RANGE_ITERATOR<d> node_iterator(RANGE<T_INDEX>(cell_index, cell_index+1));node_iterator.Valid();node_iterator.Next(), vertex_number++)
                {
                    const T_INDEX& node_index = node_iterator.Index();
                    PHYSBAM_ASSERT( vertex_number < T_CELL::vertices_per_cell );
                    
                    cell.indices[vertex_number] = flat_index;
                    root_grid_node_mapping.Append( node_index );
                    node_root_grid_mapping( node_index ).Append( flat_index );
                    PHYSBAM_ASSERT( flat_index == root_grid_node_mapping.m );
                    flat_index++;                                                   
                }
            root_cells.Append( cell );
            root_grid_cell_mapping.Append( cell_index ); 
            root_index_to_root_id.Insert( cell_index, root_grid_cell_mapping.m );
        }

    PHYSBAM_ASSERT( root_cells.m == root_grid_cell_mapping.m );

}

template<class T, int d> void
REGION_GENERATOR<T,d>::BuildIncidentLists()
{
    subcell_incident_list.Resize( 0 );
    QUEUE< int > touching_subcells( 1 << d );
    // This is really inefficent. Need a better approach.
    for( int scell=1; scell <= sub_cells.m; scell++){
        T_CELL& subcell = sub_cells(scell);
        for( int j = 0; j < T_CELL::vertices_per_cell; j++){ 
            for( int scell_touching=1; scell_touching<=node_to_subcell(subcell.indices[j]).m; scell_touching++){
                touching_subcells.Safe_Enqueue( scell_touching );
            }
        }
        subcell_incident_list.Append(ARRAY<int>());
        while(!touching_subcells.Empty())
            subcell_incident_list( scell ).Append( touching_subcells.Dequeue() );
    }
}


template class REGION_GENERATOR<float, 2>;
template class REGION_GENERATOR<float, 3>;
template class REGION_GENERATOR<double, 2>;
template class REGION_GENERATOR<double, 3>;


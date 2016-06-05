#include "CUTTER_STRATEGY.h"
#include <PhysBAM_Tools/Data_Structures/QUEUE.h>
#include <PhysBAM_Tools/Data_Structures/OPERATION_HASH.h>

using namespace PhysBAM;

template<class T, int d> void
CUTTER_STRATEGY<T,d>::Reset( ){
    subcells->Clean_Memory();
    linkage_list->Clean_Memory();
    subcell_to_root_cell->Clean_Memory();
    root_cell_to_subcell->Clean_Memory();
    subcell_incident_list->Clean_Memory();
    subcell_face_neighbors->Clean_Memory();
    node_to_subcell->Clean_Memory();
    equivalence_sets.Clear_Connectivity();
    restartRequired = true;
}



template<class T, int d> void
CUTTER_STRATEGY<T,d>::Generate_Subcells(){
    mp_instance->MaterialFragments(*mesh,*linkage_list,*subcells,*subcell_to_root_cell);

    root_cell_to_subcell->Clean_Memory();
    root_cell_to_subcell->Resize(mesh->elements.m);
    for( HASHTABLE_ITERATOR<int,int> iterator(*subcell_to_root_cell); iterator.Valid(); iterator.Next()){
        (*root_cell_to_subcell)( iterator.Data() ).Append( iterator.Key() );
    }

    equivalence_sets.Clear_Connectivity();
    equivalence_sets.Initialize( (*subcells).m * T_MESH::vertices_per_cell );
}

template<class T, int d> void
CUTTER_STRATEGY<T,d>::Renumber_And_Build_Mesh(  ) {
    LOG::SCOPE scope("Relabeling cells to form connected meshes.");

    int max_index=0;
    HASHTABLE< ARRAY<int>, int> unique_cells;
    UNION_FIND<int> collapsed_cells;
    collapsed_cells.Initialize( (*subcells).m );

    for( int i = 1; i <= (*subcells).m; i++){
        
        T_CELL& cell = (*subcells)(i);
        int root_cell = subcell_to_root_cell->Get(i);
        //const T_CELL& root = mesh.elements(root_cell);
        
        // Remap
        for( int j = 1; j <= T_MESH::vertices_per_cell; j++){
            cell.vertices(j) = equivalence_sets.Find( cell.vertices(j) );
        }
        // Check
        if( unique_cells.Contains( ARRAY<int>(cell.vertices) ) ){
            //LOG::cout << "SubCell " << i << " has collapsed into subcell " << unique_cells.Get(  ARRAY<int>(cell.vertices)) << std::endl;
            collapsed_cells.Union( i, unique_cells.Get(  ARRAY<int>(cell.vertices)) );
        }
        else{
            //LOG::cout << "SubCell " << i << " final vertices " << cell.vertices << std::endl;
            unique_cells.Set(  ARRAY<int>(cell.vertices), i );
            
            // Process
            for( int j = 1; j <= T_MESH::vertices_per_cell; j++){
                max_index = max( max_index, cell.vertices(j) );
                
                if( (*node_to_subcell).m < cell.vertices(j))
                    (*node_to_subcell).Resize( cell.vertices(j), true, true, ARRAY<int>() ); 
                
                //if( regions.vertices.m < cell.vertices(j))
                //    regions.vertices.Resize( cell.vertices(j), true, true, TV() );
                
                //LOG::cout << root_grid_node_mapping(root.indices[j]) << std::endl;
                //LOG::cout << i << std::endl;
                //LOG::cout << cell << std::endl;
                //LOG::cout << std::endl;
                
                (*node_to_subcell)( cell.vertices(j) ).Append( i );
                (*node_to_subcell)( cell.vertices(j) ).Prune_Duplicates();
                //regions.vertices(cell.vertices(j)) = mesh.Node(root_cell,j);
            }          
        }
        
    }

    {
        LOG::SCOPE scope("Building Incident Lists.");
        BuildIncidentLists();        
    }

    {
        LOG::SCOPE scope("Collapsing material fragments.");
        mp_instance->MergeSubcellMaterial(collapsed_cells);
    }
    
}

template<class T, int d> void
CUTTER_STRATEGY<T,d>::BuildIncidentLists( ) {

    {
        LOG::SCOPE scope( "Building Incident List..." );
        (*subcell_incident_list).Resize( 0 );
        QUEUE< int > touching_subcells( 1 << d );
        // This is really inefficent. Need a better approach.
        for( int scell=1; scell <= (*subcells).m; scell++){
            T_CELL& subcell = (*subcells)(scell);
            for( int j = 1; j <= T_MESH::vertices_per_cell; j++){ 
                for( int scell_touching=1; scell_touching<=(*node_to_subcell)(subcell.vertices(j)).m; scell_touching++){
                    touching_subcells.Safe_Enqueue( (*node_to_subcell)(subcell.vertices(j))(scell_touching) );
                }
            }
            (*subcell_incident_list).Append(ARRAY<int>());
            while(!touching_subcells.Empty())
                (*subcell_incident_list)( scell ).Append( touching_subcells.Dequeue() );
        }
    }

    {
        LOG::SCOPE scope( "Building Face Neighbor List..." );
        OPERATION_HASH<int> canidates;
        canidates.Initialize( (*subcells).m ); 
        (*subcell_face_neighbors).Resize( (*subcells).m );
        for( int scell = 1; scell <= (*subcells).m; scell++ ){
            T_CELL& subcell = (*subcells)(scell);
            for( int f = 1; f <= T_MESH::faces_per_cell; f++){
                int axis = ((f-1) / 2)+1;
                bool direction = bool( ((f-1) % 2) );
                VECTOR<int,T_MESH::vertices_per_face> face = subcell.face( axis, direction );
                
                ARRAY<int> &final_pass_faces = (*subcell_face_neighbors)(scell)(f);
                // Progressively mark subcells that have subsequent vertices AND prior ones
                for( int v = 1; v <= T_MESH::vertices_per_face; v++){
                    canidates.Next_Operation();
                    if( v == 1 ){
                        for( int scell_touching=1; scell_touching <= (*node_to_subcell)(face(1)).m; scell_touching++ ){
                            if( (*node_to_subcell)(face(1))(scell_touching) == scell )
                                continue;
                            canidates.Mark( (*node_to_subcell)(face(1))(scell_touching) );
                        }  
                    }
                    else if( v != T_MESH::vertices_per_face) {
                        for( int scell_touching=1; scell_touching<=(*node_to_subcell)(face(v)).m; scell_touching++ ){
                            if(canidates.operations( (*node_to_subcell)(face(v))(scell_touching) ) == canidates.current_operation-1)
                                canidates.Mark( (*node_to_subcell)(face(v))(scell_touching) );
                            
                        }
                    }
                    else{
                        for( int scell_touching=1; scell_touching<=(*node_to_subcell)(face(v)).m; scell_touching++ ){
                            if(canidates.operations( (*node_to_subcell)(face(v))(scell_touching) ) == canidates.current_operation-1)
                                final_pass_faces.Append( (*node_to_subcell)(face(v))(scell_touching) );
                            
                        }
                    }
                }           
            }
        }
    }

}


template<class T, int d> void
CUTTER_STRATEGY<T,d>::Collect_Regions( ) {

/*
    ARRAY<int> region_parents;
    {
        LOG::SCOPE scope("Forming all Regions.");
        
        for( int s = 1; s <= (*subcells).m; s++ )
            {
                //const int rootcell = subcell_to_root_cell.Get(s);
                //const T_INDEX& grid_cell = mesh.elements(rootcell).index;
                
                int parent = region_sets.Find( s );
                int insert_region;
                
                if( region_parents.Find( parent, insert_region) ){
                    regions.regions(insert_region).Append( (*subcells)(s) );
                    regions.region_ids(insert_region).Append(s);
                    //regions.grid_regions(insert_region).Append(grid_cell);
                }
                else{
                    region_parents.Append(parent);
                    regions.regions.Append( typename CUTTER_REGIONS<T,d>::T_REGION() );
                    regions.region_ids.Append( typename CUTTER_REGIONS<T,d>::T_REGION_IDS() );
                    //regions.grid_regions.Append( typename CUTTER_REGIONS<T,d>::T_GRID_REGION() );
                    regions.regions.Last().Append( (*subcells)(s) );
                    regions.region_ids.Last().Append( s );
                    //regions.grid_regions.Last().Append(grid_cell);
                }
            }
        
    }
*/  

}


template class CUTTER_STRATEGY<float,2>;
template class CUTTER_STRATEGY<float,3>;
template class CUTTER_STRATEGY<double,2>;
template class CUTTER_STRATEGY<double,3>;

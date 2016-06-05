#include "NONMANIFOLD_LEVELSET_MESH.h"

// #include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR.h>
// #include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_NODE.h>
// #include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_EDGE.h>

// #include <PhysBAM_Tools/Data_Structures/HASHTABLE.h>
// #include <PhysBAM_Tools/Data_Structures/HASHTABLE_ITERATOR.h>
// #include <PhysBAM_Tools/Data_Structures/OPERATION_HASH.h>

namespace PhysBAM{

    // template< class T, int d>
    // NONMANIFOLD_LEVELSET_MESH<T,d>::NONMANIFOLD_LEVELSET_MESH(T_GRID& grid){
    //     Initialize_From_Grid(grid);
    // }
        
    // template< class T, int d>
    // NONMANIFOLD_LEVELSET_MESH<T,d>::~NONMANIFOLD_LEVELSET_MESH<T,d>(){
    // }
    
    
    // template<class T, int d>
    // NONMANIFOLD_LEVELSET_MESH<T,d>::Initialize_From_Grid( T_GRID& grid ){
    //     nodes.Clear();
    //     edges.Clear();
        
    //     grid_base = grid;
        
    //     HASHTABLE<T_INDEX, int> node_map;

    //     for( UNIFORM_GRID_ITERATOR_NODE<TV> iterator(grid_base); iterator.Valid(); iterator.Next()){
    //         const T_INDEX& node_index = iterator.Node_Index();
    //         T_NODE nm_node;
    //         nm_node.grid_index = node_index;
    //         int p = nodes.Append( nm_node );
    //         node_map.Get_Or_Insert( node_index, p );            
    //     }

    //     RANGE<T_INDEX> node_domain = grid_base.Node_Indices();
    //     for( int iterator = 1; iterator <= nodes.m; iterator++){
    //         T_NODE& node = nodes(iterator);
    //         const T_INDEX& node_index = node.grid_index;
            
    //         for( int i = 1; i <= grid_base.number_of_neighbors_per_node; i++){
    //             const T_INDEX& index = grid_base.Node_Neighbor(node_index, i);
    //             if(node_domain.Lazy_Inside(index)){
    //                 int neighbor = node_map.Get(index);
    //                 node.neighbors.Append(neighbor);
    //             }
    //         }
    //     }

    //     for( UNIFORM_GRID_ITERATOR_EDGE<TV> iterator(grid_base); iterator.Valid(); iterator.Next()){
    //         const EDGE_INDEX<d>& edge_index = iterator.Full_Index();
    //         const T_INDEX nodeA_index = edge_index.First_Node_Index();
    //         const T_INDEX nodeB_index = edge_index.Second_Node_Index();

    //         int nodeA_id = node_map.Get(nodeA_index);
    //         int nodeB_id = node_map.Get(nodeB_index);

    //         T_EDGE nm_edge;
    //         nm_edge.first_node = nodeA_id;
    //         nm_edge.second_node = nodeB_id;

    //         int p = edges.Append(nm_edge);

    //         T_NODE& nodeA = nodes(nodeA_id);
    //         T_NODE& nodeB = nodes(nodeB_id);
    //         nodeA.edges.Append_Unique(p);
    //         nodeB.edges.Append_Unique(p);
    //     }

    //     OPERATION_HASH<int,int> edge_set;
    //     edge_set.Initialize( edges.m );
    //     for( int iterator = 1; iterator <= edges.m; iterator++){
    //         T_NODE& edge = edges(iterator);
    //         edge_set.Next_Operation();

    //         for( int j = 1; j <= nodes(edge.first_node).edges.m; j++)
    //             if(nodes(edge.first_node).edges(j) != iterator )
    //                 edge.neighbors.Append_Unique( nodes(edge.first_node).edges(j)

    //         for( int j = 1; j <= nodes(edge.second_node).edges.m; j++)
    //             if(nodes(edge.second_node).edges(j) != iterator )
    //                 edge.neighbors.Append_Unique( nodes(edge.second_node).edges(j)

    //     }
    

    // }

};

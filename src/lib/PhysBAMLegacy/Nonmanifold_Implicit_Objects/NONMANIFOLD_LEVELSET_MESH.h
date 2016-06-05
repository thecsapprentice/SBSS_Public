#ifndef _NONMANIFOLD_LEVELSET_MESH_H_
#define _NONMANIFOLD_LEVELSET_MESH_H_

#include <PhysBAM_Tools/Data_Structures/PAIR.h>
#include <PhysBAM_Tools/Data_Structures/TRIPLE.h>
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <Common_Tools/Math_Tools/RANGE_ITERATOR.h>
#include <PhysBAM_Geometry/Spatial_Acceleration/BOX_HIERARCHY.h>
#include <PhysBAM_Tools/Data_Structures/OPERATION_HASH.h>
#include <PhysBAM_Tools/Data_Structures/HASHTABLE.h>
#include <PhysBAM_Tools/Data_Structures/HASHTABLE_ITERATOR.h>
#include <Common_Tools/Arrays/HYBRID_ARRAY.h>
#include <Common_Tools/Arrays/HYBRID_ARRAY_ITERATOR.h>
#include <PhysBAM_Tools/Math_Tools/Hash.h>
#include <stdint.h>

namespace PhysBAM{

    template <int d> class NONMANIFOLD_LEVELSET_MESH_POLICY;

    template <class T, int d>        
        struct NONMANIFOLD_LEVELSET_MESH_CELL{
            enum { nodes_per_cell=1<<d,
                   nodes_per_face=1<<(d-1),
                   faces_per_cell=2*d
            };
            
            typedef T SCALAR;
            typedef VECTOR<T,d> TV;
            typedef VECTOR<int,d> T_INDEX;
            typedef VECTOR<int,nodes_per_face> T_FACE;
            typedef PAIR<T_INDEX, int> HYBRID_INDEX;
            typedef HYBRID_INDEX T_NODE;
            
        T_INDEX cell_index;
        VECTOR<T_NODE,nodes_per_cell> nodes;
        // PAIR<index,array_index>
        /* index        == 0 if cell does not have a transition face, otherwise points to entry into the transition_faces   */
        /* array_index  == points to entry in the array left_cells/right_cells                                              */
        VECTOR<PAIR<int,int>,2*d> transition;
        // ARRAY<int>
        /* All of the face neighbors which are NOT part of a transition   
         */
        VECTOR<ARRAY<int>,2*d> cell_neighbors;
        
        static T_FACE face_indices(int axis, bool direction) {
            VECTOR<int,nodes_per_face> f;
            int flat_index=1;
            int face_index=1;
            for(RANGE_ITERATOR<d> iter(RANGE<T_INDEX>(T_INDEX(),T_INDEX()+1)); iter.Valid(); iter.Next()){
                const T_INDEX& index=iter.Index();
                if( index(axis) == direction ){
                    f(face_index) = flat_index;
                    face_index++;
                }
                flat_index++;
            }
            PHYSBAM_ASSERT(face_index == nodes_per_face+1);
            return f;
        };
        
        T_FACE face(int axis, bool direction) const {
            VECTOR<int,nodes_per_face> f;
            int flat_index=1;
            int face_index=1;
            for(RANGE_ITERATOR<d> iter(RANGE<T_INDEX>(T_INDEX(),T_INDEX()+1)); iter.Valid(); iter.Next()){
                const T_INDEX& index=iter.Index();
                if( index(axis) == direction ){
                    f(face_index) = nodes(flat_index);
                    face_index++;
                }
                flat_index++;
            }
            PHYSBAM_ASSERT(face_index == nodes_per_face+1);
            return f;
        };
        
        VECTOR<T_FACE,faces_per_cell> faces() const {
            VECTOR<T_FACE,faces_per_cell> f;
            for( int i=0; i<d; i++)
                for( int a=0; a<2; a++)
                    f( (2*i + a) +1 ) = face(i+1, a);
            return f;
        };
    };

    template <class T, int d>
    struct NONMANIFOLD_LEVELSET_MESH_TRANSITION_FACE{
            uint32_t type;
        // entries are ordered in the same order as the bit mask in type (lower bits being lower index)
        ARRAY<PAIR<int,int> > left_cells,right_cells;       // PAIR<cell_index,face_index>
    };

    template <class T, int d>
    class NONMANIFOLD_LEVELSET_MESH {
        STATIC_ASSERT((d==2)||(d==3));
    public:
        enum { nodes_per_cell=1<<d,
               nodes_per_face=1<<(d-1),
               faces_per_cell=2*d
        };
        typedef T SCALAR;
        typedef VECTOR<T,d> TV;
        typedef VECTOR<int,d> T_INDEX;
        typedef VECTOR<int,nodes_per_face> T_FACE;
        typedef ARRAY<bool> T_FLAG_ARRAY;

        typedef NONMANIFOLD_LEVELSET_MESH_POLICY<d> T_POLICY;

        typedef PAIR<T_INDEX, int> HYBRID_INDEX;
        typedef HYBRID_INDEX T_NODE;
        //          <     0 , >0 > -- Mesh Mapped
        //          < index , 0  > -- Grid Mapped
        
        typedef ARRAY<T_NODE> NEIGHBOR_ARRAY;
        // typedef VECTOR<NEIGHBOR_ARRAY, 2*d > AXIS_ALIGNED_NEIGHBORS;
        typedef VECTOR<NEIGHBOR_ARRAY,2*d> AXIS_ALIGNED_NEIGHBORS;
        
        typedef NONMANIFOLD_LEVELSET_MESH_CELL<T,d> T_CELL;
        typedef NONMANIFOLD_LEVELSET_MESH_TRANSITION_FACE<T,d> TRANSITION_FACE;
                
        T_FLAG_ARRAY cell_is_grid;
        
        HASHTABLE<T_NODE, AXIS_ALIGNED_NEIGHBORS> node_neighbors;
        AXIS_ALIGNED_NEIGHBORS empty_neighbors;
        
        ARRAY<T_CELL> cells;
        ARRAY<TRANSITION_FACE> transition_faces;
        bool initialized;

        RANGE<T_INDEX> node_domain;
        int node_mesh_count;

        // Location Structures
        T dx;
        BOX_HIERARCHY<TV> aabb_hierarchy;
        HASHTABLE<T_NODE,TV> node_locations;

        NONMANIFOLD_LEVELSET_MESH(T dx) : dx(dx)
        {
            cells.Remove_All();
            cell_is_grid.Resize( 0, true, false);
            cell_is_grid.Fill(false);

            node_domain = RANGE<T_INDEX>(T_INDEX(),T_INDEX());
            node_mesh_count = 0;

            //ARRAY<RANGE<TV> > leaf_boxes;
            //aabb_hierarchy.Set_Leaf_Boxes(leaf_boxes, true);
            aabb_hierarchy.Clean_Memory();
            //Update_Node_Neighbors();
            initialized = false;
        }

        template<class T_GRID>
        NONMANIFOLD_LEVELSET_MESH(T_GRID& grid) : dx(grid.dX(1))
        {
            cells.Remove_All();
            cell_is_grid.Resize( 0, true, false);

            node_domain = grid.Node_Indices();
            node_mesh_count = 0;

            ARRAY<RANGE<TV> > leaf_boxes;

            LOG::cout << (unsigned int)(Hash (dx ) ) << std::endl;
            LOG::cout << (unsigned int)(Hash (grid.dX(1)) ) << std::endl;

            for(RANGE_ITERATOR<d> cell_iterator(grid.Cell_Indices());cell_iterator.Valid();cell_iterator.Next()){
                const T_INDEX& cell_index=cell_iterator.Index();
                T_CELL cell;
                cell.cell_index = cell_index;
                int index=0;
                for(RANGE_ITERATOR<d> node_iterator(RANGE<T_INDEX>(cell_index,cell_index+1));node_iterator.Valid();node_iterator.Next()){
                    cell.nodes(++index).x=node_iterator.Index();
                    cell.nodes(index).y=0;
                    node_locations.Set( cell.nodes(index), grid.X(node_iterator.Index()));
                }
                cells.Append(cell);
                cell_is_grid.Append(true);
                //leaf_boxes.Append( grid.Cell_Domain(cell_index) );
                leaf_boxes.Append( RANGE<TV>( grid.X(cell_index), grid.X(cell_index+1)) );
            }
            aabb_hierarchy.Set_Leaf_Boxes(leaf_boxes, true);
            
            Update_Node_Neighbors();
            initialized = true;
        }

        template<class T_GRID> 
        void Append_Cells_From_Grid(T_GRID& grid, ARRAY<bool, T_INDEX> active_cells, HASHTABLE<T_INDEX,HYBRID_INDEX>* gridnode_map = 0 ){
            ARRAY<RANGE<TV> > new_leaf_boxes;
            for(int i = 1; i <= aabb_hierarchy.leaves; i++)
                new_leaf_boxes.Append( aabb_hierarchy.box_hierarchy(i) );
            
            if( !initialized ){
                node_domain = grid.Node_Indices();
            }
            else{
                RANGE<T_INDEX> new_node_domain = node_domain;
                new_node_domain.Enlarge_To_Include_Box( grid.Node_Indices() );
                node_domain = new_node_domain;
            }

            int new_node_mesh_count = node_mesh_count;
            
            HASHTABLE<T_INDEX,int> range_map;
            VECTOR<T_INDEX, 4> inv_range_map;
            int flat_index = 1;
            for(RANGE_ITERATOR<d> iterator(RANGE<T_INDEX>(T_INDEX(),T_INDEX()+1));
                iterator.Valid();iterator.Next(),flat_index++){
                range_map.Insert(iterator.Index(), flat_index);
                inv_range_map(flat_index) = iterator.Index();
            }

            // Build acceleration structure for mesh/grid mapping.
            HASHTABLE<T_INDEX,ARRAY<int> > cell_bins;
            for( int i = 1; i <= cells.m; i++ ){
                ARRAY<int> bin = cell_bins.Get_Or_Insert(cells(i).cell_index);
                bin.Append( i );
                cell_bins.Set( cells(i).cell_index, bin );
            }


            // Build acceleration structure for mesh/grid mapping.
            HASHTABLE<T_INDEX,int > node_bins;
            for( int i = 1; i <= cells.m; i++ ){
                for( int j = 1; j <= 4; j++){
                    int bin = node_bins.Get_Or_Insert(cells(i).cell_index + inv_range_map(j));
                    bin++;
                    node_bins.Set( cells(i).cell_index + inv_range_map(j), bin );
                }
            }

            // Insert new cells, marking mesh coord. as -1 to indicate unknown assignment thus far
            for(RANGE_ITERATOR<d> cell_iterator(grid.Cell_Indices());cell_iterator.Valid();cell_iterator.Next()){
                const T_INDEX& cell_index=cell_iterator.Index();
                const bool& active_cell = active_cells(cell_iterator.Index());                
                if(active_cell){
                    T_CELL cell;
                    cell.cell_index = cell_index;
                    int index=0;
                    ARRAY<int> bin = cell_bins.Get_Or_Insert( cell_index );
                    LOG::cout << "Building new cell (" << cell_index << "):" << std::endl;
                    for(RANGE_ITERATOR<d> node_iterator(RANGE<T_INDEX>(cell_index,cell_index+1));node_iterator.Valid();node_iterator.Next()){                        
                        cell.nodes(++index).x=node_iterator.Index();
                        cell.nodes(index).y=-1;
                        LOG::cout << "\t(" << index << ") " << node_iterator.Index() << std::endl;
                    }
                    int p = cells.Append(cell);
                    //LOG::cout << "Inserting new active cell "<< cell_index << " : " << p << std::endl;
                    cell_is_grid.Append(false);
                    bin.Append(p);
                    //LOG::cout << "Cell bin count: " << bin.m << std::endl;
                    //leaf_boxes.Append( grid.Cell_Domain(cell_index) );
                    new_leaf_boxes.Append( RANGE<TV>( grid.X(cell_index), grid.X(cell_index+1)) );
                    cell_bins.Set( cell_index, bin );
                }
            }

            for( int i = 1; i <= cells.m; i++){
                const T_INDEX& cell_index = cells(i).cell_index;
                T_CELL& cell = cells(i);
                int index=1;
                LOG::cout << "Processing remaps for cell " << i << std::endl;
                for(RANGE_ITERATOR<d> node_iterator(RANGE<T_INDEX>(cell_index,cell_index+1));node_iterator.Valid();node_iterator.Next(), index++){ 
                    int bin = node_bins.Get_Default( node_iterator.Index(), 0 );
                    HYBRID_INDEX original_node = cell.nodes(index);
                    if( cell.nodes(index).y == -1 && bin > 0 ){
                        cell.nodes(index).y = ++new_node_mesh_count;
                        cell.nodes(index).x = T_INDEX();
                    }

                    if( cell.nodes(index).y == -1 && bin == 0 ){
                        cell.nodes(index).y = 0;
                    }

                    if( bin < 0 )
                        PHYSBAM_FATAL_ERROR();

                    HYBRID_INDEX new_node = cell.nodes(index);
                    node_locations.Set( new_node, grid.X(original_node.index) );

                    if(original_node != new_node) LOG::cout << "Remapping node " << original_node << " to " << new_node << std::endl;
                    if( gridnode_map && original_node != new_node) gridnode_map->Set( original_node.x, new_node );

                    if(original_node != new_node )
                        for( RANGE_ITERATOR<d> cell_iterator(RANGE<T_INDEX>(original_node.x-1, original_node.x));
                             cell_iterator.Valid(); cell_iterator.Next()){
                            ARRAY<int> candidate_cells = cell_bins.Get_Or_Insert( cell_iterator.Index() );
                            for( int c = 1; c <= candidate_cells.m; c++){
                                if( candidate_cells(c) == i ) continue;
                                int candidate_node_index = range_map.Get( original_node.x-cell_iterator.Index());
                                if( cells( candidate_cells(c) ).nodes( candidate_node_index ) == original_node)
                                    cells( candidate_cells(c) ).nodes( candidate_node_index ) = new_node;
                            }
                        }
                }
            }

            node_mesh_count = new_node_mesh_count;
            aabb_hierarchy.Set_Leaf_Boxes(new_leaf_boxes, true);            
            Update_Node_Neighbors();
            initialized = true;
        }
        
        const TV Node( const T_NODE& node ){
            return node_locations.Get_Default( node, TV() );                        
        }


        const AXIS_ALIGNED_NEIGHBORS Node_Neighbors( const T_NODE& node) const {
            return node_neighbors.Get_Default(node, empty_neighbors );
        }

        const RANGE<T_INDEX>& Node_Domain() const {
            return node_domain;
        }

        const int& Node_Mesh_Counts() const {
            return node_mesh_count;
        }

        TV Location(const int cell_index,const TV& weights) const
        {return dx*(TV(cells(cell_index).cell_index-1)+weights);}

        RANGE<TV> Cell_Domain(const int cell_index) const
        {return RANGE<TV>(Location(cell_index,TV()),Location(cell_index,TV::All_Ones_Vector()));}
    
        void Update_Node_Neighbors(){

            VECTOR<T_INDEX,2*d> neighbor_stencil;
            for(int i = 0; i<d; i++)
                for(int j = 1; j<=2; j++)
                    neighbor_stencil(i*2+j) = T_INDEX::Axis_Vector(i+1) * ( j == 1 ? 1 : -1 );

            VECTOR<int,2*d> neighbor_inverse;
            for(int i = 0; i<d; i++)
                for(int j = 1; j<=2; j++)
                    neighbor_inverse(i*2+j) = (i*2) + (3-j);

            HASHTABLE<T_INDEX, int> range_map;
            int flat_index = 1;
            for(RANGE_ITERATOR<d> iterator(RANGE<T_INDEX>(T_INDEX(),T_INDEX()+1));
                iterator.Valid();iterator.Next(),flat_index++)
                range_map.Insert(iterator.Index(), flat_index);
            
            for( int i = 1; i <= cells.m ; i++ ){
                const T_CELL& cell = cells(i);
                int flat_index=1;
                //LOG::cout << "Considering cell " << i << " with geo-index " << cell.cell_index << std::endl;
                RANGE<T_INDEX> node_range(T_INDEX(),T_INDEX()+1);
                for(RANGE_ITERATOR<d> node_iterator(node_range);node_iterator.Valid();node_iterator.Next(),flat_index++){
                    
                    T_INDEX node_index = node_iterator.Index();
                    const HYBRID_INDEX& node_hindex = cells(i).nodes(flat_index);
                    AXIS_ALIGNED_NEIGHBORS neighbors = node_neighbors.Get_Or_Insert( node_hindex );

                    for(int j = 1; j <= (d*2); j++){
                        T_INDEX neighbor_index = node_index + neighbor_stencil(j);
                        if(node_range.Lazy_Inside( neighbor_index ) ){                            
                            HYBRID_INDEX neighbor_hindex = cell.nodes(range_map.Get(neighbor_index));
                            AXIS_ALIGNED_NEIGHBORS neighbors_neighbors = node_neighbors.Get_Or_Insert( neighbor_hindex );
                            //LOG::cout << "Considering connection between <" << node_hindex.x << ","<< node_hindex.y << "> and <" << neighbor_hindex.x << ","<< neighbor_hindex.y << ">" << std::endl;
                            neighbors(j).Append_Unique( neighbor_hindex);
                            neighbors_neighbors( neighbor_inverse(j) ).Append_Unique( node_hindex );
                            node_neighbors.Set( neighbor_hindex, neighbors_neighbors );
                        }
                    }
                    node_neighbors.Set( node_hindex, neighbors );
                }
            }                     

            // Run to check neighbor counts for consistency
            /*
            HYBRID_ARRAY<int, d> neighbor_counts;
            neighbor_counts.Resize(node_domain, node_mesh_count); 

            VECTOR<int,10> counts;
            for(HYBRID_ARRAY_ITERATOR< d > iterator(neighbor_counts); iterator.Valid(); iterator.Next() ){
                const HYBRID_INDEX& node_hindex = iterator.Index();
                if(!node_neighbors.Contains(node_hindex))
                    continue;
                neighbor_counts(node_hindex) = 0;
                for(int j = 1; j <= (d*2); j++){
                    neighbor_counts(node_hindex) += node_neighbors.Get(node_hindex)(j).m;
                }
                counts(neighbor_counts(node_hindex))++;
            }
            //for(int i = 1; i<=10; i++)
            //    LOG::cout << "Number of nodes with " << i << " neighbors: " << counts(i) << std::endl;
            //PHYSBAM_FATAL_ERROR();
            */


        }

        void Compute_Acceleration_Structures_For_Transition_Faces()
        {
            // first initialize all transitions to zero.
            for(int i=1;i<=cells.m;++i) for(int j=1;j<=faces_per_cell;++j) cells(i).transition(j)=PAIR<int,int>(0,0);
            
            // now populate all the transitions.
            for(int i=1;i<=transition_faces.m;++i){const TRANSITION_FACE& t_face=transition_faces(i);
                for(int j=1;j<=t_face.left_cells.m;++j) cells(t_face.left_cells(j).x).transition(t_face.left_cells(j).y)=PAIR<int,int>(i,j);
                for(int j=1;j<=t_face.right_cells.m;++j) cells(t_face.right_cells(j).x).transition(t_face.right_cells(j).y)=PAIR<int,int>(i,j);}
        }
    };

};

#endif

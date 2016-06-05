#include "MATERIAL_PREDICATE_TESSELATED_VOLUME.h"
#include <PhysBAM_Geometry/Spatial_Acceleration/BOX_HIERARCHY.h>
#include <PhysBAM_Tools/Arrays/INDIRECT_ARRAY.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Vectors/VECTOR_2D.h>
#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>
#include <PhysBAM_Tools/Data_Structures/UNION_FIND.h>

#include <Common_Geometry/Nonmanifold_Topology_Generation/GEOMETRY_PREDICATES.h>

using namespace PhysBAM;


//=====================================================================================
//
//                                 Predicate Methods
//
//=====================================================================================

template<class T, int d>
MATERIAL_PREDICATE_TESSELATED_VOLUME<T,d>::MATERIAL_PREDICATE_TESSELATED_VOLUME(T_VOLUME& volume )
: volume(volume) 
{
    volume.mesh.Initialize_Boundary_Nodes();
    boundary = &volume.Get_Boundary_Object();

    ARRAY< RANGE<TV> > vbbs;
    for( int i = 1; i <= volume.mesh.elements.m; i++ ){
        RANGE<TV> vbb( volume.Get_Element(i).Bounding_Box() );
        //LOG::cout << "Volume Element " << i << " has a bounding box of: " << vbb << std::endl;
        vbbs.Append( vbb );
    }
    aabb_hierarchy.Set_Leaf_Boxes( vbbs, true );
}

template<class T, int d> void
MATERIAL_PREDICATE_TESSELATED_VOLUME<T,d>::MaterialFragments(const T_MESH& mesh,
                                                             const HASHTABLE<PAIR<VECTOR<int,2>,int> >& linkage_list,
                                                             ARRAY<T_CELL>& sub_cells, 
                                                             HASHTABLE<int, int>& subcell_to_root_cell){
    int vertex_counter = 1;
    sub_cells.Clean_Memory();
    subcell_to_root_cell.Clean_Memory();
    material_fragments.Clean_Memory();
    corner_fragments.Clean_Memory();

    int last_progress = 0;
#pragma omp parallel for shared(last_progress, mesh, linkage_list, vertex_counter, sub_cells, subcell_to_root_cell)
    for( int r = 1; r <= mesh.elements.m; r++ ){

        //LOG::cout << "Checking root cell: " << r << std::endl;
        ARRAY<int> intersecting_mat_elements;
        // Generate All intersecting fragments

        RANGE<TV> root_box( mesh.Node(mesh.elements(r).index), mesh.Node(mesh.elements(r).index+1) );
        //LOG::cout << "\tRoot box: " << root_box << std::endl;
        ARRAY<int> intersects;

#pragma omp critical
        aabb_hierarchy.Intersection_List(root_box, intersects,0);

        //LOG::cout << "\tPreliminary Volume Intersections: "<< intersects.m << std::endl;
       
        for( int i=1; i <= intersects.m; i++){
            //LOG::cout << "Testing Volume Element: " << intersects(i) << std::endl;
            if( GEOMETRY_PREDICATES<T>::TestIntersection(root_box, volume.Get_Element(intersects(i)))) {
                intersecting_mat_elements.Append( intersects(i) );        
                //LOG::cout << "Intersects " << intersects(i) << std::endl;
            }
        }

        //LOG::cout << "Root cell intersects " << intersecting_mat_elements.m << " material elements." << std::endl;

        // Join all connected fragments into connected components

        UNION_FIND<int> connected_components;
        connected_components.Initialize( intersecting_mat_elements.m );
        for( int i = 1; i <= intersecting_mat_elements.m; i++ )
            for( int j = 1; j <= intersecting_mat_elements.m; j++ ){
                if( i==j ) continue;
                VECTOR<int,d> neighbor_edge;
                bool neighbors = GEOMETRY_PREDICATES<T>::Element_Neighbors( volume.mesh, intersecting_mat_elements(i), intersecting_mat_elements(j), neighbor_edge); 
                int smaller_simplex_index=min(intersecting_mat_elements(i),intersecting_mat_elements(j));
                int larger_simplex_index=max(intersecting_mat_elements(i),intersecting_mat_elements(j));
                bool is_present_in_linkage_list=linkage_list.Contains(PAIR<VECTOR<int,2>,int>(VECTOR<int,2>(smaller_simplex_index,larger_simplex_index),r));   // if true: assumes these triangles intersect the bounds of the cell
                if(is_present_in_linkage_list)
                    connected_components.Union( i, j ); // Then they are connected within the cell;
                else
                    if(neighbors){ // Do these elements touch at all?
                        T_BOUNDARY_ELEMENT nee( volume.particles.X.Subset(neighbor_edge) );
                        if( GEOMETRY_PREDICATES<T>::TestIntersection( root_box, nee ) ) // Do these elements touch within the bounds of the cell?
                            connected_components.Union( i, j ); // Then they are connected within the cell;
                }
            }


        /*
         * NOT DOING THIS ANY MORE AS PER REVISED ALGORITHM!!
         *
        // For each corner, merge all fragment elements which contain it
        {
            VECTOR<ARRAY<int>, T_MESH::vertices_per_cell> cfragments;
            for( int v = 1; v <= mesh.vertices_per_cell; v++){
                for( int f = 1; f <= intersecting_mat_elements.m; f++){
                    const T_ELEMENT element = volume.Get_Element( intersecting_mat_elements(f) );
                    if( element.Inside( mesh.Node(r, v) ) )
                        cfragments(v).Append( f );
                }
                
                if(cfragments(v).m >= 1){
                    int root_ele = cfragments(v)(1);
                    for( int c = 1; c <= cfragments(v).m; c++){
                    connected_components.Union( root_ele, cfragments(v)(c) );
                    }
                }
            }
        }
        */
       
        // Collect components
        
        HASHTABLE<int,ARRAY<int> > components;
        for( int i = 1; i <= intersecting_mat_elements.m; i++ ){
            ARRAY<int> c = components.Get_Or_Insert( connected_components.Find(i), ARRAY<int>() );
            c.Append(intersecting_mat_elements(i));
            components.Set( connected_components.Find(i), c );
        }
       
        // Assemble subcells

#pragma omp critical 
        {
            int num_components = 0;
            for( HASHTABLE_ITERATOR<int,ARRAY<int> > iterator(components); iterator.Valid(); iterator.Next()){
                num_components++;
                T_CELL new_subcell;            
                new_subcell.index = mesh.elements(r).index;
                VECTOR<ARRAY<int>, T_MESH::vertices_per_cell> cfragments;
                for( int v = 1; v <= mesh.vertices_per_cell; v++){
                    new_subcell.vertices(v) = vertex_counter++;
                    for( int f = 1; f <= iterator.Data().m; f++){
                        const T_ELEMENT element = volume.Get_Element( iterator.Data()(f) );
                        if( GEOMETRY_PREDICATES<T>::Inside( element, mesh.Node(r, v) ) )
                            cfragments(v).Append( iterator.Data()(f) );
                    }
                }
                
                //LOG::cout << new_subcell.index << " , " << new_subcell.vertices << std::endl;
                //LOG::cout << "  Fragments: "<< std::endl;
                //for( int f = 1; f <= iterator.Data().m; f++){
                //    LOG::cout << "    " << iterator.Data()(f) << std::endl;
                //}
                
                int s = sub_cells.Append( new_subcell );
                subcell_to_root_cell.Set(s, r);
                material_fragments.Append( iterator.Data() );
                corner_fragments.Append(cfragments);
                
            }
            if( (int)(((T)(last_progress+1) / mesh.elements.m) * 100) % 5 == 0  &&
                (int)(((T)(last_progress) / mesh.elements.m) * 100) % 5 != 0 ){
                LOG::cout << "Progress: "<< (int)(((T)(last_progress+1) / mesh.elements.m) * 100) << "%" << std::endl;
            }
            last_progress++;
        }             
        //LOG::cout << "Generating " << num_components << " subcells from this root cell." << std::endl;
    }

    //Build_Boundary_Fragments( mesh, sub_cells );

}

template<class T, int d> bool
MATERIAL_PREDICATE_TESSELATED_VOLUME<T,d>::IsMaterialContinuous(const int axis,
                                                                const int root_cell_high, const int sub_cell_high,
                                                                const int root_cell_low,  const int sub_cell_low) const{

    //VECTOR<int,T_MESH::vertices_per_face> face_high = T_MESH::T_CELL::face_indices(axis, 0 );
    //VECTOR<int,T_MESH::vertices_per_face> face_low = T_MESH::T_CELL::face_indices(axis, 1 );

    //for( int k=1; k <= T_MESH::vertices_per_face; k++)
    //    for( int i = 1; i <= corner_fragments( sub_cell_high )( face_high(k) ).m; i++)
    //        for( int j = 1; j <= corner_fragments( sub_cell_low )( face_low(k) ).m; j++)
    //            if( corner_fragments(sub_cell_high)(face_high(k))(i) == corner_fragments(sub_cell_low)(face_low(k))(j) )
    //                return true;

    for(int i=1;i<=material_fragments(sub_cell_high).m;++i)
        for(int j=1;j<=material_fragments(sub_cell_low).m;++j)
            if(material_fragments(sub_cell_high)(i)==material_fragments(sub_cell_low)(j)) return true;

    return false;
}

template<class T, int d> void MATERIAL_PREDICATE_TESSELATED_VOLUME<T,d>::
MergeSubcellMaterial( const UNION_FIND<int>& merge_map ){

    for( int i = 1; i <= material_fragments.m; i++ ){
        ARRAY<int> master_fragment_list = material_fragments( merge_map.Find(i) );
        ARRAY<int> sub_fragment_list = material_fragments( i );
        for( int j = 1; j <= sub_fragment_list.m; j++)
            master_fragment_list.Append( sub_fragment_list(j) );
        master_fragment_list.Prune_Duplicates();
        material_fragments( merge_map.Find(i) ) = master_fragment_list;
    }
    for( int i = 1; i <= material_fragments.m; i++ ){
        material_fragments(i) = material_fragments( merge_map.Find(i) );
    }
}



template<class T, int d> void MATERIAL_PREDICATE_TESSELATED_VOLUME<T,d>::
ComputeNodalDistances( const ARRAY<TV>& nodal_positions, int subcell, ARRAY<T>& distances) const
{
    PHYSBAM_ASSERT( nodal_positions.m == T_MESH::vertices_per_cell );
    PHYSBAM_ASSERT( distances.m == T_MESH::vertices_per_cell );
    
    // Generate Boundary Primatives
    HASHTABLE<int, ARRAY< VECTOR<int,d> > > bp = GEOMETRY_PREDICATES<T>::GenerateBoundaryPrimatives( volume.mesh, material_fragments(subcell) );
    
    for( int n = 1; n <= T_MESH::vertices_per_cell; n++){
        //LOG::cout << "Computing distances for node "<< n << std::endl;
        //bool adjust_sign = false;
        int initial_sign = distances(n) < 0 ? -1: 1;
        bool inside = false;
        for( int f = 1; f <= material_fragments(subcell).m; f++){
            int e_id = material_fragments(subcell)(f);
            const T_ELEMENT element = volume.Get_Element( material_fragments(subcell)(f) );
            //LOG::cout << element << std::endl;
            const VECTOR<int,d+1> element_nodes = volume.mesh.elements( material_fragments(subcell)(f) );
            inside = inside || GEOMETRY_PREDICATES<T>::Inside( element, nodal_positions(n) );
            const ARRAY< VECTOR<int, d> > boundary_element_nodes = bp.Get( e_id );
            
            //LOG::cout << "\tChecking Fragment element "<< e_id << " with " << boundary_element_nodes.m << " boundary faces." << std::endl;

            for( int b=1; b <= boundary_element_nodes.m; b++){                
                T_BOUNDARY_ELEMENT be( volume.particles.X.Subset( boundary_element_nodes(b) ) );
                T dist = GEOMETRY_PREDICATES<T>::DistanceToElement( be, nodal_positions(n) );
                if( abs(distances(n)) > abs(dist) ){
                    //LOG::cout << "New potential distance selected: " << dist << std::endl;
                    distances(n) = abs(dist);
                }
            }           

        }
        if( inside || initial_sign == -1)
            distances(n) = -1 * abs(distances(n));
        else
            distances(n) = abs( distances(n) );
        
    }
}

template<class T, int d> bool MATERIAL_PREDICATE_TESSELATED_VOLUME<T,d>::
InsideMaterial(const int subcell_index,const int node_index) const
{
    return (corner_fragments(subcell_index)(node_index).m>0);
}

template<class T, int d> int MATERIAL_PREDICATE_TESSELATED_VOLUME<T,d>::
Node_Material_Representative(const int subcell_index,const int node_index) const
{
    if(corner_fragments(subcell_index)(node_index).m==0) return -1;
    return corner_fragments(subcell_index)(node_index)(1);
}

template<class T, int d> int MATERIAL_PREDICATE_TESSELATED_VOLUME<T,d>::
Subcell_Material_Representative(const int subcell_index) const
{
    if(material_fragments(subcell_index).m==0) return -1;
    return material_fragments(subcell_index)(1);
}

template class MATERIAL_PREDICATE_TESSELATED_VOLUME<float, 2>;
template class MATERIAL_PREDICATE_TESSELATED_VOLUME<float, 3>;
template class MATERIAL_PREDICATE_TESSELATED_VOLUME<double, 2>;
template class MATERIAL_PREDICATE_TESSELATED_VOLUME<double, 3>;

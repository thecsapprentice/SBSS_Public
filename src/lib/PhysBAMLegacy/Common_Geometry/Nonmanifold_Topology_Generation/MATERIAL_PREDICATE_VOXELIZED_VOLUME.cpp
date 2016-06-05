
#include "MATERIAL_PREDICATE_VOXELIZED_VOLUME.h"

using namespace PhysBAM;

template<class T, int d>
MATERIAL_PREDICATE_VOXELIZED_VOLUME<T,d>::
MATERIAL_PREDICATE_VOXELIZED_VOLUME(const int refinement_input,
                                    const ARRAY< VOXMAP > voxmap_input ) :
    refinement( refinement_input ), voxmap( voxmap_input ),
    voxmap_domain( T_INDEX::All_Ones_Vector(),
                   T_INDEX::All_Ones_Vector()*refinement_input ) {
}

template<class T, int d> void 
MATERIAL_PREDICATE_VOXELIZED_VOLUME<T,d>::
MaterialFragments(const T_MESH& mesh,
                  const HASHTABLE<PAIR<VECTOR<int,2>,int> >& linkage_list,
                  ARRAY<T_CELL>& sub_cells, 
                  HASHTABLE<int, int>& subcell_to_root_cell){
    

    int vertex_counter = 1;
    sub_cells.Clean_Memory();
    subcell_to_root_cell.Clean_Memory();
    material_fragments.Clean_Memory();
    corner_fragments.Clean_Memory();

    ARRAY< ARRAY< PAIR<T_INDEX,T_INDEX> > > root_cell_linkages;
    root_cell_linkages.Resize( mesh.elements.m );
    for( HASHTABLE_ITERATOR< PAIR<VECTOR<int,2>,int > > iterator(linkage_list); iterator.Valid(); iterator.Next() )
        root_cell_linkages( iterator.Key().y ).Append( PAIR<T_INDEX,T_INDEX>( FlatIndexToVoxel(iterator.Key().x(1)) , FlatIndexToVoxel(iterator.Key().x(2)) ) );
    

    PHYSBAM_ASSERT( mesh.elements.m == voxmap.m );
    int last_progress = 0;
#pragma omp parallel for shared(last_progress, mesh, linkage_list, vertex_counter, sub_cells, subcell_to_root_cell)
    for( int r = 1; r <= mesh.elements.m; r++ ){
        const VOXMAP& root_voxmap = voxmap( r );

        // Initialize a flood map for this cell
        ARRAY<int,T_INDEX> flood_map;
        flood_map.Resize( voxmap_domain );

        flood_map.Fill( -1 );
        for( RANGE_ITERATOR<d> iterator(voxmap_domain); iterator.Valid(); iterator.Next() ) {
            if( root_voxmap(iterator.Index() )){
                flood_map(iterator.Index()) = 0;
            } 
            else{
                flood_map(iterator.Index()) = -1; 
            }
        }

        // Label all disconnected regions
        ARRAY<T_INDEX> flood_fill_queue;
        int region_count = 0;
        for( RANGE_ITERATOR<d> iterator(voxmap_domain); iterator.Valid(); iterator.Next() ){
            const T_INDEX& index = iterator.Index();
            if(flood_map(index) == 0){
                flood_fill_queue.Insert(index,1);
                region_count++;
                //LOG::cout << "Found new zero entry at " << linear_index << ". Adding to queue." << std::endl;
                {
                    //LOG::cout << "Considering region " << region_count << std::endl;
                while( flood_fill_queue.Size() != 0){
                    
                    T_INDEX current = flood_fill_queue.Pop();
                    //LOG::cout << current << std::endl;

                    // Queue a linkage if one is present

                    for( int l = 1; l <= root_cell_linkages(r).m; l++){
                        if( current == root_cell_linkages(r)(l).x
                            && flood_map( root_cell_linkages(r)(l).y ) != region_count )
                            flood_fill_queue.Insert( root_cell_linkages(r)(l).y, 1 );
                        if( current == root_cell_linkages(r)(l).y 
                            && flood_map( root_cell_linkages(r)(l).x ) != region_count )
                            flood_fill_queue.Insert( root_cell_linkages(r)(l).x, 1 );
                    }


                    flood_map( current ) = region_count;
                    for( int dir = 1; dir <= d; dir++){
                        T_INDEX direction;
                        direction(dir)=1;
                        if( voxmap_domain.Lazy_Inside(current+direction) && 
                            flood_map( current+direction ) == 0){
                            flood_fill_queue.Insert( current+direction, 1 );}
                        direction(dir)=-1;
                        if( voxmap_domain.Lazy_Inside(current+direction) && 
                            flood_map( current+direction ) == 0){
                            flood_fill_queue.Insert( current+direction, 1 );}
                    }
                    flood_fill_queue.Prune_Duplicates();
                }
                }
            }              
        }

        //  Build Subcell Voxmap fragments
        
        //LOG::cout << "Discovered " << region_count << " regions." << std::endl;
        for( int i = 1; i <= region_count; i++ ) {
            
            VOXMAP subcell_fragment;
            subcell_fragment.Resize( voxmap_domain );
            subcell_fragment.Fill( false );

            for( RANGE_ITERATOR<d> iterator(voxmap_domain); iterator.Valid(); iterator.Next() )
                if( flood_map(iterator.Index()) == i )
                    subcell_fragment(iterator.Index()) = true;
         
#pragma omp critical 
            {
                T_CELL new_subcell;
                VECTOR<T_INDEX, T_MESH::vertices_per_cell> cfragments;

                new_subcell.index = mesh.elements(r).index;
                int v = 1;
                for( RANGE_ITERATOR<d> iterator(RANGE<T_INDEX>(T_INDEX(), T_INDEX()+1)); iterator.Valid(); iterator.Next(), v++ ){
                    new_subcell.vertices(v) = vertex_counter++;
                    T_INDEX index = iterator.Index();
                    index = index * (refinement-1);
                    index = index + 1;
                    if( subcell_fragment( index ) )
                        cfragments( v ) = index;
                    else
                        cfragments( v ) = T_INDEX::All_Ones_Vector()*(-1);
                }
                
                int s = sub_cells.Append( new_subcell );
                subcell_to_root_cell.Set(s, r);
                material_fragments.Append( subcell_fragment );
                corner_fragments.Append(cfragments);
            }
        }

#pragma omp critical 
        {
            if( (int)(((T)(last_progress+1) / mesh.elements.m) * 100) % 5 == 0  &&
                (int)(((T)(last_progress) / mesh.elements.m) * 100) % 5 != 0 ){
                LOG::cout << "Progress: "<< (int)(((T)(last_progress+1) / mesh.elements.m) * 100) << "%" << std::endl;
            }
            last_progress++;
        }
     
    }

}

template<class T, int d> bool 
MATERIAL_PREDICATE_VOXELIZED_VOLUME<T,d>::
IsMaterialContinuous(const int axis,
                     const int root_cell_high, const int sub_cell_high,
                     const int root_cell_low,  const int sub_cell_low) const {
    
    typedef VECTOR<int, d-1> T_FACE_INDEX;
    RANGE<T_FACE_INDEX> face_domain( T_FACE_INDEX::All_Ones_Vector(),
                                     T_FACE_INDEX::All_Ones_Vector()*refinement);

    T_FACE_INDEX permuteBase;
    for( int i = 1; i < d; i++)
        permuteBase(i) = i;

    bool isContinuous = false;
    for( RANGE_ITERATOR<d-1> iterator(face_domain); iterator.Valid(); iterator.Next() ) {
        T_INDEX base(iterator.Index());
        T_INDEX baseHigh = base; baseHigh(d) = 1;
        T_INDEX baseLow = base; baseLow(d) = refinement;
        T_INDEX permute = permuteBase.Insert( 3, axis );
        T_INDEX High = baseHigh.Permute( permute );
        T_INDEX Low = baseLow.Permute( permute );
        if( material_fragments(sub_cell_high)(High) && material_fragments(sub_cell_low)(Low) ){
            isContinuous = true;
            break;
        }
    }
    return isContinuous;
}


template<class T, int d> void 
MATERIAL_PREDICATE_VOXELIZED_VOLUME<T,d>::
MergeSubcellMaterial( const UNION_FIND<int>& merge_map ) {

    for( int i = 1; i <= material_fragments.m; i++ ){
        if(  merge_map.Find(i) == i )
            continue;
        VOXMAP master_fragment = material_fragments( merge_map.Find(i) );
        VOXMAP sub_fragment = material_fragments( i );
        for( RANGE_ITERATOR<d> iterator(voxmap_domain); iterator.Valid(); iterator.Next() )
            master_fragment(iterator.Index()) |= sub_fragment(iterator.Index());
        material_fragments( merge_map.Find(i) ) = master_fragment;
    }
    for( int i = 1; i <= material_fragments.m; i++ ){
        if(  merge_map.Find(i) == i )
            continue;
        material_fragments(i) = material_fragments( merge_map.Find(i) );
    }
}

template<class T, int d> void 
MATERIAL_PREDICATE_VOXELIZED_VOLUME<T,d>::
ComputeNodalDistances( const ARRAY<TV>& nodal_positions,
                       int subcell, ARRAY<T>& distances) const {

    PHYSBAM_ASSERT( nodal_positions.m == T_MESH::vertices_per_cell );
    PHYSBAM_ASSERT( distances.m == T_MESH::vertices_per_cell );

    PHYSBAM_NOT_IMPLEMENTED( "Computing Nodal Distances is not Supported for Voxelized Material Descriptions." );
}


template<class T, int d> bool
MATERIAL_PREDICATE_VOXELIZED_VOLUME<T,d>::
InsideMaterial(const int subcell_index,const int node_index) const {
    return (corner_fragments(subcell_index)(node_index) != (T_INDEX::All_Ones_Vector()*-1) );
}

template<class T, int d> int
MATERIAL_PREDICATE_VOXELIZED_VOLUME<T,d>::
Node_Material_Representative(const int subcell_index,const int node_index) const {
    if( corner_fragments(subcell_index)(node_index) == (T_INDEX::All_Ones_Vector() *(-1)) ) return -1;
    return VoxelIndexToFlat( corner_fragments(subcell_index)(node_index) );
}

template<class T, int d> int
MATERIAL_PREDICATE_VOXELIZED_VOLUME<T,d>::
Subcell_Material_Representative(const int subcell_index) const {
    T_INDEX first_active = T_INDEX::All_Ones_Vector()*(-1);
    for( RANGE_ITERATOR<d> iterator(voxmap_domain); iterator.Valid(); iterator.Next() )
        if( material_fragments(subcell_index)(iterator.Index()) ){
            first_active = iterator.Index();
            break;}
    if( first_active == (T_INDEX::All_Ones_Vector() *(-1)) ) return -1;
    return VoxelIndexToFlat( first_active );
}


template<class T, int d> const typename MATERIAL_PREDICATE_VOXELIZED_VOLUME<T,d>::VOXMAP& 
MATERIAL_PREDICATE_VOXELIZED_VOLUME<T,d>::SubcellVoxMaterialFragment(const int subcell_index ) const {
    return material_fragments(subcell_index);
}


template<class T, int d> int
MATERIAL_PREDICATE_VOXELIZED_VOLUME<T,d>::
VoxelIndexToFlat(const T_INDEX voxel_index) const {
    int flat_index = 0;
    for( int w = 1; w <= d; w++ ){
        int scale=1;
        for( int t = 0; t < d-w; t++ )
            scale *= refinement;
        flat_index += scale*(voxel_index(w)-1);
    }
    return flat_index;
}


template<class T, int d> typename MATERIAL_PREDICATE_VOXELIZED_VOLUME<T,d>::T_INDEX
MATERIAL_PREDICATE_VOXELIZED_VOLUME<T,d>::
FlatIndexToVoxel(const int flat_index) const {
    T_INDEX voxel_index;
    int q, r;
    int f = flat_index;
    for( int w = 1; w <= d; w++ ){
        int scale=1;
        for( int t = 0; t < d-w; t++ )
            scale *= refinement;
        
        q = f / scale;
        r = f % scale;
        voxel_index(w) = q;
        f = r;
    }   
    voxel_index+=1;
    return voxel_index;
}



template class MATERIAL_PREDICATE_VOXELIZED_VOLUME<float, 2>;
template class MATERIAL_PREDICATE_VOXELIZED_VOLUME<float, 3>;
template class MATERIAL_PREDICATE_VOXELIZED_VOLUME<double, 2>;
template class MATERIAL_PREDICATE_VOXELIZED_VOLUME<double, 3>;

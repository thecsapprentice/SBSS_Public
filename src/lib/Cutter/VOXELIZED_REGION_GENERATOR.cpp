#include "OVERRIDES.h"


#include "VOXELIZED_REGION_GENERATOR.h"
#include "PhysBAM_Tools/Data_Structures/HASHTABLE.h"
#include "RANGE_ITERATOR.h"
#include <cstring>
#include <cstdlib>

#ifdef ENABLE_PHYSBAM_IO
#include "Write_Output_3D.cpp"

#include <PhysBAM_Geometry/Geometry_Particles/REGISTER_GEOMETRY_READ_WRITE.h>
#include <PhysBAM_Geometry/Read_Write/Geometry/READ_WRITE_TRIANGULATED_SURFACE.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#endif

using namespace PhysBAM;


template <class T, int d> 
VOXELIZED_REGION_GENERATOR<T,d>::VOXELIZED_REGION_GENERATOR(int refinement_factor_input,
                                                            const T_GRID& fine_grid_input,
                                                            const T_FLAG_ARRAY& voxmap_input)
    : BASE( T_GRID( fine_grid_input.Counts() / refinement_factor_input + 1,
                    fine_grid_input.Domain() )),
      refinement_factor( refinement_factor_input),
      fine_grid( fine_grid_input ),
      voxmap (voxmap_input),
      vregions( BASE::regions )
{
    // Do some checks here to confirm grid sizes are consistent.
}



template <class T, int d> void
VOXELIZED_REGION_GENERATOR<T,d>::ResolveCuts()
{
    LOG::SCOPE scope("Voxelized_Region_Generator: Resolving Cuts");

    int vertex_counter = 1;

    // For each cell, 
    for( int cell = 1; cell <= root_cells.m; cell++)
        {
            ARRAY<int,T_INDEX> flood_map;
            const T_INDEX& coarse_index = root_grid_cell_mapping( cell );

            T_INDEX start = (coarse_index - 1) * refinement_factor+1;
            T_INDEX end   = (coarse_index - 1) * refinement_factor + refinement_factor;

            RANGE<T_INDEX> sub_cell_domain( start, end );
            RANGE<T_INDEX> flood_map_domain( T_INDEX::All_Ones_Vector(),
                                             T_INDEX::All_Ones_Vector()*refinement_factor);
            flood_map.Resize(flood_map_domain );
            for( RANGE_ITERATOR<d> iterator(sub_cell_domain), flood_iterator(flood_map_domain);
                 iterator.Valid() && flood_iterator.Valid(); iterator.Next(), flood_iterator.Next() )
                {
                    if( voxmap(iterator.Index() ))
                        flood_map(flood_iterator.Index()) = 0; 
                    else
                        flood_map(flood_iterator.Index()) = -1; 
                }
            
            ARRAY<T_INDEX> flood_fill_queue;
            int region_count = 0;
            
            for( RANGE_ITERATOR<d> iterator(flood_map_domain); iterator.Valid(); iterator.Next() ){
                const T_INDEX& linear_index = iterator.Index();
                if(flood_map(linear_index) == 0){
                    flood_fill_queue.Insert(linear_index,1);
                    region_count++;
                    //LOG::cout << "Found new zero entry at " << linear_index << ". Adding to queue." << std::endl;
                    while( flood_fill_queue.Size() != 0){
                        T_INDEX current = flood_fill_queue.Pop();
                        //LOG::cout << "Dequeing " << current << std::endl;
                        flood_map( current ) = region_count;
                        for( int dir = 1; dir <= d; dir++){
                            T_INDEX direction;
                            direction(dir)=1;
                            if( flood_map_domain.Lazy_Inside(current+direction) && 
                                flood_map( current+direction ) == 0){
                                flood_fill_queue.Insert( current+direction, 1 );}
                            direction(dir)=-1;
                            if( flood_map_domain.Lazy_Inside(current+direction) && 
                                flood_map( current+direction ) == 0){
                                flood_fill_queue.Insert( current+direction, 1 );}
                        }
                        flood_fill_queue.Prune_Duplicates();
                    }
                }              
            }

            //LOG::cout << "Regions found for cell " << coarse_index << " are " << region_count << std::endl;

            for( int i = 1; i <= region_count; i++ )
                {
                    T_CELL subcell;
                    for( int j=0; j<T_CELL::vertices_per_cell; j++){
                        subcell.indices[j] = vertex_counter++;
                        rootnode_to_subnode( root_cells(cell).indices[j] ).Append( subcell.indices[j] );
                    }

                    sub_cells.Append( subcell );
                    root_sub_mapping.Append( cell );
                    T_FLAG_ARRAY subcell_voxmap;
                    subcell_voxmap.Resize( flood_map_domain );
                    for( RANGE_ITERATOR<d> iterator(flood_map_domain); iterator.Valid(); iterator.Next() )
                            subcell_voxmap(iterator.Index()) = (flood_map(iterator.Index()) == i ? true : false);
                    subcell_voxmaps.Append( subcell_voxmap );
                    PHYSBAM_ASSERT( (sub_cells.m == root_sub_mapping.m) && (sub_cells.m == subcell_voxmaps.m));
                }

        }

    root_cell_to_subcell.Resize( root_cells.m );
    for( int i=1; i <= sub_cells.m; i++)
        root_cell_to_subcell( root_sub_mapping( i ) ).Append( i );
   
}

template <int d> VECTOR<int,d> FacePermute(VECTOR<int,d> base, AXIS axis, int offset)
{PHYSBAM_FATAL_ERROR();}

template <> VECTOR<int,1> FacePermute<1>(VECTOR<int,1> base, AXIS axis, int offset)
{
    typedef VECTOR<int,1> T_INDEX;
    
    T_INDEX permuted_base;
    switch( axis ){
    case A_X:
        permuted_base = base.Permute( T_INDEX(1) );
        permuted_base(1) = offset;
        break;
    case A_mX:
        permuted_base = base.Permute( T_INDEX(1) );
        break;
    default:
        PHYSBAM_FATAL_ERROR( "No Y-face, Z-face in one dimension." );               
    }
    
    return permuted_base;
}

template <> VECTOR<int,2> FacePermute<2>(VECTOR<int,2> base, AXIS axis, int offset)
{
    typedef VECTOR<int,2> T_INDEX;
    
    T_INDEX permuted_base;
    switch( axis ){
    case A_X:
        permuted_base = base.Permute( T_INDEX(2,1) );
        permuted_base(1) = offset;
        break;
    case A_mX:
        permuted_base = base.Permute( T_INDEX(2,1) );
        break;
    case A_Y:
        permuted_base = base.Permute( T_INDEX(1,2) );
        permuted_base(2) = offset;
        break;
    case A_mY:
        permuted_base = base.Permute( T_INDEX(1,2) );
        break;
    default:
        PHYSBAM_FATAL_ERROR( "No Z-face in two dimensions." );               
    }
    
    return permuted_base;
}

template <> VECTOR<int,3> FacePermute<3>(VECTOR<int,3> base, AXIS axis, int offset)
{
    typedef VECTOR<int,3> T_INDEX;
    
    T_INDEX permuted_base;
    switch( axis ){
    case A_X:
        permuted_base = base.Permute( T_INDEX(3,1,2) );
        permuted_base(1) = offset;
        break;
    case A_mX:
        permuted_base = base.Permute( T_INDEX(3,1,2) );
        break;
    case A_Y:
        permuted_base = base.Permute( T_INDEX(1,3,2) );
        permuted_base(2) = offset;
        break;
    case A_mY:
        permuted_base = base.Permute( T_INDEX(1,3,2) );
        break;
    case A_Z:
        permuted_base = base.Permute( T_INDEX(1,2,3) );
        permuted_base(3) = offset;
        break;
    case A_mZ:
        permuted_base = base.Permute( T_INDEX(1,2,3) );
        break;
    }
    
    return permuted_base;
}


template <class T, int d> bool
VOXELIZED_REGION_GENERATOR<T,d>::IsMaterialContinous(const int root_cellA, const int sub_cellA, const AXIS faceA,
                                                     const int root_cellB, const int sub_cellB, const AXIS faceB) const
{
    const T_INDEX& root_cell_A_index = root_grid_cell_mapping( root_cellA );
    const T_INDEX& root_cell_B_index = root_grid_cell_mapping( root_cellB );   
    const T_FLAG_ARRAY& subcell_A_voxmap = subcell_voxmaps(sub_cellA);
    const T_FLAG_ARRAY& subcell_B_voxmap = subcell_voxmaps(sub_cellB);
    
    const T_INDEX subcell_A_start = (root_cell_A_index - 1) * refinement_factor+1;
    const T_INDEX subcell_B_start = (root_cell_B_index - 1) * refinement_factor+1;

    typedef VECTOR<int, d-1> T_FACE_INDEX;
    RANGE<T_FACE_INDEX> refinement_domain( T_FACE_INDEX::All_Ones_Vector(),
                                           T_FACE_INDEX::All_Ones_Vector()*refinement_factor);
    bool isContinous = false;

    for( RANGE_ITERATOR<d-1> iterator(refinement_domain); iterator.Valid(); iterator.Next() )
        {
            T_INDEX base(iterator.Index());
            base(3) = 1; //Correct the last position.
            T_INDEX permuted_baseA = FacePermute(base, faceA, refinement_factor);
            T_INDEX permuted_baseB = FacePermute(base, faceB, refinement_factor);
            T_INDEX subcell_A_index = permuted_baseA;
            T_INDEX subcell_B_index = permuted_baseB;
         
            if( subcell_A_voxmap(subcell_A_index) && subcell_B_voxmap( subcell_B_index ) )
                isContinous = true;
        }

    return isContinous;
}



template <class T, int d> void
VOXELIZED_REGION_GENERATOR<T,d>::Generate()
{
    BASE::Generate();
    vregions = VOXELIZED_REGIONS<T,d>( (BASE::regions) );

    // We get away with this as the ordering of the regions shouldn't be different
    ARRAY<int> region_parents;

    for( int s = 1; s <= sub_cells.m; s++ )
        {
            int parent = region_sets.Find( s );
            int insert_region;

            //We don't need to insert new regions as there should not be any new ones.

            if( region_parents.Find( parent, insert_region) ){
                vregions.voxmap_regions(insert_region).Append(subcell_voxmaps(s));
            }
            else{
                region_parents.Append(parent);
                vregions.voxmap_regions.Append( typename VOXELIZED_REGIONS<T,d>::T_VOXMAP_REGION() );
                vregions.voxmap_regions.Last().Append(subcell_voxmaps(s));
            }
        }

    PHYSBAM_ASSERT( vregions.voxmap_regions.m == vregions.regions.m );
    for( int i = 1; i <= vregions.regions.m; i++)
        PHYSBAM_ASSERT( vregions.voxmap_regions(i).m == vregions.regions(i).m );

//    typedef float RW;
//    Initialize_Geometry_Particle();Initialize_Read_Write_Structures();
//    RW rw=RW();STREAM_TYPE stream_type(rw);
//    Write_Output( stream_type, *this, vregions, "output", 0 );
}


template <class T, int d> const VOXELIZED_REGIONS<T,d>* 
VOXELIZED_REGION_GENERATOR<T,d>::GetRegionData() const
{
    return &vregions;
};



template <class T, int d> const typename VOXELIZED_REGION_GENERATOR<T,d>::T_FLAG_ARRAY& 
VOXELIZED_REGION_GENERATOR<T,d>::GetSubcellVoxmap(const T_INDEX& cell_index, int duplicate) const
{
    int root_cell_index = BASE::root_index_to_root_id.Get( cell_index );
    return subcell_voxmaps( root_cell_to_subcell( root_cell_index )(duplicate) ); 
}








template class VOXELIZED_REGION_GENERATOR<float,2>;
template class VOXELIZED_REGION_GENERATOR<float,3>;
template class VOXELIZED_REGION_GENERATOR<double,2>;
template class VOXELIZED_REGION_GENERATOR<double,3>;

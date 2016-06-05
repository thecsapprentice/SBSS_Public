#include <PhysBAM_Geometry/Basic_Geometry/TRIANGLE_3D.h>
#include <PhysBAM_Geometry/Basic_Geometry_Intersections/BOX_TRIANGLE_3D_INTERSECTION.h>

#include "EMBEDDING_TOOLS.h"

using namespace PhysBAM;

//#####################################################################
// Function Get_Grid_Bounds
//#####################################################################
template< class T, int d>
void EMBEDDINGTOOLS<T,d>::Get_Grid_Bounds( const ARRAY<VECTOR<T,d> >& vertices, const ARRAY<VECTOR<int,d> >& triangles, const T dx, VECTOR<int, d>& cell_bounds, VECTOR<T, d>& min_corner)
{
    LOG::SCOPE scope("Embedding_Tools::Get_Grid_Bounds()");

    typedef VECTOR<T,d> TV;
    typedef VECTOR<int,d> T_INDEX;
    
    // Compute (padded) bounding box
    
    RANGE<TV> box=RANGE<TV>::Empty_Box();
    for(int v=1;v<=vertices.m;v++)
                box.Enlarge_To_Include_Point(vertices(v));
    box.Change_Size(.05*dx);
    
    // Subdivide box (potentially enlarging it) so that the cells along each dimension == 3 (mod 4)
    
    for(int v=1;v<=3;v++)
        cell_bounds(v)=((1+(int)(box.Edge_Lengths()(v)/dx))/4)*4+3;
    box.max_corner=box.min_corner+TV(cell_bounds)*dx;        
    min_corner = box.min_corner;              
    return;
}


//#####################################################################
// Function Rasterize
//#####################################################################

template<class T, int d>
void EMBEDDINGTOOLS<T,d>::Rasterize( const ARRAY<VECTOR<T,d> >& vertices, const ARRAY<VECTOR<int,d> >& triangles, const GRID<VECTOR<T,d> >& domain, const RANGE<VECTOR<int,d> >& unpadded_domain, const RANGE<VECTOR<int,d> >& padded_domain,  ARRAY< bool, VECTOR<int, d> >& voxmap,  ARRAY< bool, VECTOR<int, d> >& voxmap_node)
{
    LOG::SCOPE scope("Embedding_Tools::Rasterize()");

    typedef VECTOR<T,d> TV;
    typedef VECTOR<int,d> T_INDEX;
    
    const RANGE<T_INDEX>& unpadded_cell_domain = unpadded_domain;
    const RANGE<T_INDEX>& padded_cell_domain = padded_domain;

    ARRAY< char, T_INDEX > flood_voxmap(padded_cell_domain);
    flood_voxmap.Fill(-1);
    voxmap.Fill( false );
    const GRID<TV>& grid=domain;

    {
        LOG::SCOPE scope("Marking all grid cells intersecting mesh.");
        
        for(int t=1;t<=triangles.m;t++){
            TRIANGLE_3D<T> triangle(vertices(triangles(t)(1)),vertices(triangles(t)(2)),vertices(triangles(t)(3)));
            RANGE<TV> triangle_box=triangle.Bounding_Box();
            RANGE<T_INDEX> cell_range;
            cell_range.min_corner=grid.Cell(triangle_box.min_corner,0);
            cell_range.max_corner=grid.Cell(triangle_box.max_corner,0);

            if(!cell_range.Lazy_Intersection(unpadded_cell_domain))
                continue;
            cell_range = RANGE<T_INDEX>::Intersect(cell_range,unpadded_cell_domain);

            for(RANGE_ITERATOR<d> iterator(cell_range);iterator.Valid();iterator.Next()){
                const T_INDEX& cell_index=iterator.Index();
                RANGE<TV> cell=grid.Cell_Domain(cell_index);
                if(INTERSECTION::Intersects<T>(cell,triangle,0))
                    flood_voxmap(cell_index)=1;}}
    }
    // Flood fill the exterior region up to the boundary
    
    {
        LOG::SCOPE scope("Flood filling all nodes to determine outside region.");
    
        QUEUE<T_INDEX> queue(unpadded_cell_domain.Surface_Area());
        for(RANGE_ITERATOR<d> iterator(padded_cell_domain);iterator.Valid();iterator.Next()){
            const T_INDEX& cell_index=iterator.Index();
            if(unpadded_cell_domain.Lazy_Outside(cell_index)){
            flood_voxmap(cell_index)=0;
            queue.Safe_Enqueue(cell_index);}}
        while(!queue.Empty()){
            const T_INDEX cell_index=queue.Dequeue();
            for(int q = 1; q <=d; q++){
                for(int r = -1; r <=1; r++){
                    T_INDEX neighbor_index=cell_index;
                    neighbor_index(q) += r;
                    if(padded_cell_domain.Lazy_Inside(neighbor_index) && flood_voxmap(neighbor_index)==-1){
	                    flood_voxmap(neighbor_index)=0;
                        queue.Safe_Enqueue(neighbor_index);}}}}

        for(RANGE_ITERATOR<d> iterator(padded_cell_domain);iterator.Valid();iterator.Next()){
            const T_INDEX& cell_index=iterator.Index();
            if(flood_voxmap(cell_index)==-1 || flood_voxmap(cell_index)==1){
                flood_voxmap(cell_index)=1;
                voxmap(cell_index) = true;
            }
        }
    }

    voxmap_node.Fill( 0 );

    {
        LOG::SCOPE scope("Marking all remaining cells as inside.");
        {int interior_count=0;
            for(RANGE_ITERATOR<d> iterator(padded_cell_domain);iterator.Valid();iterator.Next())
                if(flood_voxmap(iterator.Index())==1){
                    interior_count++;
                    T_INDEX nodes[8];
                    domain.Nodes_In_Cell_From_Minimum_Corner_Node( iterator.Index(), nodes );
                    for( int i=0; i<8; i++) voxmap_node( nodes[i] )=1;
                }
            LOG::cout<<"Interior fine cells : "<<interior_count<<std::endl;}
        
        {int interior_count=0;
            for(RANGE_ITERATOR<d> iterator(padded_cell_domain);iterator.Valid();iterator.Next())
                if(flood_voxmap(iterator.Index())==0)
                interior_count++;
            LOG::cout<<"Exterior fine cells : "<<interior_count<<std::endl;}
        
    }
   
    
}


//#####################################################################
// Function Coarsen
//#####################################################################

template<class T, int d>
void EMBEDDINGTOOLS<T,d>::Coarsen( const GRID<VECTOR<T,d> >& fine_domain, const RANGE<VECTOR<int,d> >& fine_unpadded_domain,
                              const RANGE<VECTOR<int,d> >& fine_padded_domain, const ARRAY< bool, VECTOR<int, d> >& fine_voxmap,
                              const GRID<VECTOR<T,d> >& coarse_domain, const RANGE<VECTOR<int,d> >& coarse_unpadded_domain,
                              const RANGE<VECTOR<int,d> >& coarse_padded_domain, ARRAY< CELL_TYPE , VECTOR<int, d> >& coarse_voxmap)
{
    LOG::SCOPE scope("Embedding_Tools::Coarsen()");

    typedef VECTOR<T,d> TV;
    typedef VECTOR<int,d> T_INDEX;
    
    int interior_count=0;
    int exterior_count=0;
    int boundary_count=0;

    // We start at the same place.
    PHYSBAM_ASSERT( fine_domain.Xmin() == coarse_domain.Xmin() );
    int coarsen_ratio = coarse_domain.Maximum_Edge_Length() / fine_domain.Maximum_Edge_Length();
    
    coarse_voxmap.Fill( EXTERIOR_CELL_TYPE );
    for( RANGE_ITERATOR<d> coarse_iterator( coarse_unpadded_domain ); coarse_iterator.Valid(); coarse_iterator.Next() ){
        bool is_full = true;
        bool is_empty = true;
        
        const T_INDEX& index = coarse_iterator.Index();
        T_INDEX fine_start = (index-1)*coarsen_ratio+1;
        T_INDEX fine_end = (index-1)*coarsen_ratio+T_INDEX::All_Ones_Vector()*coarsen_ratio;

        for( RANGE_ITERATOR<d> fine_iterator( RANGE<T_INDEX>(fine_start, fine_end)); fine_iterator.Valid(); fine_iterator.Next()){
                if( fine_voxmap( fine_iterator.Index() ) == true)
                    is_empty = false;
                if( fine_voxmap( fine_iterator.Index() ) == false)
                    is_full = false;
            }
            
        if( is_empty ){
            coarse_voxmap( index ) = EXTERIOR_CELL_TYPE;
            exterior_count++;
        }
        if( is_full ){
            coarse_voxmap( index ) = INTERIOR_CELL_TYPE;
            interior_count++;
        }
        if( !is_full && !is_empty ){
            coarse_voxmap( index ) = INTERIOR_CELL_TYPE; 
            interior_count++;
        }
    }

    LOG::cout << "Coarse Grid has " << interior_count << " INTERIOR_CELLS" << std::endl;
    LOG::cout << "Coarse Grid has " << exterior_count << " EXTERIOR_CELLS" << std::endl;
    LOG::cout << "Coarse Grid has " << boundary_count << " BOUNDARY_CELLS" << std::endl;
}


//#####################################################################
// Function Coarsen
//#####################################################################

template<class T, int d>
void EMBEDDINGTOOLS<T,d>::Coarsen_Density( const GRID<VECTOR<T,d> >& fine_domain, const RANGE<VECTOR<int,d> >& fine_unpadded_domain,
                              const RANGE<VECTOR<int,d> >& fine_padded_domain, const ARRAY< bool, VECTOR<int, d> >& fine_voxmap,
                              const GRID<VECTOR<T,d> >& coarse_domain, const RANGE<VECTOR<int,d> >& coarse_unpadded_domain,
                              const RANGE<VECTOR<int,d> >& coarse_padded_domain, ARRAY< float , VECTOR<int, d> >& coarse_voxmap)
{
    LOG::SCOPE scope("Embedding_Tools::Coarsen_Density()");

    typedef VECTOR<T,d> TV;
    typedef VECTOR<int,d> T_INDEX;
    
    int interior_count=0;
    int exterior_count=0;
    int boundary_count=0;

    // We start at the same place.
    PHYSBAM_ASSERT( fine_domain.Xmin() == coarse_domain.Xmin() );
    int coarsen_ratio = coarse_domain.Maximum_Edge_Length() / fine_domain.Maximum_Edge_Length();
    
    coarse_voxmap.Fill( EXTERIOR_CELL_TYPE );
    for( RANGE_ITERATOR<d> coarse_iterator( coarse_unpadded_domain ); coarse_iterator.Valid(); coarse_iterator.Next() ){
        int is_full = 0;
        int is_empty = 0;
        
        const T_INDEX& index = coarse_iterator.Index();
        T_INDEX fine_start = (index-1)*coarsen_ratio+1;
        T_INDEX fine_end = (index-1)*coarsen_ratio+T_INDEX::All_Ones_Vector()*coarsen_ratio;

        for( RANGE_ITERATOR<d> fine_iterator( RANGE<T_INDEX>(fine_start, fine_end)); fine_iterator.Valid(); fine_iterator.Next()){
                if( fine_voxmap( fine_iterator.Index() ) == true)
                    is_full++;
                else
                    is_empty++;                    
            }
            
        coarse_voxmap( index ) = (float)(is_full)/((float)(is_full + is_empty));
    }
}


//#####################################################################
// Function Multilinear_Interpolation_Stencil
//#####################################################################
template<class T,int d>
STENCIL<T,d> EMBEDDINGTOOLS<T,d>::
Multilinear_Interpolation_Stencil(const VECTOR<int, d>& cell_index,
                                  const VECTOR<T,d>& multilinear_coordinates)
{
    typedef STENCIL<T,d> T_STENCIL;
    typedef VECTOR<T,d> TV;
    typedef VECTOR<int,d> T_INDEX;

    T_STENCIL interpolation_stencil;
    for(RANGE_ITERATOR<d> iterator(RANGE<T_INDEX>(cell_index,cell_index+1));iterator.Valid();iterator.Next()){
        T_INDEX node_index=iterator.Index();
        T weight=(T)1.;
        for(int v=1;v<=d;v++) weight*=node_index(v)==cell_index(v)?(T)1.-multilinear_coordinates(v):multilinear_coordinates(v);
        interpolation_stencil.Insert(node_index,weight);}
    return interpolation_stencil;
}


//#####################################################################
// Function Deformation
//#####################################################################
template<class T,int d>
VECTOR<T,d> EMBEDDINGTOOLS<T,d>::
Deformation(const GRID<VECTOR<T,d> >& domain,
            const VECTOR< ARRAY<T,VECTOR<int,d> >, d>& displacement,
            const VECTOR<int,d>& cell_index,
            const VECTOR<T,d>& multilinear_coordinates,
            bool displace_only)
{
    typedef STENCIL<T,d> T_STENCIL;
    typedef VECTOR<T,d> TV;
    typedef VECTOR<int,d> T_INDEX;

    T_STENCIL interpolation_stencil=Multilinear_Interpolation_Stencil(cell_index,multilinear_coordinates);
    TV result;
    for(int v=1;v<=d;v++) result(v)=interpolation_stencil*displacement(v);
    if( displace_only )
        return result;
    else
        return result+domain.Node(cell_index)+domain.dX*multilinear_coordinates;
}

//#####################################################################
// Function Multilinear_Interpolation
//#####################################################################
template<class T,int d>
VECTOR<T,d> EMBEDDINGTOOLS<T,d>::
Multilinear_Interpolation(const GRID<VECTOR<T,d> >& domain,
                          const VECTOR<int,d>& cell_index,
                          const VECTOR<T,d>& multilinear_coordinates)
{
    typedef STENCIL<T,d> T_STENCIL;
    typedef VECTOR<T,d> TV;
    typedef VECTOR<int,d> T_INDEX;

    return domain.Node(cell_index)+domain.dX*multilinear_coordinates;
}

template struct EMBEDDINGTOOLS<float, 3>;
template struct EMBEDDINGTOOLS<double, 3>;

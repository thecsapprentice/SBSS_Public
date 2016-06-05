#define DEBUG_LEVELSET

#include <PhysBAM_Tools/Log/DEBUG_UTILITIES.h>
#include <vector>
#include <PhysBAM_Geometry/Basic_Geometry/BOX.h>
#include <PhysBAM_Geometry/Basic_Geometry/TORUS.h>
#include <PhysBAM_Geometry/Grids_Uniform_Computations/LEVELSET_MAKER_UNIFORM.h>
#include <PhysBAM_Geometry/Tessellation/TORUS_TESSELLATION.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/FREE_PARTICLES.h>

#ifdef DEBUG_LEVELSET
#include <PhysBAM_Geometry/Solids_Geometry/DEFORMABLE_GEOMETRY_COLLECTION.h>
#include <PhysBAM_Geometry/Geometry_Particles/REGISTER_GEOMETRY_READ_WRITE.h>
#include <PhysBAM_Tools/Read_Write/Grids_Uniform/READ_WRITE_GRID.h>
#include <PhysBAM_Tools/Read_Write/Grids_Uniform_Arrays/READ_WRITE_ARRAYS.h>
#endif

#include "PHYSBAM_LEVELSET_COLLISION.h"

using namespace PhysBAM;

template<class T, int d>
PHYSBAM_LEVELSET_COLLISION<T,d>::PHYSBAM_LEVELSET_COLLISION(const std::vector<std::vector<float> > vertices, const std::vector<std::vector<int> > triangles, const float refinement) : grid(NULL), phi(NULL), levelset(NULL)
{
    Initialize_Geometry_Particle();Initialize_Read_Write_Structures();

    STREAM_TYPE stream_type((float)0);
    typedef VECTOR<T,d> TV;
    typedef VECTOR<int,d> T_INDEX;
    TRIANGULATED_SURFACE<T>* surface=TRIANGULATED_SURFACE<T>::Create();
    for( int v=1;v<=vertices.size();v++){
        int p = surface->particles.array_collection->Add_Element();
        for(int w=1;w<=d;w++)
            surface->particles.X(p)(w)=vertices[v-1][w-1];
    }
    for( int t=1;t<=triangles.size();t++){
        VECTOR<int,d> tri;
        for(int v=1; v<=d; v++)
            tri(v) = triangles[t-1][v-1]+1;
        surface->mesh.elements.Append(tri);
    }
    LOG::cout << "Vertices: " << surface->particles.X.m << std::endl;
    surface->Update_Number_Nodes();

    #ifdef DEBUG_LEVELSET
    DEFORMABLE_GEOMETRY_COLLECTION<TV> collection(surface->particles);
    collection.Add_Structure(surface);
    FILE_UTILITIES::Create_Directory("collision/0");
    collection.Write(stream_type,"collision",0,0,true);
    #endif

    LOG::cout<<"Building levelset"<<std::endl;

    T dx=.02f,margin=2.5f*dx;
    surface->Update_Bounding_Box();
    RANGE<TV> range=surface->bounding_box->Bounding_Box().Thickened(margin);
    T_INDEX counts=T_INDEX(range.Edge_Lengths()/dx)+1;
    range.max_corner=range.min_corner+TV(counts-1)*dx;
    grid = new GRID<TV>(counts,range);
    LOG::cout<<"Grid = "<<*grid<<std::endl;
    phi = new ARRAY<T,T_INDEX>(grid->Domain_Indices());
     

    LEVELSET_MAKER_UNIFORM<T> levelset_maker;
    levelset_maker.Verbose_Mode(true);
    bool success=levelset_maker.Compute_Level_Set(*surface,*grid,*phi);
    LOG::cout<<"Compute_Level_Set returned value : "<<success<<std::endl;

    #ifdef DEBUG_LEVELSET
    FILE_UTILITIES::Write_To_File(stream_type,"collision/0/levelset",*grid,*phi);
    #endif

    levelset = new LEVELSET_IMPLICIT_OBJECT<TV>( *grid, *phi );

};

template<class T, int d>  
PHYSBAM_LEVELSET_COLLISION<T,d>::~PHYSBAM_LEVELSET_COLLISION()
{
    if(levelset)
        delete levelset;
    if(grid)
        delete grid;
    if(phi)
        delete phi;
}

template<class T, int d> T 
PHYSBAM_LEVELSET_COLLISION<T,d>::Phi_Implementation( const T world_location[d]) const
{
    PHYSBAM_ASSERT(levelset!=NULL);
    VECTOR<T,d> wl;
    for(int i=1;i<=d;i++)
        wl(i) = world_location[i-1];
    return levelset->Extended_Phi( wl );
};

template<class T, int d> void
PHYSBAM_LEVELSET_COLLISION<T,d>::Normal_Implementation( const T world_location[d], T normal[d]) const
{
    PHYSBAM_ASSERT(levelset!=NULL);
    VECTOR<T,d> n;
    VECTOR<T,d> wl;
    for(int i=1;i<=d;i++)
        wl(i) = world_location[i-1];
    n = levelset->Extended_Normal(wl);
    for(int i=1;i<=d;i++)
        normal[i-1] = n(i);
};

template class PHYSBAM_LEVELSET_COLLISION<float, 3>;
template class PHYSBAM_LEVELSET_COLLISION<double, 3>;
//template class PHYSBAM_LEVELSET_COLLISION<float, 2>;
//template class PHYSBAM_LEVELSET_COLLISION<double, 2>;

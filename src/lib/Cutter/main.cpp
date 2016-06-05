//#####################################################################
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Geometry/Geometry_Particles/REGISTER_GEOMETRY_READ_WRITE.h>
#include <PhysBAM_Geometry/Read_Write/Geometry/READ_WRITE_TRIANGULATED_SURFACE.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>

#include <fstream>
#include <string>

#include <PhysBAM_Tools/Grids_Uniform/GRID.h>

#include "VOXELIZED_REGION_GENERATOR.h"
#include "RANGE_ITERATOR.h"

#include "Write_Output_3D.cpp"

using namespace PhysBAM;

void Initialize( int size, int refinement_factor,
                 GRID< VECTOR<float,3> >& coarse_grid,  GRID< VECTOR<float,3> >& fine_grid,
                 ARRAY< bool, VECTOR<int,3> >& voxmap )
{
    // Set up an example by initializing grid & voxmap
    typedef float T;
    const int d=3;
    typedef VECTOR<float,d> TV;
    typedef VECTOR<int,d> T_INDEX;
    typedef GRID<TV> T_GRID;

    int coarse_size = size;
    int fine_size = size*refinement_factor;

    coarse_grid.Initialize(T_INDEX::All_Ones_Vector()*coarse_size+1,RANGE<TV>::Unit_Box());
    fine_grid.Initialize(T_INDEX::All_Ones_Vector()*fine_size+1,RANGE<TV>::Unit_Box());
    voxmap.Resize(fine_grid.Cell_Indices());
    
    RANGE<T_INDEX> fine_grid_domain( fine_grid.Cell_Indices() );

    for(RANGE_ITERATOR<3> iterator(fine_grid.Cell_Indices());iterator.Valid();iterator.Next()){
#if 1
        const T_INDEX& index=iterator.Index();
        TV X_canonical=fine_grid.Center(index);
        TV dX_canonical=X_canonical-.5;
        int i=argmax(abs(dX_canonical(1)),abs(dX_canonical(2)),abs(dX_canonical(3)));
        T sj=dX_canonical(i%3+1)/(1e-4+abs(dX_canonical(i)));
        T sk=dX_canonical((i+1)%3+1)/(1e-4+abs(dX_canonical(i)));
        if(abs(dX_canonical(i))<1e-4) sj=sk=0;
        sj=.5*(sj+sin(half_pi*sj));
        sk=.5*(sk+sin(half_pi*sk));
        T offset=(1.-cos(two_pi*sj))*(1.-cos(two_pi*sk));
        T radius=(dX_canonical).Magnitude();
        if(radius<.17+0.02*offset && fine_grid_domain.Lazy_Inside(index))
            voxmap(iterator.Index())=true;
        else
            voxmap(iterator.Index())=false;
#else
        if((fine_grid.Center(iterator.Index())-.5).Magnitude()<.3)
            voxmap(iterator.Index())=true;
#endif
        
#if 0
        if( iterator.Index()(1) < iterator.Index()(2)+1 && iterator.Index()(1) > iterator.Index()(2)-1 &&
            iterator.Index()(2) < iterator.Index()(1)+1 && iterator.Index()(2) > iterator.Index()(1)-1 )
            voxmap( iterator.Index())=false;

        if( iterator.Index()(1) < iterator.Index()(3)+1 && iterator.Index()(1) > iterator.Index()(3)-1 &&
            iterator.Index()(3) < iterator.Index()(1)+1 && iterator.Index()(3) > iterator.Index()(1)-1 )
            voxmap( iterator.Index())=false;

        if((fine_grid.Center(iterator.Index())-.5).Magnitude()<.2 && 
           (fine_grid.Center(iterator.Index())-.5).Magnitude()>.19)
            voxmap( iterator.Index())=false;
#else
        if((fine_grid.Center(iterator.Index()))(1) > 0.49 &&
           (fine_grid.Center(iterator.Index()))(1) <= 0.51 &&
           (fine_grid.Center(iterator.Index()))(2) > .5)
            voxmap( iterator.Index())=false;



#endif

    }
}


int main(int argc,char* argv[])
{
    typedef float T;
    typedef float RW;
    RW rw=RW();STREAM_TYPE stream_type(rw);
    static const int d=3;

    LOG::Initialize_Logging();

    Initialize_Geometry_Particle();Initialize_Read_Write_Structures();

    PARSE_ARGS parse_args;
    //parse_args.Add_String_Argument("-input","","file","Input .tri filename");
    //parse_args.Add_String_Argument("-output","","file","Input .tri_raw filename");
    parse_args.Parse(argc,argv);

    typedef VECTOR<float,3> TV;
    typedef VECTOR<int,3> T_INDEX;
    typedef GRID<TV> T_GRID;
    
    T_GRID fine_grid;
    T_GRID coarse_grid;
    ARRAY< bool, VECTOR<int,3> > voxmap;
    
    const int size = 5;
    const int refinement_factor = 20;

    Initialize( size, refinement_factor, 
                coarse_grid, fine_grid, 
                voxmap );

    VOXELIZED_REGION_GENERATOR<T,d> embedded_voxelized_object( refinement_factor,
                                                               fine_grid,
                                                               voxmap);
    
    VOXELIZED_REGIONS<T,d>* regions = embedded_voxelized_object.Generate();   
    Write_Output( stream_type, embedded_voxelized_object, *regions, "output", 0 );
    delete regions;

    LOG::Finish_Logging();

    return 0;
}
//#####################################################################

//#####################################################################
// Copyright 2013, Nathan Mitchell, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class EMBEDDED_VOXELIZED_OBJECT
//#####################################################################
#ifndef __EMBEDDED_VOXELIZED_OBJECT__
#define __EMBEDDED_VOXELIZED_OBJECT__

#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include "RANGE_ITERATOR.h"

namespace PhysBAM{

template<class T,int d>
class EMBEDDED_VOXELIZED_OBJECT
{
    STATIC_ASSERT(d==2||d==3);

    enum{ vertices_per_cell=(d==2)?4:8 };

    typedef VECTOR<T,d> TV;
    typedef GRID<TV> T_GRID;
    typedef VECTOR<int,d> T_INDEX;
    typedef bool T_FLAG;
    typedef ARRAY<T_FLAG,T_INDEX> T_FLAG_ARRAY;
    typedef VECTOR<int,vertices_per_cell> T_ELEMENT;

    // ???????
    typedef VECTOR<int,d-1> T_FACE_INDEX;
    typedef ARRAY<T_FLAG,T_FACE_INDEX> T_FLAG_FACE_ARRAY;

    typedef ARRAY<bool> T_FACE_FLAGS;
    typedef ARRAY<unsigned char> T_FACE_FLAGS_COMPACT;

public:

    T_GRID grid;         // This is a grid of *nodes*, i.e. it's sized (nx+1,ny+1,nz+1) if (nx,ny,nz) are cell counts
    T_FLAG_ARRAY voxmap; // This is the high-resolution voxmap that defines our object

    T_GRID coarse_grid;                    // This is a grid "step" times coarser than grid, with the same min_corner
    ARRAY<T_ELEMENT> coarse_mesh;          // This is a quad/hex mesh, each element of which contains the vertex IDs for its 4/8 vertices
                                           // Vertices are listed in RANGE_ITERATOR ordering
    ARRAY<T_INDEX> coarse_nodes;           // This maps a vertex ID to an (i,j,k) location in the coarse_grid
    ARRAY<T_INDEX> coarse_cell_of_element; // This has the size=coarse_mesh.m, maps elements to CELL (i,j,k) indices

    void Initialize(/* ???? */)
    {
        const int size=100;
        grid=T_GRID(T_INDEX::All_Ones_Vector()*size+1,RANGE<TV>::Unit_Box());
        voxmap.Resize(grid.Cell_Indices());

        for(RANGE_ITERATOR<d> iterator(grid.Cell_Indices());iterator.Valid();iterator.Next())
            if((grid.Center(iterator.Index())-.5).Magnitude()<.3)
                voxmap(iterator.Index())=true;

        // Set up an example by initializing grid & voxmap
    }

    void Embed(const int step)
    {

        // make coarse_grid

        for(RANGE_ITERATOR<d> iterator(coarse_grid.Cell_Indices());iterator.Valid();iterator.Next())
        
    }
    
    

//#####################################################################
//#####################################################################
};
}
#endif

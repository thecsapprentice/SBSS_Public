//#####################################################################
// Copyright 2014, Nathan Mitchell.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class READ_WRITE_NONMANIFOLD_LEVELSET_MESH
//#####################################################################
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
#ifndef __READ_WRITE_NONMANIFOLD_LEVELSET_MESH__
#define __READ_WRITE_NONMANIFOLD_LEVELSET_MESH__

#include <PhysBAM_Tools/Read_Write/Utilities/READ_WRITE.h>
#include <PhysBAM_Tools/Read_Write/Arrays/READ_WRITE_ARRAY.h>
#include <PhysBAM_Tools/Read_Write/Data_Structures/READ_WRITE_PAIR.h>
#include <PhysBAM_Tools/Read_Write/Data_Structures/READ_WRITE_TRIPLE.h>
#include <PhysBAM_Tools/Read_Write/Data_Structures/READ_WRITE_HASHTABLE.h>
#include <Nonmanifold_Implicit_Objects/NONMANIFOLD_LEVELSET_MESH.h>

#include <iostream>
namespace PhysBAM{



    template<class RW,class T,int d>
    class Read_Write<NONMANIFOLD_LEVELSET_MESH_TRANSITION_FACE<T,d>,RW>
{
public:
    static void Read(std::istream& input,NONMANIFOLD_LEVELSET_MESH_TRANSITION_FACE<T,d>& object)
    {
        Read_Binary<RW>(input,object.type);
        Read_Binary<RW>(input,object.left_cells);
        Read_Binary<RW>(input,object.right_cells);
    }

    static void Write(std::ostream& output,const NONMANIFOLD_LEVELSET_MESH_TRANSITION_FACE<T,d>& object)
    {
        Write_Binary<RW>(output,object.type);
        Write_Binary<RW>(output,object.left_cells);
        Write_Binary<RW>(output,object.right_cells);
    }
};


    template<class RW,class T,int d>
    class Read_Write<NONMANIFOLD_LEVELSET_MESH_CELL<T,d>,RW>
{
public:
    static void Read(std::istream& input,NONMANIFOLD_LEVELSET_MESH_CELL<T,d>& object)
    {
        Read_Binary<RW>(input,object.cell_index);
        Read_Binary<RW>(input,object.nodes);
        Read_Binary<RW>(input,object.transition);
        Read_Binary<RW>(input,object.cell_neighbors);
    }

    static void Write(std::ostream& output,const NONMANIFOLD_LEVELSET_MESH_CELL<T,d>& object)
    {
        Write_Binary<RW>(output,object.cell_index);
        Write_Binary<RW>(output,object.nodes);
        Write_Binary<RW>(output,object.transition);
        Write_Binary<RW>(output,object.cell_neighbors);
    }
};


    template<class RW,class T,int d>
    class Read_Write<NONMANIFOLD_LEVELSET_MESH<T,d>,RW>
{
public:
    static void Read(std::istream& input,NONMANIFOLD_LEVELSET_MESH<T,d>& object)
    {
        Read_Binary<RW>(input,object.dx);
        Read_Binary<RW>(input,object.node_mesh_count);
        Read_Binary<RW>(input,object.node_domain);
        Read_Binary<RW>(input,object.initialized);
        Read_Binary<RW>(input,object.transition_faces);
        Read_Binary<RW>(input,object.cells);
        Read_Binary<RW>(input,object.cell_is_grid);
        Read_Binary<RW>(input,object.node_locations);
        ARRAY<RANGE<VECTOR<T,d> > > leaves;
        Read_Binary<RW>(input,leaves);
        object.aabb_hierarchy.Clean_Memory();
        if( leaves.m > 0 )
            object.aabb_hierarchy.Set_Leaf_Boxes(leaves, true);
    }

    static void Write(std::ostream& output,const NONMANIFOLD_LEVELSET_MESH<T,d>& object)
    {
        Write_Binary<RW>(output,object.dx);
        Write_Binary<RW>(output,object.node_mesh_count);
        Write_Binary<RW>(output,object.node_domain);
        Write_Binary<RW>(output,object.initialized);
        Write_Binary<RW>(output,object.transition_faces);
        Write_Binary<RW>(output,object.cells);
        Write_Binary<RW>(output,object.cell_is_grid);
        Write_Binary<RW>(output,object.node_locations);
        ARRAY<RANGE<VECTOR<T,d> > > leaves;
        leaves.Resize( object.aabb_hierarchy.leaves );
        for( int l=1; l<= object.aabb_hierarchy.leaves; l++ )
            leaves(l) =  object.aabb_hierarchy.box_hierarchy(l);
        LOG::cout << "leaves: "<< leaves.m << std::endl;
        Write_Binary<RW>(output,leaves);
    }
};
}
#endif
#endif

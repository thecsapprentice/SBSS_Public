//#####################################################################
// Copyright 2013, Nathan Mitchell, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENVDB_LEVELSET_COLLISION
//#####################################################################
#ifndef __OPENVDB_LEVELSET_COLLISION_H__
#define __OPENVDB_LEVELSET_COLLISION_H__

#include "COLLISION_INTERFACE.h"
//#include <openvdb/openvdb.h>
#include <vector>

namespace PhysBAM{

    template<class T, int d>
    class OPENVDB_LEVELSET_COLLISION : public COLLISION_INTERFACE<T,d>{
    public:
        OPENVDB_LEVELSET_COLLISION(const std::vector<std::vector<float> > vertices,
                                   const std::vector<std::vector<int> > triangles,
                                   const float refinement);      
    private:

/*
        typedef typename openvdb::tree::Tree4<T, 5,4,3 >::Type PhiTree;
        typedef typename openvdb::tree::Tree4<openvdb::math::Vec3<T>, 5,4,3 >::Type NormalTree;

        typedef openvdb::Grid<PhiTree> PhiGrid;
        typedef openvdb::Grid<NormalTree> NormalGrid;
        typedef T PhiType;
        typedef openvdb::math::Vec3<T> NormalType;
        typedef openvdb::math::Vec3<T> WorldCoord;
*/
        virtual T Phi_Implementation( const T world_location[d]) const;
        virtual void Normal_Implementation( const T world_location[d], T normal[d]) const;
/*
        typename PhiGrid::Ptr phi_grid;
        typename NormalGrid::Ptr normal_grid;
        openvdb::math::Transform::Ptr xtrans;
*/
    };
}

#endif

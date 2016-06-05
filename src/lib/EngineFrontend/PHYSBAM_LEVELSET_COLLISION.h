//#####################################################################
// Copyright 2013, Nathan Mitchell, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class PHYSBAM_LEVELSET_COLLISION
//#####################################################################
#ifndef __PHYSBAM_LEVELSET_COLLISION_H__
#define __PHYSBAM_LEVELSET_COLLISION_H__

#include "COLLISION_INTERFACE.h"
#include <vector>
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Geometry/Implicit_Objects_Uniform/LEVELSET_IMPLICIT_OBJECT.h>

namespace PhysBAM{

    template<class T, int d>
    class PHYSBAM_LEVELSET_COLLISION : public COLLISION_INTERFACE<T,d>{
    public:
        PHYSBAM_LEVELSET_COLLISION(const std::vector<std::vector<float> > vertices,
                                   const std::vector<std::vector<int> > triangles,
                                   const float refinement);      
        ~PHYSBAM_LEVELSET_COLLISION();
        
    private:

        typedef VECTOR<T,d> TV;
        typedef VECTOR<int,d> T_INDEX;
        
        LEVELSET_IMPLICIT_OBJECT<TV>* levelset;
        GRID<TV>* grid;
        ARRAY<T,T_INDEX>* phi;

        virtual T Phi_Implementation( const T world_location[d]) const;
        virtual void Normal_Implementation( const T world_location[d], T normal[d]) const;

    };
}

#endif

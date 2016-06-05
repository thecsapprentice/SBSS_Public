//#####################################################################
// Copyright 2013, Nathan Mitchell, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class PHYSBAM_RIGID_BODY_LEVELSET_COLLISION
//#####################################################################
#ifndef __PHYSBAM_RIGID_BODY_LEVELSET_COLLISION_H__
#define __PHYSBAM_RIGID_BODY_LEVELSET_COLLISION_H__

#include "COLLISION_INTERFACE.h"
#include <vector>
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Geometry/Implicit_Objects_Uniform/LEVELSET_IMPLICIT_OBJECT.h>
#include <PhysBAM_Geometry/Implicit_Objects/IMPLICIT_OBJECT_TRANSFORMED.h>

#include <PhysBAM_Tools/Matrices/FRAME.h>


namespace PhysBAM{

    template<class T, int d>
    class PHYSBAM_RIGID_BODY_LEVELSET_COLLISION : public COLLISION_INTERFACE<T,d>{
    public:
        typedef VECTOR<T,d> TV;
        typedef VECTOR<int,d> T_INDEX;

        PHYSBAM_RIGID_BODY_LEVELSET_COLLISION(const std::vector<std::vector<float> > vertices,
                                   const std::vector<std::vector<int> > triangles,
                                   const float refinement);      
        ~PHYSBAM_RIGID_BODY_LEVELSET_COLLISION();
        
        FRAME<TV> transform;
    private:
        
        LEVELSET_IMPLICIT_OBJECT<TV>* levelset;
        IMPLICIT_OBJECT_TRANSFORMED<TV, FRAME<TV> >* levelset_transformed;

        GRID<TV>* grid;
        ARRAY<T,T_INDEX>* phi;


        virtual T Phi_Implementation( const T world_location[d]) const;
        virtual void Normal_Implementation( const T world_location[d], T normal[d]) const;

    };
}

#endif

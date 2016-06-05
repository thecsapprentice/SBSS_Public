//#####################################################################
// Copyright 2010-2012, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class PARAMETER_ANIMATION
//#####################################################################
#ifndef __PARAMETER_ANIMATION_h__
#define __PARAMETER_ANIMATION_h__

#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>

namespace PhysBAM{

class PARAMETER_ANIMATION{

    typedef float T;
    typedef VECTOR<T,3> TV;

public:

    virtual ~PARAMETER_ANIMATION() {}

    virtual void Set_Bone_Frames(ARRAY<FRAME<TV> >& bone_frames,const T time) const
    {PHYSBAM_NOT_IMPLEMENTED();}

    virtual void Set_Bone_Twists(ARRAY<TWIST<TV> >& bone_twists,const T time) const
    {PHYSBAM_NOT_IMPLEMENTED();}

    virtual void Set_Muscle_Activations(ARRAY<T>& muscle_activations,const T time) const
    {PHYSBAM_NOT_IMPLEMENTED();}
    
//#####################################################################    
//#####################################################################    
};
}
#endif

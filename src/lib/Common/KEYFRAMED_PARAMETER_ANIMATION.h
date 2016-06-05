//#####################################################################
// Copyright 2010-2012, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class KEYFRAMED_PARAMETER_ANIMATION
//#####################################################################
#ifndef __KEYFRAMED_PARAMETER_ANIMATION_h__
#define __KEYFRAMED_PARAMETER_ANIMATION_h__

#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Matrices/FRAME.h>
#include <PhysBAM_Tools/Vectors/TWIST.h>
#include "PARAMETER_ANIMATION.h"

namespace PhysBAM{

class KEYFRAMED_PARAMETER_ANIMATION:public PARAMETER_ANIMATION{

    typedef float T;
    typedef VECTOR<T,3> TV;

    GRID<VECTOR<T,1> > time_grid;
    ARRAY<ARRAY<FRAME<TV> > > bone_keyframes;
    ARRAY<ARRAY<T> > muscle_keyframes;
    
//#####################################################################    
public:
    void Initialize_Bone_Keyframes(const std::string& bone_keyframes_filename,const int number_of_bones);
    void Initialize_Muscle_Keyframes(const std::string& muscle_keyframes_filename,const int number_of_muscles);
    void Set_Bone_Frames(ARRAY<FRAME<TV> >& bone_frames,const T time) const;
    void Set_Bone_Twists(ARRAY<TWIST<TV> >& bone_twists,const T time) const;
    void Set_Muscle_Activations(ARRAY<T>& muscle_activations,const T time) const;
//#####################################################################    
};
}
#endif

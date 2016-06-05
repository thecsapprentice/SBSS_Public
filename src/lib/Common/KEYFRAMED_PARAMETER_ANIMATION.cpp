//#####################################################################
// Copyright 2010-2012, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class KEYFRAMED_PARAMETER_ANIMATION
//#####################################################################
#include <PhysBAM_Tools/Matrices/FRAME.h>
#include <PhysBAM_Tools/Vectors/TWIST.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <fstream>
#include "KEYFRAMED_PARAMETER_ANIMATION.h"
//#define ENABLE_LOG_MESSAGES
using namespace PhysBAM;
//#####################################################################
// Function Initialize_Bone_Keyframes
//#####################################################################
void KEYFRAMED_PARAMETER_ANIMATION::
Initialize_Bone_Keyframes(const std::string& bone_keyframes_filename,const int number_of_bones)
{
#ifdef ENABLE_LOG_MESSAGES
    LOG::SCOPE scope("KEYFRAMED_PARAMETER_ANIMATION::Initialize_Bone_Keyframes()");
#endif

    std::ifstream input(bone_keyframes_filename.c_str());
    if(!input.is_open())
        PHYSBAM_FATAL_ERROR("Could not read bone keyframes from file : "+bone_keyframes_filename);
#ifdef ENABLE_LOG_MESSAGES
    LOG::cout<<"Reading bone keyframes from file : "<<bone_keyframes_filename<<std::endl;
#endif
    T tmin,tmax;int number_of_frames;input >> tmin >> tmax >> number_of_frames;
    time_grid=GRID<VECTOR<T,1> >(number_of_frames,tmin,tmax);
#ifdef ENABLE_LOG_MESSAGES
    LOG::cout<<"Time grid : "<<time_grid<<std::endl;
#endif
    bone_keyframes.Resize(number_of_frames);
    for(int frame=1;frame<=number_of_frames;frame++){
        bone_keyframes(frame).Resize(number_of_bones);
        for(int bone=1;bone<=number_of_bones;bone++){
            TV t;input >> t.x >> t.y >> t.z;
            QUATERNION<T> q;input >> q.v.x >> q.v.y >> q.v.z >> q.s;
            FRAME<TV> f(t,ROTATION<TV>::From_Quaternion(q));
            bone_keyframes(frame)(bone)=f;}}
    input.close();    
}
//#####################################################################
// Function Initialize_Muscle_Keyframes
//#####################################################################
void KEYFRAMED_PARAMETER_ANIMATION::
Initialize_Muscle_Keyframes(const std::string& muscle_keyframes_filename,const int number_of_muscles)
{
#ifdef ENABLE_LOG_MESSAGES
	LOG::SCOPE scope("KEYFRAME_PARAMETER_ANIMATION::Initialize_Muscle_Keyframes()");
#endif
	std::ifstream input(muscle_keyframes_filename.c_str());
	if(!input.is_open())
		PHYSBAM_FATAL_ERROR("Could not read muscle keyframes from file : "+muscle_keyframes_filename);
#ifdef ENABLE_LOG_MESSAGES
	LOG::cout<<"Reading muscle keyframes from file : "<<muscle_keyframes_filename<<std::endl;
#endif
	T tmin, tmax; int number_of_frames;input >> tmin >> tmax >> number_of_frames;
	time_grid = GRID<VECTOR<T,1> >(number_of_frames, tmin, tmax);   //shall we compare this with the one computed from bone keyframe file?
#ifdef ENABLE_LOG_MESSAGES
	LOG::cout<<"Time grid : "<<time_grid<<std::endl;
#endif
	muscle_keyframes.Resize(number_of_frames);
	for(int frame=1;frame<=number_of_frames;frame++){
		muscle_keyframes(frame).Resize(number_of_muscles);
		for(int muscle=1;muscle<=number_of_muscles;muscle++){
			T a; input >> a;
			muscle_keyframes(frame)(muscle)=a;}}
	input.close();
	
}
//#####################################################################
// Function Set_Bone_Frames
//#####################################################################
void KEYFRAMED_PARAMETER_ANIMATION::
Set_Bone_Frames(ARRAY<FRAME<TV> >& bone_frames,const T time) const
{
    int previous_frame=time_grid.Clamped_Index(VECTOR<T,1>(time)).x;
    int next_frame=min(previous_frame+1,time_grid.counts(1));
    T ratio=(time-time_grid.Node(previous_frame).x)/time_grid.dX(1);
    if(ratio<0.) ratio=0.;
    if(ratio>1.) ratio=1.;
    PHYSBAM_ASSERT(bone_frames.m==bone_keyframes(previous_frame).m);
    for(int bone=1;bone<=bone_frames.m;bone++)
        bone_frames(bone)=FRAME<TV>::Interpolation(bone_keyframes(previous_frame)(bone),bone_keyframes(next_frame)(bone),ratio);
}
//#####################################################################
// Function Set_Bone_Twists
//#####################################################################
void KEYFRAMED_PARAMETER_ANIMATION::
Set_Bone_Twists(ARRAY<TWIST<TV> >& bone_twists,const T time) const
{
    int previous_frame=time_grid.Clamped_Index(VECTOR<T,1>(time)).x;
    int next_frame=min(previous_frame+1,time_grid.counts(1));
    T ratio=(time-time_grid.Node(previous_frame).x)/time_grid.dX(1);
    if(ratio<0.) ratio=0.;
    if(ratio>1.) ratio=1.;
    PHYSBAM_ASSERT(bone_twists.m==bone_keyframes(previous_frame).m);
	T one_over_dt = 1.0/time_grid.dX(1);
	for(int bone=1;bone<=bone_twists.m;bone++)
	{
		bone_twists(bone).linear = one_over_dt*(bone_keyframes(next_frame)(bone).t-bone_keyframes(previous_frame)(bone).t);
		bone_twists(bone).angular = one_over_dt*(bone_keyframes(next_frame)(bone).r*bone_keyframes(previous_frame)(bone).r.Inverse()).Rotation_Vector();
	}
	
}
//#####################################################################
// Function Set_Muscle_Activations
//#####################################################################
void KEYFRAMED_PARAMETER_ANIMATION::
Set_Muscle_Activations(ARRAY<T>& muscle_activations,const T time) const
{
    int previous_frame=time_grid.Clamped_Index(VECTOR<T,1>(time)).x;
    int next_frame=min(previous_frame+1,time_grid.counts(1));
    T ratio=(time-time_grid.Node(previous_frame).x)/time_grid.dX(1);
    if(ratio<0.) ratio=0.;
    if(ratio>1.) ratio=1.;
    PHYSBAM_ASSERT(muscle_activations.m==muscle_keyframes(previous_frame).m);
	for(int muscle=1;muscle<=muscle_activations.m;muscle++)
		muscle_activations(muscle) = muscle_keyframes(previous_frame)(muscle)*(1.0-ratio)+muscle_keyframes(next_frame)(muscle)*ratio;

    
}
//#####################################################################

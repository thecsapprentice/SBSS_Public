//#####################################################################
// Copyright 2014, Nathan Mitchell, Raj Setaluri.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class NONMANIFOLD_LEVELSET_POLICY_UNIFORM 
//#####################################################################
#ifndef __NONMANIFOLD_LEVELSET_POLICY_UNIFORM__
#define __NONMANIFOLD_LEVELSET_POLICY_UNIFORM__

namespace PhysBAM{

template<class T,int d> class NONMANIFOLD_LEVELSET_POLICY;

template<class T,int d> class NONMANIFOLD_LEVELSET_UNIFORM;
template<class T> class NONMANIFOLD_LEVELSET_2D;
template<class T> class NONMANIFOLD_LEVELSET_3D;

template<class T> struct NONMANIFOLD_LEVELSET_POLICY<T,2>
{
    typedef NONMANIFOLD_LEVELSET_2D<T> NONMANIFOLD_LEVELSET;
};

template<class T> struct NONMANIFOLD_LEVELSET_POLICY<T,3>
{
    typedef NONMANIFOLD_LEVELSET_3D<T> NONMANIFOLD_LEVELSET;
};
}
#endif

//#####################################################################
// Copyright 2014, Raj Setaluri.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class NONMANIFOLD_LEVELSET_ISOCONTOUR_3D
//#####################################################################
#ifndef __NONMANIFOLD_LEVELSET_ISOCONTOUR_3D__
#define __NONMANIFOLD_LEVELSET_ISOCONTOUR_3D__

#include <Nonmanifold_Implicit_Objects/NONMANIFOLD_LEVELSET_3D.h>
namespace PhysBAM{

template<class T> class TRIANGULATED_SURFACE;

template<class T>
class NONMANIFOLD_LEVELSET_ISOCONTOUR_3D
{
    typedef VECTOR<T,3> TV;typedef VECTOR<int,3> TV_INT;
public:
//#####################################################################
    static TRIANGULATED_SURFACE<T>* Create_Triangulated_Surface_From_Levelset(const NONMANIFOLD_LEVELSET_3D<T>& levelset,const T isovalue);
//#####################################################################
};
}
#endif

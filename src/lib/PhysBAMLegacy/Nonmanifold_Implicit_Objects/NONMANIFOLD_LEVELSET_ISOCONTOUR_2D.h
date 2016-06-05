//#####################################################################
// Copyright 2014, Raj Setaluri.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class NONMANIFOLD_LEVELSET_ISOCONTOUR_2D
//#####################################################################
#ifndef __NONMANIFOLD_LEVELSET_ISOCONTOUR_2D__
#define __NONMANIFOLD_LEVELSET_ISOCONTOUR_2D__

#include <Nonmanifold_Implicit_Objects/NONMANIFOLD_LEVELSET_2D.h>
namespace PhysBAM{

template<class T> class SEGMENTED_CURVE_2D;

template<class T>
class NONMANIFOLD_LEVELSET_ISOCONTOUR_2D
{
    typedef VECTOR<T,2> TV;typedef VECTOR<int,2> TV_INT;
    typedef NONMANIFOLD_LEVELSET_MESH<T,2> T_NONMANIFOLD_MESH;
    enum{faces_per_cell=4};
public:
//#####################################################################
    static SEGMENTED_CURVE_2D<T>* Create_Segmented_Curve_From_Levelset(const NONMANIFOLD_LEVELSET_2D<T>& levelset,const T isovalue);
//#####################################################################
};
}
#endif

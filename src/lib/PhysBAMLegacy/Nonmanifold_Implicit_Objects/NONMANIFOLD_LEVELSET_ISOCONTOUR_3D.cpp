//#####################################################################
// Copyright 2014, Raj Setaluri.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class NONMANIFOLD_LEVELSET_ISOCONTOUR_3D
//#####################################################################
#include <Nonmanifold_Implicit_Objects/NONMANIFOLD_LEVELSET_ISOCONTOUR_3D.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
using namespace PhysBAM;
//#####################################################################
// Function Create_Triangulated_Surface_From_Levelset
//#####################################################################
template<class T> TRIANGULATED_SURFACE<T>* NONMANIFOLD_LEVELSET_ISOCONTOUR_3D<T>::
Create_Triangulated_Surface_From_Levelset(const NONMANIFOLD_LEVELSET_3D<T>& levelset,const T isovalue)
{
    PHYSBAM_NOT_IMPLEMENTED();
}
//#####################################################################
template class NONMANIFOLD_LEVELSET_ISOCONTOUR_3D<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class NONMANIFOLD_LEVELSET_ISOCONTOUR_3D<double>;
#endif

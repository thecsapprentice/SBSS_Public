//#####################################################################
// Copyright 2014, Nathan Mitchell, Raj Setaluri.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class NONMANIFOLD_LEVELSET_3D
//#####################################################################
#include <Nonmanifold_Implicit_Objects/NONMANIFOLD_LEVELSET_3D.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T> NONMANIFOLD_LEVELSET_3D<T>::
NONMANIFOLD_LEVELSET_3D(T_MESH& mesh_input,T_ARRAY_SCALAR& phi_input)
    :BASE(mesh_input,phi_input)
{}
//#####################################################################
// Destructor
//#####################################################################
template<class T> NONMANIFOLD_LEVELSET_3D<T>::
~NONMANIFOLD_LEVELSET_3D()
{
}
//#####################################################################
// Function Compute_Normals
//#####################################################################
template<class T> void NONMANIFOLD_LEVELSET_3D<T>::
Compute_Normals()
{
    PHYSBAM_NOT_IMPLEMENTED();
}
//#####################################################################
// Function Lazy_Normal
//#####################################################################
template<class T> typename NONMANIFOLD_LEVELSET_3D<T>::TV NONMANIFOLD_LEVELSET_3D<T>::
Lazy_Normal(const int cell_index,const TV& weights) const
{
    TV normal;
    T phi_node[T_MESH::nodes_per_cell];
    for(int i=1;i<=T_MESH::nodes_per_cell;i++) // 000 001 010 011 100 101 110 111
        phi_node[i-1]=phi(mesh.cells(cell_index).nodes(i));
    normal.x=((T)1.-weights.y)*((T)1.-weights.z)*(phi_node[4]-phi_node[0]) + 
             ((T)1.-weights.y)*(      weights.z)*(phi_node[5]-phi_node[1]) +
             (      weights.y)*((T)1.-weights.z)*(phi_node[6]-phi_node[2]) +
             (      weights.y)*(      weights.z)*(phi_node[7]-phi_node[3]);
    normal.y=((T)1.-weights.x)*((T)1.-weights.z)*(phi_node[2]-phi_node[0]) + 
             ((T)1.-weights.x)*(      weights.z)*(phi_node[3]-phi_node[1]) +
             (      weights.x)*((T)1.-weights.z)*(phi_node[6]-phi_node[4]) +
             (      weights.x)*(      weights.z)*(phi_node[7]-phi_node[5]);
    normal.z=((T)1.-weights.x)*((T)1.-weights.y)*(phi_node[1]-phi_node[0]) + 
             ((T)1.-weights.x)*(      weights.y)*(phi_node[3]-phi_node[2]) +
             (      weights.x)*((T)1.-weights.y)*(phi_node[5]-phi_node[4]) +
             (      weights.x)*(      weights.y)*(phi_node[7]-phi_node[6]);
    normal*=((T)1./mesh.dx); // NOTE: not strictly needed since we are normalizing
    return -normal.Normalized();
}
//#####################################################################
template class NONMANIFOLD_LEVELSET_3D<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class NONMANIFOLD_LEVELSET_3D<double>;
#endif

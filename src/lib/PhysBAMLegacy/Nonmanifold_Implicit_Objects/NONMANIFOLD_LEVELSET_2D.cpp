//#####################################################################
// Copyright 2014, Nathan Mitchell, Raj Setaluri.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class NONMANIFOLD_LEVELSET_2D
//#####################################################################
#include <Nonmanifold_Implicit_Objects/NONMANIFOLD_LEVELSET_2D.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T> NONMANIFOLD_LEVELSET_2D<T>::
NONMANIFOLD_LEVELSET_2D(T_MESH& mesh_input,T_ARRAY_SCALAR& phi_input)
    :BASE(mesh_input,phi_input)
{}
//#####################################################################
// Destructor
//#####################################################################
template<class T> NONMANIFOLD_LEVELSET_2D<T>::
~NONMANIFOLD_LEVELSET_2D()
{
}
//#####################################################################
// Function Compute_Normals
//#####################################################################
template<class T> void NONMANIFOLD_LEVELSET_2D<T>::
Compute_Normals()
{
    PHYSBAM_NOT_IMPLEMENTED();
}
//#####################################################################
// Function Lazy_Normal
//#####################################################################
template<class T> typename NONMANIFOLD_LEVELSET_2D<T>::TV NONMANIFOLD_LEVELSET_2D<T>::
Lazy_Normal(const int cell_index,const TV& weights) const
{
    TV normal;
    T phi_node[T_MESH::nodes_per_cell];
    for(int i=1;i<=T_MESH::nodes_per_cell;i++) // 00 01 10 11
        phi_node[i-1]=phi(mesh.cells(cell_index).nodes(i));
    normal.x=((T)1.-weights.y)*(phi_node[2]-phi_node[0])+(weights.y)*(phi_node[3]-phi_node[1]);
    normal.y=((T)1.-weights.x)*(phi_node[1]-phi_node[0])+(weights.x)*(phi_node[3]-phi_node[2]);
    normal*=((T)1./mesh.dx); // NOTE: not strictly needed since we are normalizing
    return -normal.Normalized();
}
//#####################################################################
template class NONMANIFOLD_LEVELSET_2D<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class NONMANIFOLD_LEVELSET_2D<double>;
#endif

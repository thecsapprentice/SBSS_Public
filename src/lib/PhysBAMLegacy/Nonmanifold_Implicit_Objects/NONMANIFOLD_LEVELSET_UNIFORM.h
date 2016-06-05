//#####################################################################
// Copyright 2014, Nathan Mitchell, Raj Setaluri.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class NONMANIFOLD_LEVELSET_UNIFORM
//#####################################################################
#ifndef __NONMANIFOLD_LEVELSET_UNIFORM__
#define __NONMANIFOLD_LEVELSET_UNIFORM__

#include <Nonmanifold_Implicit_Objects/NONMANIFOLD_LEVELSET_MESH.h>

namespace PhysBAM{

    template< class T, int d >
    class NONMANIFOLD_LEVELSET_UNIFORM {

    public:
        typedef VECTOR<T,d> TV;
        typedef VECTOR<int,d> TV_INT;
        typedef NONMANIFOLD_LEVELSET_MESH<T,d> T_MESH;
        typedef HYBRID_ARRAY<T,d> T_ARRAY_SCALAR;
        typedef HYBRID_ARRAY<TV,d> T_ARRAY_VECTOR;
        typedef typename T_ARRAY_SCALAR::INDEX T_INDEX;

        T_MESH& mesh;
        T_ARRAY_SCALAR& phi;
        T_ARRAY_VECTOR* normals;

        const T dx,one_over_dx;

        NONMANIFOLD_LEVELSET_UNIFORM( T_MESH& mesh_input, T_ARRAY_SCALAR& phi_input );
        virtual ~NONMANIFOLD_LEVELSET_UNIFORM();

        T Lazy_Phi(const int cell_index,const TV& weights) const;
        T Distance_To_Closest_Point(const int input_cell_index,const TV& input_weights,int& output_cell_index,TV& output_weights,TV& normal) const;

        virtual TV Lazy_Normal(const int cell_index,const TV& weights) const=0;
    };
};

#endif

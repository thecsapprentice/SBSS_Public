//#####################################################################
// Copyright 2014, Nathan Mitchell, Raj Setaluri.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class NONMANIFOLD_LEVELSET_2D  
//##################################################################### 
#ifndef __NONMANIFOLD_LEVELSET_2D__
#define __NONMANIFOLD_LEVELSET_2D__

#include <Nonmanifold_Implicit_Objects/NONMANIFOLD_LEVELSET_UNIFORM.h>

namespace PhysBAM {

    template< class T >
    class NONMANIFOLD_LEVELSET_2D : public NONMANIFOLD_LEVELSET_UNIFORM<T,2> {
        
        
    public:
        typedef NONMANIFOLD_LEVELSET_UNIFORM<T,2> BASE;
        typedef typename BASE::TV TV;
        typedef typename BASE::T_MESH T_MESH;
        typedef typename BASE::T_ARRAY_SCALAR T_ARRAY_SCALAR;
        typedef typename BASE::T_ARRAY_VECTOR T_ARRAY_VECTOR;

        using BASE::mesh;
        using BASE::phi;
        using BASE::normals;

        NONMANIFOLD_LEVELSET_2D( T_MESH& mesh_input, T_ARRAY_SCALAR& phi_input );
        virtual ~NONMANIFOLD_LEVELSET_2D();

        void Compute_Normals();
        virtual TV Lazy_Normal(const int cell_index,const TV& weights) const;
    };
};

#endif




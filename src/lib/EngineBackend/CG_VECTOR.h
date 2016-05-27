//#####################################################################
// Copyright 2011, Taylor Patterson, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class CG_VECTOR
//#####################################################################
#ifndef __CG_VECTOR__
#define __CG_VECTOR__

#include "OVERRIDES.h"

#include <PhysBAM_Tools/Krylov_Solvers/KRYLOV_VECTOR_BASE.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Vectors/VECTOR_2D.h>
#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>

#include "CG_POLICY.h"

namespace PhysBAM{

    //template<class T,int d> class NONLINEAR_ELASTICITY;

    template<class POLICY>
    class CG_VECTOR:public KRYLOV_VECTOR_BASE<typename POLICY::ELASTICITY::SCALAR>
{
    typedef typename POLICY::ELASTICITY::SCALAR T;
    typedef typename POLICY::ELASTICITY T_ELASTICITY;
    typedef KRYLOV_VECTOR_BASE<T> BASE;
    typedef typename POLICY::STATE T_STATE;
    const T_ELASTICITY& ne;
    static const int d = T_ELASTICITY::DIM;

    T_STATE& state;

    int X_StrideNode;
    int Y_StrideNode;

    int X_StrideCell;
    int Y_StrideCell;

public:
 CG_VECTOR(const T_ELASTICITY& elasticity_input,T_STATE& state_input)
     :ne(elasticity_input),state(state_input) 
    {
        VECTOR<int,3> node_grid_dims = ne.padded_node_domain.max_corner - ne.padded_node_domain.min_corner +1 ;
        X_StrideNode = node_grid_dims(2)*node_grid_dims(3);
        Y_StrideNode =                   node_grid_dims(3);

        VECTOR<int,3> cell_grid_dims = ne.padded_cell_domain.max_corner - ne.padded_cell_domain.min_corner +1;
        X_StrideCell = cell_grid_dims(2)*cell_grid_dims(3);
        Y_StrideCell =                   cell_grid_dims(3);
    }

    static T_STATE& State(BASE& base)
    {return (dynamic_cast<CG_VECTOR&>(base)).state;}

    static const T_STATE& State(const BASE& base)
    {return (dynamic_cast<const CG_VECTOR&>(base)).state;}

//#####################################################################
    BASE& operator+=(const BASE& bv);
    BASE& operator-=(const BASE& bv);
    BASE& operator*=(const T a);
    BASE& operator*=(const BASE& bv);
    void Copy(const T c,const BASE& bv);
    void Copy(const T c1,const BASE& bv1,const BASE& bv2);
    int Raw_Size() const;
    T& Raw_Get(int i);
    void Clear();
    void Set(const T a);
//#####################################################################   
};
}
#endif

//#####################################################################
// Copyright 2011, Taylor Patterson, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class CG_SYSTEM
//#####################################################################
#ifndef __CG_SYSTEM__
#define __CG_SYSTEM__

#include "OVERRIDES.h"

#include <PhysBAM_Tools/Krylov_Solvers/KRYLOV_SYSTEM_BASE.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Vectors/VECTOR_2D.h>
#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>
#include <Common/ALIGNED_ARRAY.h>
#include <EngineInterface/CELL_TYPE.h>

namespace PhysBAM{

//template<class T,int d> class NONLINEAR_ELASTICITY;

template<class T_ELASTICITY>
    class CG_SYSTEM:public KRYLOV_SYSTEM_BASE<typename T_ELASTICITY::SCALAR>
{
    typedef typename T_ELASTICITY::SCALAR T;
    typedef KRYLOV_SYSTEM_BASE<T> BASE;
    typedef KRYLOV_VECTOR_BASE<T> VECTOR_BASE;
    typedef typename T_ELASTICITY::T_STATE T_STATE;
    static const int d = T_ELASTICITY::DIM;

    typedef VECTOR<int,d> T_INDEX;
    typedef ARRAY_VIEW<T,T_INDEX> T_SCALAR_VARIABLE_VIEW;
    typedef VECTOR<T_SCALAR_VARIABLE_VIEW,d> T_VECTOR_VARIABLE_VIEW;
    typedef ARRAY_VIEW<const T,T_INDEX> T_SCALAR_VARIABLE_VIEW_CONST;
    typedef VECTOR<T_SCALAR_VARIABLE_VIEW_CONST,d> T_VECTOR_VARIABLE_VIEW_CONST;
    typedef ARRAY<int,T_INDEX> T_FLAG;
    typedef ARRAY<CELL_TYPE,T_INDEX> T_CELL_TYPE_FIELD;    

    const T_ELASTICITY& elasticity;

    int X_StrideNode;
    int Y_StrideNode;

    int X_StrideCell;
    int Y_StrideCell;

    // For Preconditioner
    const VECTOR_BASE* D;
    const VECTOR_BASE* D_primal;
    const VECTOR_BASE* D_dual;
    const VECTOR_BASE* omega;

//#####################################################################
public:
    CG_SYSTEM(const T_ELASTICITY& elasticity_input);
    void Multiply(const VECTOR_BASE& v,VECTOR_BASE& result) const;
    double Inner_Product(const VECTOR_BASE& x,const VECTOR_BASE& y) const;
    T Convergence_Norm(const VECTOR_BASE& x) const;
    void Project(VECTOR_BASE& x) const;
    void Set_Boundary_Conditions(VECTOR_BASE& x) const;
    void Project_Nullspace(VECTOR_BASE& x) const;
    void Apply_Preconditioner(const VECTOR_BASE& r, VECTOR_BASE& z) const;

    void SetDiagonalMatrix(const VECTOR_BASE* D_input) { D = D_input; };
    void SetDiagonalPrimalMatrix(const VECTOR_BASE* D_input) { D_primal = D_input; };
    void SetDiagonalDualMatrix(const VECTOR_BASE* D_input) { D_dual = D_input; };
    void SetOmega(const VECTOR_BASE* omega_input) { omega = omega_input; };
 //#####################################################################
};
}
#endif

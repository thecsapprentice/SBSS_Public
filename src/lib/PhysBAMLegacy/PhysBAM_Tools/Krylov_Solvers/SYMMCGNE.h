//#####################################################################
// Copyright 2006-2008, Geoffrey Irving, Craig Schroeder, Andrew Selle, Tamar Shinar, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SYMMCGNE
//#####################################################################
#ifndef __SYMMCGNE__
#define __SYMMCGNE__

#include <PhysBAM_Tools/Krylov_Solvers/KRYLOV_SOLVER.h>
namespace PhysBAM{

// see Golub and Van Loan, 10.2.6, p. 546 for details

template<class T>
class SYMMCGNE:public KRYLOV_SOLVER<T>
{
    typedef KRYLOV_SOLVER<T> BASE;
public:
    using BASE::restart_iterations;using BASE::residual_magnitude_squared;using BASE::iterations_used;using BASE::print_diagnostics;using BASE::print_residuals;
    using BASE::nullspace_measure;using BASE::nullspace_tolerance;using BASE::Solve;

    SYMMCGNE()
    {}

    virtual ~SYMMCGNE();

//#####################################################################
    bool Solve(const KRYLOV_SYSTEM_BASE<T>& system,KRYLOV_VECTOR_BASE<T>& x,const KRYLOV_VECTOR_BASE<T>& b,KRYLOV_VECTOR_BASE<T>& q,KRYLOV_VECTOR_BASE<T>& s,KRYLOV_VECTOR_BASE<T>& r,
        KRYLOV_VECTOR_BASE<T>& k,KRYLOV_VECTOR_BASE<T>& z,const T tolerance,const int min_iterations,const int max_iterations) PHYSBAM_OVERRIDE;
//#####################################################################
};
}
#endif

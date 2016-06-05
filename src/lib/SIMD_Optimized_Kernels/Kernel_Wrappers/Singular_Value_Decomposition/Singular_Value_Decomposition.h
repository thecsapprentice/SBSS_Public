//#####################################################################
//  Copyright (c) 2011-2013 Nathan Mitchell, Eftychios Sifakis.
//  This file is covered by the FreeBSD license. Please refer to the 
//  license.txt file for more information.
//#####################################################################

//#####################################################################
// Function Singular_Value_Decomposition
//#####################################################################
#ifndef __Singular_Value_Decomposition_Wrapper__
#define __Singular_Value_Decomposition_Wrapper__

#include "../../Kernels/Singular_Value_Decomposition/Singular_Value_Decomposition.h"

namespace PhysBAM{
namespace{

template<class T,int d> void
Singular_Value_Decomposition(const MATRIX<T,d>& A, MATRIX<T,d>& U, DIAGONAL_MATRIX<T,d>& Sigma, MATRIX<T,d>& V)
{
    A.Fast_Singular_Value_Decomposition(U,Sigma,V);
}

template<> void
Singular_Value_Decomposition(const MATRIX<float,3>& A, MATRIX<float,3>& U, DIAGONAL_MATRIX<float,3>& Sigma, MATRIX<float,3>& V)
{
    typedef const float (&refConstMat3)[9];
    typedef float (&refMat3)[9];
    typedef float (&refDMat3)[3];

    ::Singular_Value_Decomposition<float,float,int>(
        reinterpret_cast<refConstMat3>(A),
        reinterpret_cast<refMat3>(U),
        reinterpret_cast<refDMat3>(Sigma),
        reinterpret_cast<refMat3>(V)
    );
}

}
}

#endif

//#####################################################################
//  Copyright (c) 2011-2013 Nathan Mitchell, Eftychios Sifakis.
//  This file is covered by the FreeBSD license. Please refer to the 
//  license.txt file for more information.
//#####################################################################

//#####################################################################
// Function Penalty_Measure_Gradient
//#####################################################################
#ifndef __Penalty_Measure_Gradient_Wrapper__
#define __Penalty_Measure_Gradient_Wrapper__

#include "../../Kernels/Penalty_Measure_Gradient/Penalty_Measure_Gradient.h"

namespace PhysBAM{
namespace{

template<class T_MATERIAL> void
Penalty_Measure_Gradient(const DIAGONAL_MATRIX<typename T_MATERIAL::SCALAR,T_MATERIAL::dim>& Sigma,
    DIAGONAL_MATRIX<typename T_MATERIAL::SCALAR,T_MATERIAL::dim>& Q_hat)
{
    Q_hat=MATERIAL_MODEL<T_MATERIAL>::Q_hat(Sigma);
}

/*
template<> void
Penalty_Measure_Gradient<NEOHOOKEAN<float,3> >(const DIAGONAL_MATRIX<float,3>& Sigma,DIAGONAL_MATRIX<float,3>& Q_hat)
{
    typedef const float (&refConstDMat3)[3];
    typedef float (&refDMat3)[3];

    ::Penalty_Measure_Gradient<NEOHOOKEAN_TAG,float,float,int>::Run(reinterpret_cast<refConstDMat3>(Sigma),reinterpret_cast<refDMat3>(Q_hat));
}
*/

}
}

#endif

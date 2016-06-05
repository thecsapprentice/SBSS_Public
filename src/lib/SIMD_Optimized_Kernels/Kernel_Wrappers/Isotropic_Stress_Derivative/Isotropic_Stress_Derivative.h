//#####################################################################
//  Copyright (c) 2011-2013 Nathan Mitchell, Eftychios Sifakis.
//  This file is covered by the FreeBSD license. Please refer to the 
//  license.txt file for more information.
//#####################################################################

//#####################################################################
// Function Isotropic_Stress_Derivative
//#####################################################################
#ifndef __Isotropic_Stress_Derivative_Wrapper__
#define __Isotropic_Stress_Derivative_Wrapper__

#include "../../Kernels/Isotropic_Stress_Derivative/Isotropic_Stress_Derivative.h"

namespace PhysBAM{
namespace{

template<class T_MATERIAL> void
Isotropic_Stress_Derivative(ROTATED_STRESS_DERIVATIVE<typename T_MATERIAL::SCALAR,T_MATERIAL::dim>& dPdF,
    const DIAGONAL_MATRIX<typename T_MATERIAL::SCALAR,T_MATERIAL::dim>& Sigma,const typename T_MATERIAL::SCALAR p,
    const typename T_MATERIAL::SCALAR mu,const typename T_MATERIAL::SCALAR kappa,
    const typename T_MATERIAL::SCALAR alpha,const bool apply_definiteness_fix)
{
    dPdF=MATERIAL_MODEL<T_MATERIAL>::Isotropic_Stress_Derivative(Sigma,p,mu,kappa,alpha,apply_definiteness_fix);
}

template<class T_MATERIAL> void
Isotropic_Stress_Derivative(ROTATED_STRESS_DERIVATIVE<typename T_MATERIAL::SCALAR,T_MATERIAL::dim>& dPdF,
    const DIAGONAL_MATRIX<typename T_MATERIAL::SCALAR,T_MATERIAL::dim>& Sigma,const typename T_MATERIAL::SCALAR p,
    const typename T_MATERIAL::SCALAR mu_10,
    const typename T_MATERIAL::SCALAR mu_01,
    const typename T_MATERIAL::SCALAR kappa,
    const typename T_MATERIAL::SCALAR alpha,const bool apply_definiteness_fix)
{
    dPdF=MATERIAL_MODEL<T_MATERIAL>::Isotropic_Stress_Derivative(Sigma,p,mu_10,mu_01,kappa,alpha,apply_definiteness_fix);
}

/*
template<> void
Isotropic_Stress_Derivative<NEOHOOKEAN<float,3> >(ROTATED_STRESS_DERIVATIVE<float,3>& dPdF,
    const DIAGONAL_MATRIX<float,3>& Sigma,const float p,
    const float mu,const float kappa,
    const float alpha,const bool apply_definiteness_fix)
{
    typedef float (&refDpDf)[12];
    typedef const float (&refConstDMat3)[3];

    ::Isotropic_Stress_Derivative<NEOHOOKEAN_TAG,float,float,int>::Run(
        reinterpret_cast<refDpDf>(dPdF),
        reinterpret_cast<refConstDMat3>(Sigma),
        p,mu,kappa,alpha,apply_definiteness_fix
    );
}
*/

}
}

#endif

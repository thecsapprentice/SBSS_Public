//#####################################################################
//  Copyright (c) 2011-2013 Nathan Mitchell, Eftychios Sifakis.
//  This file is covered by the FreeBSD license. Please refer to the 
//  license.txt file for more information.
//#####################################################################


#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Arrays_Computations/MAGNITUDE.h>
#include <PhysBAM_Tools/Matrices/MATRIX_2X2.h>
#include <PhysBAM_Tools/Matrices/MATRIX_3X3.h>
#include <PhysBAM_Tools/Matrices/SYMMETRIC_MATRIX_3X3.h>

#include "NEOHOOKEAN.h"
#include "COROTATED.h"
#include "MATERIAL_MODEL.h"
#include "ROTATED_STRESS_DERIVATIVE.h"

using namespace PhysBAM;

template<class T>
void Isotropic_Stress_Derivative_Neohookean_Reference(T dPdF[12], const T Sigma[3], const T p, const T mu, const T kappa, const T alpha, const bool apply_definiteness_fix)
{
    ROTATED_STRESS_DERIVATIVE<T,3>& rdPdF=*(ROTATED_STRESS_DERIVATIVE<T,3>*)(dPdF);
    const DIAGONAL_MATRIX<T,3>& mSigma=*(const DIAGONAL_MATRIX<T,3>*)(Sigma);    
    rdPdF=MATERIAL_MODEL<NEOHOOKEAN<T,3> >::Isotropic_Stress_Derivative(mSigma,p,mu,kappa,alpha,apply_definiteness_fix);
}

template<class T>
bool Isotropic_Stress_Derivative_Neohookean_Compare(const T dPdF[12], const T dPdF_reference[12])
{
    ARRAY_VIEW<const T> adPdF(12,dPdF);
    ARRAY_VIEW<const T> adPdF_reference(12,dPdF_reference);

    std::cout<<"Computed dPdF : "<<adPdF<<std::endl;
    std::cout<<"Reference dPdF :"<<adPdF_reference<<std::endl;
    ARRAY<T> difference(adPdF-adPdF_reference);
    std::cout<<"Difference = "<<ARRAYS_COMPUTATIONS::Maximum_Magnitude(difference)<<std::endl;

    if( ARRAYS_COMPUTATIONS::Maximum_Magnitude(difference) < 0.00001 )
        return true;
    else
        return false;
}

template void Isotropic_Stress_Derivative_Neohookean_Reference(float dPdF[12], const float Sigma[3], const float p, const float mu, const float kappa, const float alpha, const bool apply_definiteness_fix);
template bool Isotropic_Stress_Derivative_Neohookean_Compare(const float dPdF[12], const float dPdF_reference[12]);

template<class T>
void Isotropic_Stress_Derivative_Corotated_Reference(T dPdF[12], const T Sigma[3], const T p, const T mu, const T kappa, const T alpha, const bool apply_definiteness_fix)
{
    ROTATED_STRESS_DERIVATIVE<T,3>& rdPdF=*(ROTATED_STRESS_DERIVATIVE<T,3>*)(dPdF);
    const DIAGONAL_MATRIX<T,3>& mSigma=*(const DIAGONAL_MATRIX<T,3>*)(Sigma);    
    rdPdF=MATERIAL_MODEL<COROTATED<T,3> >::Isotropic_Stress_Derivative(mSigma,p,mu,kappa,alpha,apply_definiteness_fix);
}

template<class T>
bool Isotropic_Stress_Derivative_Corotated_Compare(const T dPdF[12], const T dPdF_reference[12])
{
    ARRAY_VIEW<const T> adPdF(12,dPdF);
    ARRAY_VIEW<const T> adPdF_reference(12,dPdF_reference);

    std::cout<<"Computed dPdF : "<<adPdF<<std::endl;
    std::cout<<"Reference dPdF :"<<adPdF_reference<<std::endl;
    ARRAY<T> difference(adPdF-adPdF_reference);
    std::cout<<"Difference = "<<ARRAYS_COMPUTATIONS::Maximum_Magnitude(difference)<<std::endl;

    if( ARRAYS_COMPUTATIONS::Maximum_Magnitude(difference) < 0.00001 )
        return true;
    else
        return false;
}

template void Isotropic_Stress_Derivative_Corotated_Reference(float dPdF[12], const float Sigma[3], const float p, const float mu, const float kappa, const float alpha, const bool apply_definiteness_fix);
template bool Isotropic_Stress_Derivative_Corotated_Compare(const float dPdF[12], const float dPdF_reference[12]);

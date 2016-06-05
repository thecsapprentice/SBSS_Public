//#####################################################################
//  Copyright (c) 2011-2013 Nathan Mitchell, Eftychios Sifakis.
//  This file is covered by the FreeBSD license. Please refer to the 
//  license.txt file for more information.
//#####################################################################


template<class T>
void Isotropic_Stress_Derivative_Neohookean_Reference(T dPdF[12], const T Sigma[3], const T p, const T mu, const T kappa, const T alpha, const bool apply_definiteness_fix);

template<class T>
bool Isotropic_Stress_Derivative_Neohookean_Compare(const T dPdF[12], const T dPdF_reference[12]);


template<class T>
void Isotropic_Stress_Derivative_Corotated_Reference(T dPdF[12], const T Sigma[3], const T p, const T mu, const T kappa, const T alpha, const bool apply_definiteness_fix);

template<class T>
bool Isotropic_Stress_Derivative_Corotated_Compare(const T dPdF[12], const T dPdF_reference[12]);

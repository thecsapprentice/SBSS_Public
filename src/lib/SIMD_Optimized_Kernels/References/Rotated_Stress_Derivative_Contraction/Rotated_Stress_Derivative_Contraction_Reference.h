//#####################################################################
//  Copyright (c) 2011-2013 Nathan Mitchell, Eftychios Sifakis.
//  This file is covered by the FreeBSD license. Please refer to the 
//  license.txt file for more information.
//#####################################################################


template<class T>
void Rotated_Stress_Derivative_Contraction_Reference(const T dPdF[12], const T dF_Hat[9], T dP_Hat[9]);

template<class T>
bool Rotated_Stress_Derivative_Contraction_Compare(const T dP_Hat[9], const T dP_Hat_reference[9]);

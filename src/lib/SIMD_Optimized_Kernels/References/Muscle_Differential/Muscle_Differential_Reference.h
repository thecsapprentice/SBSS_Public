//#####################################################################
//  Copyright (c) 2011-2013 Nathan Mitchell, Eftychios Sifakis.
//  This file is covered by the FreeBSD license. Please refer to the 
//  license.txt file for more information.
//#####################################################################


template<class T>
void Muscle_Differential_Reference(T dP_fiber[9], const T dF[9], const T fiber[3],  const T Ffiber[3], const T c1, const T c2);
template<class T>
bool Muscle_Differential_Compare(const T dP_fiber[9], const T dP_fiber_reference[9]);

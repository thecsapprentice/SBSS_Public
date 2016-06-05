//#####################################################################
//  Copyright (c) 2011-2013 Nathan Mitchell, Eftychios Sifakis.
//  This file is covered by the FreeBSD license. Please refer to the 
//  license.txt file for more information.
//#####################################################################


template<class T>
void Symmetric_Definite_Projection_Reference(const T A[6], T Apd[6]);

template<class T>
bool Symmetric_Definite_Projection_Compare(const T Apd[6], const T Apd_reference[6]);

//#####################################################################
//  Copyright (c) 2011-2013 Nathan Mitchell, Eftychios Sifakis.
//  This file is covered by the FreeBSD license. Please refer to the 
//  license.txt file for more information.
//#####################################################################


template<class T>
void Unweighted_Gradient_Reference(const T u[3][8], T F[9], const T one_over_h);

template<class T>
bool Unweighted_Gradient_Compare(const T F[9], const T F_reference[9]);

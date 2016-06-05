//#####################################################################
//  Copyright (c) 2011-2013 Nathan Mitchell, Eftychios Sifakis.
//  This file is covered by the FreeBSD license. Please refer to the 
//  license.txt file for more information.
//#####################################################################


template<class T>
void Unweighted_Accumulation_Reference(T u[3][8], const T F[9], const T one_over_h,const T scale);

template<class T>
bool Unweighted_Accumulation_Compare(const T u[3][8], const T u_reference[3][8]);

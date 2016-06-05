//#####################################################################
//  Copyright (c) 2011-2013 Nathan Mitchell, Eftychios Sifakis.
//  This file is covered by the FreeBSD license. Please refer to the 
//  license.txt file for more information.
//#####################################################################


template<class T>
void Force_Stabilization_Reference(const T Du[3][8], const T constant,/* const T h, const T stabilization_factor, */ T dH[3][8] );

template<class T>
bool Force_Stabilization_Compare(const T dH[3][8],  const T dH_reference[3][8]);

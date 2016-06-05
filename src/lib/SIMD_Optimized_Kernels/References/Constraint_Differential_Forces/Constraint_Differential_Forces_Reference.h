//#####################################################################
//  Copyright (c) 2011-2013 Nathan Mitchell, Eftychios Sifakis.
//  This file is covered by the FreeBSD license. Please refer to the 
//  license.txt file for more information.
//#####################################################################


template<class T>
void Constraint_Differential_Forces_Reference(T df[3][8], const T du[3][8],
    const T W[3], const T scale );

template<class T>
bool Constraint_Differential_Forces_Compare(const T df[3][8], const T df_reference[3][8]);

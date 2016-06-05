//#####################################################################
//  Copyright (c) 2011-2013 Nathan Mitchell, Eftychios Sifakis.
//  This file is covered by the FreeBSD license. Please refer to the 
//  license.txt file for more information.
//#####################################################################


template<class T>
void Muscle_Differential_Forces_Reference(T df[3][8], const T du[3][8], const T fiber[3],  const T Ffiber[3], const T c1, const T c2, const T one_over_h, const T cell_volume);
template<class T>
bool Muscle_Differential_Forces_Compare(const T df[3][8], const T df_reference[3][8]);

//#####################################################################
//  Copyright (c) 2011-2013 Nathan Mitchell, Eftychios Sifakis.
//  This file is covered by the FreeBSD license. Please refer to the 
//  license.txt file for more information.
//#####################################################################


template<class T_RAW,class T_DATA=void,class I_DATA=void>
void Matrix_Transpose_Times(const T_DATA (&A)[9], const T_DATA (&B)[9], T_DATA (&C)[9]);

//#####################################################################
//  Copyright (c) 2011-2013 Nathan Mitchell, Eftychios Sifakis.
//  This file is covered by the FreeBSD license. Please refer to the 
//  license.txt file for more information.
//#####################################################################


template<class T_MATERIAL,class Tw,class T_DATA=void,class I_DATA=void> 
struct Pressure_Force
{
    static void Run(T_DATA (&q), const T_DATA (&Sigma)[3], const T_DATA (&p), const T_DATA (&alpha), const T_DATA (&alpha_squared_over_kappa));
};

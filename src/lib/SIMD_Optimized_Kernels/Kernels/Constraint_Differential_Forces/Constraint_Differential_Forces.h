//#####################################################################
//  Copyright (c) 2011-2013 Nathan Mitchell, Eftychios Sifakis.
//  This file is covered by the FreeBSD license. Please refer to the 
//  license.txt file for more information.
//#####################################################################


template<class T_RAW,class T_DATA=void,class I_DATA=void>
    void Constraint_Differential_Forces(T_DATA (&df)[3][8], const T_DATA (&u)[3][8],
                                        const T_DATA (&W)[3], const T_DATA (&scale) );

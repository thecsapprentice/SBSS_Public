//#####################################################################
//  Copyright (c) 2011-2013 Nathan Mitchell, Eftychios Sifakis.
//  This file is covered by the FreeBSD license. Please refer to the 
//  license.txt file for more information.
//#####################################################################


template<class T_RAW,class T_DATA=void,class I_DATA=void>
    void Rotated_Stress_Derivative_Contraction(const T_DATA (&dPdF)[12], const T_DATA (&dF_Hat)[9], T_DATA (&dP_Hat)[9]);

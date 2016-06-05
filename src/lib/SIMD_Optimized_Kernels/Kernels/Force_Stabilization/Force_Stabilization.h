//#####################################################################
//  Copyright (c) 2011-2013 Nathan Mitchell, Eftychios Sifakis.
//  This file is covered by the FreeBSD license. Please refer to the 
//  license.txt file for more information.
//#####################################################################


template<class T_RAW,class T_DATA=void,class I_DATA=void>
    void Force_Stabilization(const T_DATA (&Du)[3][8], const T_DATA (&constant), T_DATA (&dH)[3][8] );

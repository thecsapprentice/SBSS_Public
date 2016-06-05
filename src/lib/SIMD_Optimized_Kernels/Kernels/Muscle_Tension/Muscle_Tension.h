//#####################################################################
//  Copyright (c) 2011-2013 Nathan Mitchell, Eftychios Sifakis.
//  This file is covered by the FreeBSD license. Please refer to the 
//  license.txt file for more information.
//#####################################################################


template<class T_RAW,class T_DATA=void,class I_DATA=void>
void Tension(T_DATA (&tension), const T_DATA (&stretch), const T_DATA (&activation), const T_DATA (&density), const T_DATA (&fiber_max_stress));

template<class T_RAW,class T_DATA=void,class I_DATA=void>
void Tension_Derivative(T_DATA (&tension_derivative), const T_DATA (&stretch), const T_DATA (&activation), const T_DATA (&density), const T_DATA (&fiber_max_stress));


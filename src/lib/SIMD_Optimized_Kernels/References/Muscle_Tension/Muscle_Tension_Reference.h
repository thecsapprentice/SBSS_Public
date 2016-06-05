//#####################################################################
//  Copyright (c) 2011-2013 Nathan Mitchell, Eftychios Sifakis.
//  This file is covered by the FreeBSD license. Please refer to the 
//  license.txt file for more information.
//#####################################################################


template<class T>
void Tension_Reference(T& tension, const T stretch, const T activation, const T density, const T fiber_max_stress);

template<class T>
void Tension_Derivative_Reference(T& tension_derivative, const T stretch, const T activation, const T density, const T fiber_max_stress);

template<class T>
bool Tension_Compare(const T tension, const T tension_reference);

template<class T>
bool Tension_Derivative_Compare(const T tension_derivative, const T tension_derivative_reference);

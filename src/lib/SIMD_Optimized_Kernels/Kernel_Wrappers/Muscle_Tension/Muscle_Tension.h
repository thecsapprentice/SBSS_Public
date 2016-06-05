//#####################################################################
//  Copyright (c) 2011-2013 Nathan Mitchell, Eftychios Sifakis.
//  This file is covered by the FreeBSD license. Please refer to the 
//  license.txt file for more information.
//#####################################################################

//#####################################################################
// Function Central_Gradient
//#####################################################################
#ifndef __Muscle_Tension_Wrapper__
#define __Muscle_Tension_Wrapper__

#include "../../Kernels/Muscle_Tension/Muscle_Tension.h"

namespace PhysBAM{

 template <class T> void
    K_Tension(T& tension, const T& stretch, const T& activation, const T& density, const T& fiber_max_stress)
{
    ::Tension<float,float,int>(tension, stretch, activation, density, fiber_max_stress);
}


 template<class T> void
    K_Tension_Derivative(T& tension, const T& stretch, const T& activation, const T& density, const T& fiber_max_stress)
{
    ::Tension_Derivative<float,float,int>(tension, stretch, activation, density, fiber_max_stress);
}



}

#endif

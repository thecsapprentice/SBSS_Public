//#####################################################################
//  Copyright (c) 2011-2013 Nathan Mitchell, Eftychios Sifakis.
//  This file is covered by the FreeBSD license. Please refer to the 
//  license.txt file for more information.
//#####################################################################


#ifndef SUBROUTINE_Muscle_Tension
#include <assert.h>
#include "KernelCommon.h"
#else
namespace {
#endif

namespace{

const float Sfiber_p1=.05f;
const float Sfiber_p2=6.6f;
const float Sfiber_cutoff=1.4f;
const float Scutoff_scaled=Sfiber_p2*(Sfiber_cutoff-1.f);
const float Sfiber_p3=Sfiber_p1*Sfiber_p2*(exp(Scutoff_scaled)-1.f);
const float Sfiber_p4=Sfiber_p1*(exp(Scutoff_scaled)*(1.f-Sfiber_p2*Sfiber_cutoff)+Sfiber_p2-1.f);
const float Sscale=25.f/6.f;

 BUILD_CONSTANT(fiber_p1,Sfiber_p1);
 BUILD_CONSTANT(fiber_p2,Sfiber_p2);
 BUILD_CONSTANT(fiber_cutoff,Sfiber_cutoff);
 BUILD_CONSTANT(fiber_p3,Sfiber_p3);
 BUILD_CONSTANT(fiber_p4,Sfiber_p4);
 BUILD_CONSTANT(scale,Sscale);
 BUILD_CONSTANT(one,1.0f);
 BUILD_CONSTANT(two,2.0f);
 BUILD_CONSTANT(two_fifths,0.4f);
 BUILD_CONSTANT(three_fifths,0.6f);
    
}

template<class Tw,class T_DATA=void,class I_DATA=void>
#ifdef SUBROUTINE_Muscle_Tension
inline
#endif
void Tension(T_DATA (&tension), const T_DATA (&stretch), const T_DATA (&activation), const T_DATA (&density), const T_DATA (&fiber_max_stress))
{
    typedef Number<Tw> Tn;

    Tn rfiber_p1,rfiber_p2,rfiber_cutoff,rfiber_p3,rfiber_p4,rscale,rone,rtwo,rtwo_fifths,rthree_fifths;
    rfiber_p1.Load_Aligned(fiber_p1);rfiber_p2.Load_Aligned(fiber_p2);rfiber_cutoff.Load_Aligned(fiber_cutoff);rfiber_p3.Load_Aligned(fiber_p3);rfiber_p4.Load_Aligned(fiber_p4);
    rscale.Load_Aligned(scale);rone.Load_Aligned(one);rtwo.Load_Aligned(two);rtwo_fifths.Load_Aligned(two_fifths);rthree_fifths.Load_Aligned(three_fifths);

    Tn rstretch,ractivation,rdensity,rfiber_max_stress;
    rstretch.Load(stretch);ractivation.Load(activation);rdensity.Load(density);rfiber_max_stress.Load(fiber_max_stress);
    
    Tn rstrain=rstretch-rone,rstrain_abs=rstrain.abs(),ractive_tension,rpassive_tension;
    
    Tn rpassive_tension1=rfiber_p1*((rfiber_p2*rstrain).exp()-rfiber_p2*rstrain-rone);
    Tn rpassive_tension2=rfiber_p3*rstretch+rfiber_p4;
    rpassive_tension=blend(rone<rstretch,rpassive_tension,rpassive_tension1);
    rpassive_tension=blend(rfiber_cutoff<rstretch,rpassive_tension,rpassive_tension2);

    Tn ractive_tension1=rtwo*rscale*ractivation*rdensity*(rstrain_abs-rthree_fifths)*(rstrain_abs-rthree_fifths);
    Tn ractive_tension2=ractivation*rdensity*(rone-rscale*rstrain*rstrain);
    ractive_tension=blend(rstrain_abs<rthree_fifths,ractive_tension,ractive_tension1);
    ractive_tension=blend(rstrain_abs<rtwo_fifths,ractive_tension,ractive_tension2);

    Tn rtension=rfiber_max_stress*(ractive_tension+rpassive_tension);
    Store(tension,rtension);
}

template<class Tw,class T_DATA=void,class I_DATA=void>
#ifdef SUBROUTINE_Muscle_Tension
inline
#endif
void Tension_Derivative(T_DATA (&tension_derivative), const T_DATA (&stretch), const T_DATA (&activation), const T_DATA (&density), const T_DATA (&fiber_max_stress))
{
    typedef Number<Tw> Tn;

    Tn rfiber_p1,rfiber_p2,rfiber_cutoff,rfiber_p3,rscale,rone,rtwo,rtwo_fifths,rthree_fifths;
    rfiber_p1.Load_Aligned(fiber_p1);rfiber_p2.Load_Aligned(fiber_p2);rfiber_cutoff.Load_Aligned(fiber_cutoff);rfiber_p3.Load_Aligned(fiber_p3);
    rscale.Load_Aligned(scale);rone.Load_Aligned(one);rtwo.Load_Aligned(two);rtwo_fifths.Load_Aligned(two_fifths);rthree_fifths.Load_Aligned(three_fifths);

    Tn rstretch,ractivation,rdensity,rfiber_max_stress;
    rstretch.Load(stretch);ractivation.Load(activation);rdensity.Load(density);rfiber_max_stress.Load(fiber_max_stress);
    
    Tn rstrain=rstretch-rone,rstrain_abs=rstrain.abs(),ractive_tension_d,rpassive_tension_d;
    
    Tn rpassive_tension1=rfiber_p1*rfiber_p2*((rfiber_p2*rstrain).exp()-rone);
    Tn rpassive_tension2=rfiber_p3;
    rpassive_tension_d=blend(rone<rstretch,rpassive_tension_d,rpassive_tension1);
    rpassive_tension_d=blend(rfiber_cutoff<rstretch,rpassive_tension_d,rpassive_tension2);

    Tn ractive_tension1=(rtwo+rtwo)*rscale*ractivation*rdensity*(rstrain-rstrain.sign()*rthree_fifths);
    Tn ractive_tension2=ractive_tension_d-(rtwo*rscale*ractivation*rdensity*rstrain);
    ractive_tension_d=blend(rstrain_abs<rthree_fifths,ractive_tension_d,ractive_tension1);
    ractive_tension_d=blend(rstrain_abs<rtwo_fifths,ractive_tension_d,ractive_tension2);

    Tn rtension_d=rfiber_max_stress*(ractive_tension_d+rpassive_tension_d);
    Store(tension_derivative,rtension_d);
}


#ifndef SUBROUTINE_Muscle_Tension
#define INSTANCE_KERNEL_Tension(WIDTH) WIDETYPE(float,WIDTH) (&tension), const WIDETYPE(float,WIDTH) (&stretch), const WIDETYPE(float,WIDTH) (&activation), const WIDETYPE(float,WIDTH) (&density), const WIDETYPE(float,WIDTH) (&fiber_max_stress)
#define INSTANCE_KERNEL_Tension_Derivative(WIDTH) WIDETYPE(float,WIDTH) (&tension_derivative), const WIDETYPE(float,WIDTH) (&stretch), const WIDETYPE(float,WIDTH) (&activation), const WIDETYPE(float,WIDTH) (&density), const WIDETYPE(float,WIDTH) (&fiber_max_stress)
INSTANCE_KERNEL(Tension);
INSTANCE_KERNEL(Tension_Derivative);
#undef INSTANCE_KERNEL_Tension
#undef INSTANCE_KERNEL_Tension_Derivative
#else
}
#endif


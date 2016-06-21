//#####################################################################
//  Copyright (c) 2011-2013 Nathan Mitchell, Eftychios Sifakis.
//  This file is covered by the FreeBSD license. Please refer to the 
//  license.txt file for more information.
//#####################################################################


#ifndef SUBROUTINE_Pressure_Force
#include <assert.h>
#include "KernelCommon.h"
#else
namespace {
#endif

#ifdef FORCE_INLINE
#include "../Volume_Preservation_Deviation/Volume_Preservation_Deviation.h"
#else
#define SUBROUTINE_Volume_Preservation_Deviation
#include "../Volume_Preservation_Deviation/Volume_Preservation_Deviation.cpp"
#undef SUBROUTINE_Augmented_Rotated_Stress_Derivative
#endif

template<class T_MATERIAL,class Tw,class T_DATA=void,class I_DATA=void> 
struct Pressure_Force
{
static 
#ifdef SUBROUTINE_Pressure_Force
inline
#endif
void Run(T_DATA (&q), const T_DATA (&Sigma)[3], const T_DATA (&p), const T_DATA (&alpha), const T_DATA (&alpha_squared_over_kappa))
{
    typedef Number<Tw> Tn;
#ifdef USE_NONMIXED_FORMULAS
    Tn rp;
    Tn rq;
    Tn rAlphaSqrd;
    rp.Load( p );
    rAlphaSqrd.Load( alpha_squared_over_kappa );
    rq = rq - (rAlphaSqrd*rp);
    Store( q, rq );
#else
    T_DATA M;
    Volume_Preservation_Deviation<T_MATERIAL,Tw,T_DATA,I_DATA>::Run(M, Sigma);
    
    Tn rp;
    Tn rq;
    Tn rM;
    Tn rAlpha;
    Tn rAlphaSqrd;

    rp.Load( p );
    rM.Load( M );
    rAlpha.Load( alpha );
    rAlphaSqrd.Load( alpha_squared_over_kappa );
    
    rq = rAlphaSqrd*rp - rAlpha*rM;

    Store( q, rq );
#endif
}
};

#ifndef SUBROUTINE_Pressure_Force
#define INSTANCE_KERNEL_Pressure_Force(WIDTH) WIDETYPE(float,WIDTH) (&q),    const WIDETYPE(float,WIDTH) (&Sigma)[3],    const WIDETYPE(float,WIDTH) (&p),    const WIDETYPE(float,WIDTH) (&alpha),    const WIDETYPE(float,WIDTH) (&alpha_squared_over_kappa)
INSTANCE_KERNEL_MATERIAL(Pressure_Force,COROTATED_TAG);
INSTANCE_KERNEL_MATERIAL(Pressure_Force,NEOHOOKEAN_TAG);
#if defined(BIPHASIC_SUPPORT)
INSTANCE_KERNEL_MATERIAL(Pressure_Force,BIPHASIC_TAG);
#endif
#undef INSTANCE_KERNEL_Pressure_Force
#else
}
#endif

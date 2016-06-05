//#####################################################################
//  Copyright (c) 2011-2013 Nathan Mitchell, Eftychios Sifakis.
//  This file is covered by the FreeBSD license. Please refer to the 
//  license.txt file for more information.
//#####################################################################


#ifndef SUBROUTINE_Pressure_Force_Differential
#include <assert.h>
#include "KernelCommon.h"
#else
namespace {
#endif

template<class Tw,class T_DATA=void,class I_DATA=void>
#ifdef SUBROUTINE_Pressure_Force_Differential
inline
#endif
void Pressure_Force_Differential(T_DATA (&dq), const T_DATA (&Q_hat)[3], const T_DATA (&dF_hat)[9],
                                 const T_DATA (&dp), const T_DATA (&alpha),
                                 const T_DATA (&alpha_squared_over_kappa))
{
    typedef enum { x11=0,x21,x31,x12,x22,x32,x13,x23,x33} Matrix_Entry;

    typedef Number<Tw> Tn;

    Tn rdq;

    Tn dF_hat_11;
    Tn dF_hat_22;
    Tn dF_hat_33;

    Tn rAlpha;
    Tn rAlphaSqrd;
    Tn rdp;
    Tn rQ_hat1,rQ_hat2,rQ_hat3;
#ifdef USE_NONMIXED_FORMULAS
    rAlphaSqrd.Load(alpha_squared_over_kappa);
    rdp.Load(dp);
    rdq = rdq - (rAlphaSqrd*rdp);
    Store(dq, rdq);
#else
    rAlpha.Load(alpha);
    rAlphaSqrd.Load(alpha_squared_over_kappa);
    rdp.Load(dp);

    rQ_hat1.Load(Q_hat[0]);
    dF_hat_11.Load(dF_hat[x11]);

    rQ_hat2.Load(Q_hat[1]);   
    dF_hat_22.Load(dF_hat[x22]);

    rQ_hat3.Load(Q_hat[2]);   
    dF_hat_33.Load(dF_hat[x33]);


    rdq = (dF_hat_11 * rQ_hat1) + (dF_hat_22 * rQ_hat2)  + (dF_hat_33 * rQ_hat3);

    rdq = (rAlphaSqrd*rdp) -  (rdq*rAlpha);

    Store(dq, rdq);
#endif
}

#ifndef SUBROUTINE_Pressure_Force_Differential
#define INSTANCE_KERNEL_Pressure_Force_Differential(WIDTH) WIDETYPE(float,WIDTH) (&dq), const WIDETYPE(float,WIDTH) (&Q_hat)[3],  const WIDETYPE(float,WIDTH) (&dF_hat)[9],  const WIDETYPE(float,WIDTH) (&dp), const WIDETYPE(float,WIDTH) (&alpha), const WIDETYPE(float,WIDTH) (&alpha_squared_over_kappa)
INSTANCE_KERNEL(Pressure_Force_Differential);
#undef INSTANCE_KERNEL_Pressure_Force_Differential
#else
}
#endif

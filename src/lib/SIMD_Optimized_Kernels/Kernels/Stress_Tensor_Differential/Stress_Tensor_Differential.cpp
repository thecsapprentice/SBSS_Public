//#####################################################################
//  Copyright (c) 2011-2013 Nathan Mitchell, Eftychios Sifakis.
//  This file is covered by the FreeBSD license. Please refer to the 
//  license.txt file for more information.
//#####################################################################


#ifndef SUBROUTINE_Stress_Tensor_Differential
#include <assert.h>
#include "KernelCommon.h"
#else
namespace {
#endif

template<class Tw,class T_DATA=void,class I_DATA=void>
#ifdef SUBROUTINE_Stress_Tensor_Differential
inline
#endif
    void Stress_Tensor_Differential(T_DATA (&dP_hat)[9], const T_DATA (&dPdF)[12], const T_DATA (&dF_hat)[9],
                                    const T_DATA (&Q_hat)[3], const T_DATA (&dp), const T_DATA (&alpha))
{
    typedef enum { x11=0,x21,x31,x12,x22,x32,x13,x23,x33} Matrix_Entry;
    typedef enum { x1111=0,x1122,x1133,x2222,x2233,x3333,x1212,x1221,x1313,x1331,x2323,x2332 } RSD_Entry;

    typedef Number<Tw> Tn;

    Tn dP_hat_11;
    Tn dP_hat_21;
    Tn dP_hat_31;
    Tn dP_hat_12;
    Tn dP_hat_22;
    Tn dP_hat_32;
    Tn dP_hat_13;
    Tn dP_hat_23;
    Tn dP_hat_33;

    Tn dF_hat_11;
    Tn dF_hat_21;
    Tn dF_hat_31;
    Tn dF_hat_12;
    Tn dF_hat_22;
    Tn dF_hat_32;
    Tn dF_hat_13;
    Tn dF_hat_23;
    Tn dF_hat_33;

    Tn a1111;
    Tn a1122;
    Tn a1133;
    Tn a2222;
    Tn a2233;
    Tn a3333;
    Tn a1212;
    Tn a1221;
    Tn a1313;
    Tn a1331;
    Tn a2323;
    Tn a2332;
    
    dF_hat_11.Load(dF_hat[x11]);
    dF_hat_21.Load(dF_hat[x21]);
    dF_hat_31.Load(dF_hat[x31]);
    dF_hat_12.Load(dF_hat[x12]);
    dF_hat_22.Load(dF_hat[x22]);
    dF_hat_32.Load(dF_hat[x32]);
    dF_hat_13.Load(dF_hat[x13]);
    dF_hat_23.Load(dF_hat[x23]);
    dF_hat_33.Load(dF_hat[x33]);

    a1111.Load(dPdF[x1111]);
    a1122.Load(dPdF[x1122]);
    a1133.Load(dPdF[x1133]);
    a2222.Load(dPdF[x2222]);
    a2233.Load(dPdF[x2233]);
    a3333.Load(dPdF[x3333]);
    a1212.Load(dPdF[x1212]);
    a1221.Load(dPdF[x1221]);
    a1313.Load(dPdF[x1313]);
    a1331.Load(dPdF[x1331]);
    a2323.Load(dPdF[x2323]);
    a2332.Load(dPdF[x2332]);

    dP_hat_11=a1111*dF_hat_11+a1122*dF_hat_22+a1133*dF_hat_33;
    dP_hat_22=a1122*dF_hat_11+a2222*dF_hat_22+a2233*dF_hat_33;
    dP_hat_33=a1133*dF_hat_11+a2233*dF_hat_22+a3333*dF_hat_33;
    dP_hat_12=a1212*dF_hat_12+a1221*dF_hat_21;
    dP_hat_21=a1221*dF_hat_12+a1212*dF_hat_21;
    dP_hat_13=a1313*dF_hat_13+a1331*dF_hat_31;
    dP_hat_31=a1331*dF_hat_13+a1313*dF_hat_31;
    dP_hat_23=a2323*dF_hat_23+a2332*dF_hat_32;
    dP_hat_32=a2332*dF_hat_23+a2323*dF_hat_32;
    
    Store(dP_hat[x12], dP_hat_12);
    Store(dP_hat[x13], dP_hat_13);
    Store(dP_hat[x21], dP_hat_21);
    Store(dP_hat[x23], dP_hat_23);
    Store(dP_hat[x31], dP_hat_31);
    Store(dP_hat[x32], dP_hat_32);

#ifdef USE_NONMIXED_FORMULAS
    // Don't do this part
#else
    Tn rAlpha;
    Tn rdp;
    Tn rQ_hat;

    rAlpha.Load(alpha);
    rdp.Load(dp);

    rQ_hat.Load(Q_hat[0]);   
    dP_hat_11 = dP_hat_11 + (rAlpha*rdp*rQ_hat);

    rQ_hat.Load(Q_hat[1]);   
    dP_hat_22 = dP_hat_22 + (rAlpha*rdp*rQ_hat);

    rQ_hat.Load(Q_hat[2]);   
    dP_hat_33 = dP_hat_33 + (rAlpha*rdp*rQ_hat);  
#endif

    Store(dP_hat[x11], dP_hat_11);
    Store(dP_hat[x22], dP_hat_22);
    Store(dP_hat[x33], dP_hat_33);

}

#ifndef SUBROUTINE_Stress_Tensor_Differential
#define INSTANCE_KERNEL_Stress_Tensor_Differential(WIDTH) WIDETYPE(float,WIDTH) (&dP_hat)[9], const WIDETYPE(float,WIDTH) (&dPdF)[12], const WIDETYPE(float,WIDTH) (&dF_hat)[9],const WIDETYPE(float,WIDTH) (&Q_hat)[3], const WIDETYPE(float,WIDTH) (&p), const WIDETYPE(float,WIDTH) (&alpha)
INSTANCE_KERNEL(Stress_Tensor_Differential);
#undef INSTANCE_KERNEL_Stress_Tensor_Differential
#else
}
#endif

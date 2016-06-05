//#####################################################################
//  Copyright (c) 2011-2013 Nathan Mitchell, Eftychios Sifakis.
//  This file is covered by the FreeBSD license. Please refer to the 
//  license.txt file for more information.
//#####################################################################


#ifndef SUBROUTINE_Rotated_Stress_Derivative_Contraction
#include <assert.h>
#include "KernelCommon.h"
#else
namespace {
#endif

template<class Tw,class T_DATA=void,class I_DATA=void>
#ifdef SUBROUTINE_Rotated_Stress_Derivative_Contraction
inline
#endif
void Rotated_Stress_Derivative_Contraction(const T_DATA (&dPdF)[12], const T_DATA (&dF_Hat)[9], T_DATA (&dP_Hat)[9])
{
    typedef enum { x11=0,x21,x31,x12,x22,x32,x13,x23,x33} Matrix_Entry;
    typedef enum { x1111=0,x1122,x1133,x2222,x2233,x3333,x1212,x1221,x1313,x1331,x2323,x2332 } RSD_Entry;


    typedef Number<Tw> Tn;

    Tn rdP_HatEntry;

    Tn rdF_HatEntry1;
    Tn rdF_HatEntry2;
    Tn rdF_HatEntry3;

    Tn rRSDEntry1;
    Tn rRSDEntry2;
    Tn rRSDEntry3;

    
    //dP_hat(1,1)=a1111*dF_hat(1,1)+a1122*dF_hat(2,2)+a1133*dF_hat(3,3);
    //dP_hat(2,2)=a1122*dF_hat(1,1)+a2222*dF_hat(2,2)+a2233*dF_hat(3,3);
    //dP_hat(3,3)=a1133*dF_hat(1,1)+a2233*dF_hat(2,2)+a3333*dF_hat(3,3);

    rdF_HatEntry1.Load(dF_Hat[x11]);
    rdF_HatEntry2.Load(dF_Hat[x22]);
    rdF_HatEntry3.Load(dF_Hat[x33]);

    rRSDEntry1.Load(dPdF[x1111]);
    rRSDEntry2.Load(dPdF[x1122]);
    rRSDEntry3.Load(dPdF[x1133]);

    rdP_HatEntry = rRSDEntry1*rdF_HatEntry1+rRSDEntry2*rdF_HatEntry2+rRSDEntry3*rdF_HatEntry3;
    Store(dP_Hat[x11],rdP_HatEntry);


    rRSDEntry1.Load(dPdF[x1122]);
    rRSDEntry2.Load(dPdF[x2222]);
    rRSDEntry3.Load(dPdF[x2233]);

    rdP_HatEntry = rRSDEntry1*rdF_HatEntry1+rRSDEntry2*rdF_HatEntry2+rRSDEntry3*rdF_HatEntry3;
    Store(dP_Hat[x22],rdP_HatEntry);


    rRSDEntry1.Load(dPdF[x1133]);
    rRSDEntry2.Load(dPdF[x2233]);
    rRSDEntry3.Load(dPdF[x3333]);

    rdP_HatEntry = rRSDEntry1*rdF_HatEntry1+rRSDEntry2*rdF_HatEntry2+rRSDEntry3*rdF_HatEntry3;
    Store(dP_Hat[x33],rdP_HatEntry);

    //dP_hat(1,2)=a1212*dF_hat(1,2)+a1221*dF_hat(2,1);
    //dP_hat(2,1)=a1221*dF_hat(1,2)+a1212*dF_hat(2,1);

    rdF_HatEntry1.Load(dF_Hat[x12]);
    rdF_HatEntry2.Load(dF_Hat[x21]);

    
    rRSDEntry1.Load(dPdF[x1212]);
    rRSDEntry2.Load(dPdF[x1221]);

    rdP_HatEntry = rRSDEntry1*rdF_HatEntry1+rRSDEntry2*rdF_HatEntry2;
    Store(dP_Hat[x12],rdP_HatEntry);

    rRSDEntry1.Load(dPdF[x1221]);
    rRSDEntry2.Load(dPdF[x1212]);

    rdP_HatEntry = rRSDEntry1*rdF_HatEntry1+rRSDEntry2*rdF_HatEntry2;
    Store(dP_Hat[x21],rdP_HatEntry);

    //dP_hat(1,3)=a1313*dF_hat(1,3)+a1331*dF_hat(3,1);
    //dP_hat(3,1)=a1331*dF_hat(1,3)+a1313*dF_hat(3,1);


    rdF_HatEntry1.Load(dF_Hat[x13]);
    rdF_HatEntry3.Load(dF_Hat[x31]);

    
    rRSDEntry1.Load(dPdF[x1313]);
    rRSDEntry3.Load(dPdF[x1331]);

    rdP_HatEntry = rRSDEntry1*rdF_HatEntry1+rRSDEntry3*rdF_HatEntry3;
    Store(dP_Hat[x13],rdP_HatEntry);

    rRSDEntry1.Load(dPdF[x1331]);
    rRSDEntry3.Load(dPdF[x1313]);

    rdP_HatEntry = rRSDEntry1*rdF_HatEntry1+rRSDEntry3*rdF_HatEntry3;
    Store(dP_Hat[x31],rdP_HatEntry);

    //dP_hat(2,3)=a2323*dF_hat(2,3)+a2332*dF_hat(3,2);
    //dP_hat(3,2)=a2332*dF_hat(2,3)+a2323*dF_hat(3,2);


    rdF_HatEntry2.Load(dF_Hat[x23]);
    rdF_HatEntry3.Load(dF_Hat[x32]);

    
    rRSDEntry2.Load(dPdF[x2323]);
    rRSDEntry3.Load(dPdF[x2332]);

    rdP_HatEntry = rRSDEntry2*rdF_HatEntry2+rRSDEntry3*rdF_HatEntry3;
    Store(dP_Hat[x23],rdP_HatEntry);

    rRSDEntry2.Load(dPdF[x2332]);
    rRSDEntry3.Load(dPdF[x2323]);

    rdP_HatEntry = rRSDEntry2*rdF_HatEntry2+rRSDEntry3*rdF_HatEntry3;
    Store(dP_Hat[x32],rdP_HatEntry);




}

#ifndef SUBROUTINE_Rotated_Stress_Derivative_Contraction
#define INSTANCE_KERNEL_Rotated_Stress_Derivative_Contraction(WIDTH) const WIDETYPE(float,WIDTH) (&dPdF)[12], const WIDETYPE(float,WIDTH) (&dF_Hat)[9], WIDETYPE(float,WIDTH) (&dP_Hat)[9]
INSTANCE_KERNEL(Rotated_Stress_Derivative_Contraction);
#undef INSTANCE_KERNEL_Rotated_Stress_Derivative_Contraction
#else
}
#endif

//#####################################################################
//  Copyright (c) 2011-2013 Nathan Mitchell, Eftychios Sifakis.
//  This file is covered by the FreeBSD license. Please refer to the 
//  license.txt file for more information.
//#####################################################################


#ifndef SUBROUTINE_Muscle_Differential
#include <assert.h>
#include "KernelCommon.h"
#else
namespace {
#endif


template<class Tw,class T_DATA=void,class I_DATA=void>
#ifdef SUBROUTINE_Muscle_Differential
inline
#endif
void Muscle_Differential(T_DATA (&dP_fiber)[9], const T_DATA (&dF)[9], const T_DATA (&fiber)[3],
                         const T_DATA (&Ffiber)[3], const T_DATA (&c1), const T_DATA (&c2) )
{
    
    #define M_I(x,y) (y * 3 + x )

    typedef Number<Tw> Tn;

    Tn fiber0;   Tn fiber1;   Tn fiber2; 
    Tn F_fiber0;   Tn F_fiber1;   Tn F_fiber2; 
    Tn dF0;   Tn dF1;   Tn dF2; 
    Tn dw0;   Tn dw1;   Tn dw2; 
    Tn q0;   Tn q1;   Tn q2; 
    Tn dw_dp; 
    Tn temp; 
    Tn C1;     Tn C2; 
    Tn dP_fiber0;  Tn dP_fiber1;  Tn dP_fiber2; 

    Vector3<Tn>        vdF_r;
    Vector3<Tn>        vfiber;
    Vector3<Tn>        vdw;
    Vector3<Tn>        vF_fiber;  
    Vector3<Tn>        vq;
    Vector3<Tn>        vdP_fiber;

    vfiber.Load_Aligned(fiber);
    
    // 
    //   TV dw = dF*fiber;
    //
    
    // Load first row of dF
    vdF_r.Load_Aligned(dF[M_I(0,0)],dF[M_I(0,1)],dF[M_I(0,2)]);    
    vdw.x = vdF_r.DotProduct(vfiber);

    // Load second row of dF
    vdF_r.Load_Aligned(dF[M_I(1,0)],dF[M_I(1,1)],dF[M_I(1,2)]);
    vdw.y = vdF_r.DotProduct(vfiber);

    // Load third row of dF
    vdF_r.Load_Aligned(dF[M_I(2,0)],dF[M_I(2,1)],dF[M_I(2,2)]);
    vdw.z = vdF_r.DotProduct(vfiber);
      
    //
    // TV q = c1 * dw + c2 * VECTOR<T,d>::Dot_Product(dw,F_fiber) * F_fiber;
    //
    
    vF_fiber.Load_Aligned(Ffiber);
    C1.Load(c1);
    C2.Load(c2);

    //--
    //-- TV q = C1 * dw + C2 * VECTOR<T,d>::Dot_Product(dw,F_fiber) * F_fiber;
    //--   

    vq = (vdw * C1) + (vF_fiber * (C2 * vdw.DotProduct(vF_fiber)));
        
    // 
    // MATRIX<T,d> dP_fiber=MATRIX<T,d>::Outer_Product(q,fiber);
    //
    
    vdP_fiber = vfiber*vq.x;
    vdP_fiber.Store( dP_fiber[M_I(0,0)], dP_fiber[M_I(0,1)], dP_fiber[M_I(0,2)] );

    vdP_fiber = vfiber*vq.y;
    vdP_fiber.Store( dP_fiber[M_I(1,0)], dP_fiber[M_I(1,1)], dP_fiber[M_I(1,2)] );

    vdP_fiber = vfiber*vq.z;
    vdP_fiber.Store( dP_fiber[M_I(2,0)], dP_fiber[M_I(2,1)], dP_fiber[M_I(2,2)] );

    
}


#ifdef SUBROUTINE_Muscle_Differential
}
#else
#define INSTANCE_KERNEL_Muscle_Differential(WIDTH) WIDETYPE(float,WIDTH) (&dP_fiber)[9], const WIDETYPE(float,WIDTH) (&dF)[9], const WIDETYPE(float,WIDTH) (&fiber)[3], const WIDETYPE(float,WIDTH) (&Ffiber)[3], const WIDETYPE(float,WIDTH) (&c1), const WIDETYPE(float,WIDTH) (&c2)
INSTANCE_KERNEL(Muscle_Differential);
#undef INSTANCE_KERNEL_Muscle_Differential
#endif

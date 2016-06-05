//#####################################################################
//  Copyright (c) 2011-2013 Nathan Mitchell, Eftychios Sifakis.
//  This file is covered by the FreeBSD license. Please refer to the 
//  license.txt file for more information.
//#####################################################################


#include <immintrin.h>
#include "KernelCommon.h"

#ifdef FORCE_INLINE
#include "../Unweighted_Accumulation/Unweighted_Accumulation.h"
#else
#define SUBROUTINE_Unweighted_Accumulation
#include "../Unweighted_Accumulation/Unweighted_Accumulation.cpp"
#undef  SUBROUTINE_Unweighted_Accumulation
#endif


template<class Tw,class T_DATA=void,class I_DATA=void>
void Muscle_Forces(T_DATA (&f)[3][8], const T_DATA (&fiber)[3],
                   const T_DATA (&Ffiber)[3], const T_DATA (&c1),
                   const T_DATA (&one_over_h), const T_DATA (&cell_volume))
{
    #ifndef M_I
    #define M_I(x,y) (y * 3 + x)
    #endif

    typedef Number<Tw> Tn;

    Tn fiber0;   Tn fiber1;   Tn fiber2; 
    Tn F_fiber0;   Tn F_fiber1;   Tn F_fiber2; 
    Tn C1;
    Tn P_fiber0;  Tn P_fiber1;  Tn P_fiber2; 
    
    Vector3<Tn> vfiber(fiber0, fiber1, fiber2);
    Vector3<Tn> vP_fiber(P_fiber0, P_fiber1, P_fiber2);

    vfiber.Load_Aligned(fiber);
    
    F_fiber0.Load(Ffiber[0]); F_fiber1.Load(Ffiber[1]); F_fiber2.Load(Ffiber[2]);
    C1.Load( c1 );
    
    T_DATA P_fiber[9];

    vP_fiber = vfiber * F_fiber0 * C1;
    vP_fiber.Store(P_fiber[M_I(0,0)], P_fiber[M_I(0,1)], P_fiber[M_I(0,2)]);

    vP_fiber = vfiber * F_fiber1 * C1;
    vP_fiber.Store(P_fiber[M_I(1,0)], P_fiber[M_I(1,1)], P_fiber[M_I(1,2)]);

    vP_fiber = vfiber * F_fiber2 * C1;
    vP_fiber.Store(P_fiber[M_I(2,0)], P_fiber[M_I(2,1)], P_fiber[M_I(2,2)]);

    T_DATA scale;
    Tn CellVolume;
    CellVolume.Load(cell_volume);
    CellVolume = Tn() - CellVolume;
    Store(scale, CellVolume);

    // Accumulation Part
    Unweighted_Accumulation<Tw,T_DATA,I_DATA>(f,
                                             P_fiber,
                                             one_over_h,
                                             scale);

}

#define INSTANCE_KERNEL_Muscle_Forces(WIDTH) \
    WIDETYPE(float,WIDTH) (&f)[3][8], \
    const WIDETYPE(float,WIDTH) (&fiber)[3], \
    const WIDETYPE(float,WIDTH) (&Ffiber)[3], \
    const WIDETYPE(float,WIDTH) (&c1), \
    const WIDETYPE(float,WIDTH) (&one_over_h), \
    const WIDETYPE(float,WIDTH) (&cell_volume) 
INSTANCE_KERNEL(Muscle_Forces);
#undef INSTANCE_KERNEL_Muscle_Forces


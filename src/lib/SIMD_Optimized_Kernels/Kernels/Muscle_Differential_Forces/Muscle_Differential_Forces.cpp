//#####################################################################
//  Copyright (c) 2011-2013 Nathan Mitchell, Eftychios Sifakis.
//  This file is covered by the FreeBSD license. Please refer to the 
//  license.txt file for more information.
//#####################################################################


#include <assert.h>

#include "KernelCommon.h"

#ifdef FORCE_INLINE
#include "../Unweighted_Gradient/Unweighted_Gradient.h"
#include "../Muscle_Differential/Muscle_Differential.h"
#include "../Unweighted_Accumulation/Unweighted_Accumulation.h"
#else
#define SUBROUTINE_Unweighted_Gradient
#include "../Unweighted_Gradient/Unweighted_Gradient.cpp"
#undef SUBROUTINE_Weighted_Gradient

#define SUBROUTINE_Muscle_Differential
#include "../Muscle_Differential/Muscle_Differential.cpp"
#undef SUBROUTINE_Muscle_Differential

#define SUBROUTINE_Weighted_Accumulation
#include "../Unweighted_Accumulation/Unweighted_Accumulation.cpp"
#undef SUBROUTINE_Weighted_Accumulation
#endif

template<class Tw,class T_DATA=void,class I_DATA=void>
    void Muscle_Differential_Forces(T_DATA (&df)[3][8], const T_DATA (&du)[3][8], const T_DATA (&fiber)[3],  const T_DATA (&Ffiber)[3], const T_DATA (&c1), const T_DATA (&c2), const T_DATA (&one_over_h), const T_DATA (&cell_volume))
{
    typedef Number<Tw> Tn;

    // Gradient Part
    BUILD_TDATA(dF,[9]);
    Unweighted_Gradient<Tw,T_DATA,I_DATA>(du,
                                         dF,
                                         one_over_h);

    // Differential Part
    BUILD_TDATA(dP_fiber,[9]);
    Muscle_Differential<Tw,T_DATA,I_DATA>(dP_fiber,
                                         dF,
                                         fiber,
                                         Ffiber,
                                         c1,
                                         c2);

    BUILD_TDATA(scale,);
    Tn CellVolume;
    CellVolume.Load(cell_volume);
    CellVolume = Tn() - CellVolume;
    Store(scale, CellVolume);

    // Accumulation Part
    Unweighted_Accumulation<Tw,T_DATA,I_DATA>(df,
                                             dP_fiber,
                                             one_over_h,
                                             scale);

  
}


#define INSTANCE_KERNEL_Muscle_Differential_Forces(WIDTH) WIDETYPE(float,WIDTH) (&df)[3][8], const WIDETYPE(float,WIDTH) (&du)[3][8], const WIDETYPE(float,WIDTH) (&fiber)[3],  const WIDETYPE(float,WIDTH) (&Ffiber)[3], const WIDETYPE(float,WIDTH) (&c1), const WIDETYPE(float,WIDTH) (&c2), const WIDETYPE(float,WIDTH) (&one_over_h), const WIDETYPE(float,WIDTH) (&cell_volume)
INSTANCE_KERNEL(Muscle_Differential_Forces);
#undef INSTANCE_KERNEL_Muscle_Differential_Forces


//#####################################################################
//  Copyright (c) 2011-2013 Nathan Mitchell, Eftychios Sifakis.
//  This file is covered by the FreeBSD license. Please refer to the 
//  license.txt file for more information.
//#####################################################################


#include <assert.h>

#include "KernelCommon.h"

#ifdef FORCE_INLINE
#include "../Muscle_Tension/Muscle_Tension.h"
#include "../Unweighted_Gradient/Unweighted_Gradient.h"
namespace {
    BUILD_CONSTANT(one,1.f);
}
#else

namespace {
    BUILD_CONSTANT(stretch_limit,3.f);
}

#define SUBROUTINE_Isotropic_Stress_Derivative
#include "../Muscle_Tension/Muscle_Tension.cpp"
#undef SUBROUTINE_Isotropic_Stress_Derivative

#define SUBROUTINE_Unweighted_Gradient
#include "../Unweighted_Gradient/Unweighted_Gradient.cpp"
#undef SUBROUTINE_Unweighted_Gradient

#endif

template<class Tw,class T_DATA=void,class I_DATA=void> 
void Muscle_Update_Position_Based_State(    const T_DATA (&u)[3][8], 
                                            const I_DATA (&muscle_id),
                                            const T_DATA (&fiber)[3],
                                            const T_DATA (&density),
                                            const T_DATA (&one_over_h),
                                            
                                            T_DATA (&c1),
                                            T_DATA (&c2),
                                            T_DATA (&F_fiber)[3],
                                            
                                            const float *activations,
                                            const float *fiber_max_stresses
                                            )
    {
        #define M_I(x,y) (y * 3 + x )
        typedef Number<Tw> Tn;
        typedef typename Number<Tw>::Mask Tm;

        BUILD_TDATA(F,[9]);
        
        Unweighted_Gradient<Tw,T_DATA,I_DATA>(u,F,one_over_h);
        
        Tn val;
        Tn rone;
        rone.Load_Aligned(one);
        for(int i=0;i<9;i+=4){
            val.Load(F[i]);
            val = val + rone;
            Store(F[i],val);
        }

        Tn fiber1, fiber2, fiber3;
        Tn Ffiber1, Ffiber2, Ffiber3;
        Tn F1, F2, F3;

        Vector3<Tn> rfiber;
        Vector3<Tn> rFfiber;
        Vector3<Tn> rF1, rF2, rF3;

        rfiber.Load_Aligned(fiber);    

        rF1.Load_Aligned(F[M_I(0, 0)], F[M_I(0, 1)], F[M_I(0, 2)] );
        rF2.Load_Aligned(F[M_I(1, 0)], F[M_I(1, 1)], F[M_I(1, 2)] );
        rF3.Load_Aligned(F[M_I(2, 0)], F[M_I(2, 1)], F[M_I(2, 2)] );

        rFfiber.x = rF1.DotProduct(rfiber);
        rFfiber.y = rF2.DotProduct(rfiber);
        rFfiber.z = rF3.DotProduct(rfiber);

        BUILD_TDATA(stretch, );

        Tn stretch_squared = rFfiber.Magnitude_Squared();
        Tn rstretch = stretch_squared.sqrt();
        Tn rstretch_limit;rstretch_limit.Load_Aligned(stretch_limit);
        rstretch=min(rstretch,rstretch_limit);
        Store(stretch, rstretch);

        BUILD_TDATA(tension, );
        BUILD_TDATA(tension_derivative, );

        BUILD_TDATA(activation, );
        BUILD_TDATA(fiber_max_stress, );
        Tn ractivation;
        Tn rfiber_max_stress;
        
        ractivation.Gather(activations, muscle_id);
        rfiber_max_stress.Gather(fiber_max_stresses, muscle_id);

        Store(activation, ractivation);
        Store(fiber_max_stress, rfiber_max_stress);

        Tension<Tw,T_DATA,I_DATA>(tension, stretch, activation, density, fiber_max_stress);
        Tension_Derivative<Tw,T_DATA,I_DATA>(tension_derivative, stretch, activation, density, fiber_max_stress);

        Tn rtension;
        Tn rtension_derivative;

        rtension.Load(tension);
        rtension_derivative.Load(tension_derivative);
        rtension_derivative = max( rtension_derivative, Tn() );

        Tn rc1 = rtension / rstretch;
        Tn rc2 = (rtension_derivative-rc1) / stretch_squared;
        
        Tn Zero;
        Tn rdensity;
        rdensity.Load(density);
        Tm IsZero = rdensity == Zero;

        rFfiber.x = blend(IsZero, rFfiber.x, Zero);
        rFfiber.y = blend(IsZero, rFfiber.y, Zero);
        rFfiber.z = blend(IsZero, rFfiber.z, Zero);
        rc1 = blend(IsZero, rc1, Zero);
        rc2 = blend(IsZero, rc2, Zero);

        rFfiber.Store(F_fiber);
        Store(c1, rtension);
        Store(c2, rtension_derivative);

};

#define INSTANCE_KERNEL_Muscle_Update_Position_Based_State(WIDTH) const WIDETYPE(float,WIDTH) (&u)[3][8], const WIDETYPE(int,WIDTH) (&muscle_id),  const WIDETYPE(float,WIDTH) (&fiber)[3],  const WIDETYPE(float,WIDTH) (&density),  const WIDETYPE(float,WIDTH) (&one_over_h),  WIDETYPE(float,WIDTH) (&c1),   WIDETYPE(float,WIDTH) (&c2),  WIDETYPE(float,WIDTH) (&F_fiber)[3], const float *activations,  const float *fiber_max_stresses
INSTANCE_KERNEL(Muscle_Update_Position_Based_State);
#undef INSTANCE_KERNEL_Muscle_Update_Position_Based_State




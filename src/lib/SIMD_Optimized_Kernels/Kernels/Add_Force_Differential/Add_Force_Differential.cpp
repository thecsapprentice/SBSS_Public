//#####################################################################
//  Copyright (c) 2011-2013 Nathan Mitchell, Eftychios Sifakis.
//  This file is covered by the FreeBSD license. Please refer to the 
//  license.txt file for more information.
//#####################################################################


#include <assert.h>
#include "KernelCommon.h"

#ifdef FORCE_INLINE
#include "../Unweighted_Gradient/Unweighted_Gradient.h"
#include "../Unweighted_Accumulation/Unweighted_Accumulation.h"
#include "../Matrix_Times_Transpose/Matrix_Times_Transpose.h"
#include "../Matrix_Transpose_Times/Matrix_Transpose_Times.h"
#include "../Stress_Tensor_Differential/Stress_Tensor_Differential.h"
#include "../Pressure_Force_Differential/Pressure_Force_Differential.h"
namespace {
    BUILD_CONSTANT(one, 1.f);
}
#else

#define SUBROUTINE_Unweighted_Gradient
#include "../Unweighted_Gradient/Unweighted_Gradient.cpp"
#undef SUBROUTINE_Unweighted_Gradient

#define SUBROUTINE_Unweighted_Accumulation
#include "../Unweighted_Accumulation/Unweighted_Accumulation.cpp"
#undef SUBROUTINE_Unweighted_Accumulation

#define SUBROUTINE_Matrix_Times_Transpose
#include "../Matrix_Times_Transpose/Matrix_Times_Transpose.cpp"
#undef SUBROUTINE_Matrix_Times_Transpose

#define SUBROUTINE_Matrix_Transpose_Times
#include "../Matrix_Transpose_Times/Matrix_Transpose_Times.cpp"
#undef SUBROUTINE_Matrix_Transpose_Times

#define SUBROUTINE_Stress_Tensor_Differential
#include "../Stress_Tensor_Differential/Stress_Tensor_Differential.cpp"
#undef SUBROUTINE_Stress_Tensor_Differential

#define SUBROUTINE_Pressure_Force_Differential
#include "../Pressure_Force_Differential/Pressure_Force_Differential.cpp"
#undef SUBROUTINE_Pressure_Force_Differential


#endif

template<class Tw,class T_DATA=void,class I_DATA=void> 
void Add_Force_Differential(const T_DATA (&du)[3][8], 
                              const T_DATA (&dp),
                              const T_DATA (&alpha_squared_over_kappa),
                              const T_DATA (&alpha),
                              const T_DATA (&one_over_h),
                              const T_DATA (&cell_volume),
                              const T_DATA (&Q_hat)[3],
                              const T_DATA (&U)[9],
                              const T_DATA (&V)[9],
                              const T_DATA (&dPdF)[12],
                              
                              T_DATA (&df)[3][8],
                              T_DATA (&dq))
{
        typedef Number<Tw> Tn;
        
        BUILD_TDATA(dF,[9]);
        BUILD_TDATA(dq_update,);

        BUILD_TDATA(neg_CellVolume,);
        Tn CellVolume;
        CellVolume.Load(cell_volume);

        
        Unweighted_Gradient<Tw,T_DATA,I_DATA>(du,dF,one_over_h);
        
        Matrix_Transpose_Times<Tw,T_DATA,I_DATA>(dF,U,dF);
        Matrix_Transpose_Times<Tw,T_DATA,I_DATA>(dF,V,dF);

        Pressure_Force_Differential<Tw,T_DATA,I_DATA>(dq_update, Q_hat, dF, dp, alpha, alpha_squared_over_kappa);
        Tn dq_old;
        Tn dq_upd;
        dq_old.Load(dq);
        dq_upd.Load(dq_update);
        dq_old = dq_old+(dq_upd*CellVolume);
        Store(dq, dq_old);


        Stress_Tensor_Differential<Tw,T_DATA,I_DATA>(dF, dPdF, dF, Q_hat, dp, alpha);

        Matrix_Times_Transpose<Tw,T_DATA,I_DATA>(V,dF,dF);
        Matrix_Times_Transpose<Tw,T_DATA,I_DATA>(U,dF,dF);
     
        CellVolume = Tn() - CellVolume;
        Store(neg_CellVolume,CellVolume);

        Unweighted_Accumulation<Tw,T_DATA,I_DATA>(df, dF, one_over_h, neg_CellVolume);


}

#define INSTANCE_KERNEL_Add_Force_Differential(WIDTH) \
                        const WIDETYPE(float,WIDTH) (&du)[3][8], \
                        const WIDETYPE(float,WIDTH) (&dp), \
                        const WIDETYPE(float,WIDTH) (&alpha_squared_over_kappa), \
                        const WIDETYPE(float,WIDTH) (&alpha), \
                        const WIDETYPE(float,WIDTH) (&one_over_h), \
                        const WIDETYPE(float,WIDTH) (&cell_volume), \
                        const WIDETYPE(float,WIDTH) (&Q_hat)[3], \
                        const WIDETYPE(float,WIDTH) (&U)[9], \
                        const WIDETYPE(float,WIDTH) (&V)[9], \
                        const WIDETYPE(float,WIDTH) (&dPdF)[12], \
                        WIDETYPE(float,WIDTH) (&df)[3][8], \
                        WIDETYPE(float,WIDTH) (&dq)
INSTANCE_KERNEL(Add_Force_Differential);
#undef INSTANCE_KERNEL_Add_Force_Differential

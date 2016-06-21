//#####################################################################
//  Copyright (c) 2011-2013 Nathan Mitchell, Eftychios Sifakis.
//  This file is covered by the FreeBSD license. Please refer to the 
//  license.txt file for more information.
//#####################################################################


#include <assert.h>

#include "KernelCommon.h"

#ifdef FORCE_INLINE
#include "../Isotropic_Stress_Derivative/Isotropic_Stress_Derivative.h"
#include "../Unweighted_Gradient/Unweighted_Gradient.h"
#include "../Singular_Value_Decomposition/Singular_Value_Decomposition.h"
#include "../Penalty_Measure_Gradient/Penalty_Measure_Gradient.h"
namespace {
    BUILD_CONSTANT(one,1.f);
}
#else

#define SUBROUTINE_Isotropic_Stress_Derivative
#include "../Isotropic_Stress_Derivative/Isotropic_Stress_Derivative.cpp"
#undef SUBROUTINE_Isotropic_Stress_Derivative

#define SUBROUTINE_Unweighted_Gradient
#include "../Unweighted_Gradient/Unweighted_Gradient.cpp"
#undef SUBROUTINE_Unweighted_Gradient

#define SUBROUTINE_Singular_Value_Decomposition
#include "../Singular_Value_Decomposition/Singular_Value_Decomposition.cpp"
#undef SUBROUTINE_Singular_Value_Decomposition

#define SUBROUTINE_Penalty_Measure_Gradient
#include "../Penalty_Measure_Gradient/Penalty_Measure_Gradient.cpp"
#undef SUBROUTINE_Penalty_Measure_Gradient

//#define SUBROUTINE_Compute_Diagonal_Contribution
//#include "../Compute_Diagonal_Contribution/Compute_Diagonal_Contribution.cpp"
//#undef SUBROUTINE_Compute_Diagonal_Contribution

#define SUBROUTINE_Compute_Cell_Matrix
#include "../Compute_Cell_Matrix/Compute_Cell_Matrix.cpp"
#undef SUBROUTINE_Compute_Cell_Matrix

#endif

template<class T_MATERIAL,class Tw,class T_DATA=void, class I_DATA=void> 
struct Update_Position_Based_State
{
    static void Run(const T_DATA (&u)[3][8], 
                    const T_DATA (&p),
                    const T_DATA (&mu),
                    const T_DATA (&mu_stab),
                    const T_DATA (&kappa),
                    const T_DATA (&alpha),
                    const T_DATA (&cutoff),
                    const T_DATA (&one_over_h),
                    const T_DATA (&cell_volume),
                    
                    T_DATA (&U)[9],
                    T_DATA (&V)[9],
                    T_DATA (&Sigma)[3],
                    T_DATA (&Q_hat)[3],
                    T_DATA (&dPdF)[12],
                    T_DATA (&d)[3][8])
    {
        typedef Number<Tw> Tn;

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

        Singular_Value_Decomposition<Tw,T_DATA,I_DATA>(F,U,Sigma,V);

        Tn co;
        co.Load(cutoff);
        for(int i=0;i<3;i++){
            val.Load(Sigma[i]);
            val = max(val,co);
            Store(Sigma[i],val);
        }
        
        Penalty_Measure_Gradient<T_MATERIAL,Tw,T_DATA,I_DATA>::Run(Sigma,Q_hat);

        Isotropic_Stress_Derivative<T_MATERIAL,Tw,T_DATA,I_DATA>::Run(dPdF,Sigma,p,mu,kappa,alpha,true);     
       
        //Compute_Diagonal_Contribution<Tw,T_DATA,I_DATA>(one_over_h, mu_stab, cell_volume, U, V, dPdF, d);
    }



    static void Run(const T_DATA (&u)[3][8], 
                    const T_DATA (&p),
                    const T_DATA (&mu),
                    const T_DATA (&mu_stab),
                    const T_DATA (&kappa),
                    const T_DATA (&alpha),
                    const T_DATA (&cutoff),
                    const T_DATA (&one_over_h),
                    const T_DATA (&cell_volume),
                    
                    T_DATA (&U)[9],
                    T_DATA (&V)[9],
                    T_DATA (&Sigma)[3],
                    T_DATA (&Q_hat)[3],
                    T_DATA (&dPdF)[12],
                    T_DATA (&d)[3][8],
                    T_DATA (&system_matrix)[300]
                    )
    {
        typedef Number<Tw> Tn;

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

        Singular_Value_Decomposition<Tw,T_DATA,I_DATA>(F,U,Sigma,V);

        Tn co;
        co.Load(cutoff);
        for(int i=0;i<3;i++){
            val.Load(Sigma[i]);
            val = max(val,co);
            Store(Sigma[i],val);
        }
        
        Penalty_Measure_Gradient<T_MATERIAL,Tw,T_DATA,I_DATA>::Run(Sigma,Q_hat);

        Isotropic_Stress_Derivative<T_MATERIAL,Tw,T_DATA,I_DATA>::Run(dPdF,Sigma,p,mu,kappa,alpha,true);     
       
        //Compute_Diagonal_Contribution<Tw,T_DATA,I_DATA>(one_over_h, mu_stab, cell_volume, U, V, dPdF, d);
        Compute_Cell_Matrix<Tw,T_DATA,I_DATA>(one_over_h, mu_stab, cell_volume, U, V, dPdF, system_matrix);
    }


};


#define INSTANCE_KERNEL_Update_Position_Based_State(WIDTH)   \
    const WIDETYPE(float,WIDTH) (&u)[3][8],                  \
    const WIDETYPE(float,WIDTH) (&p),                        \
    const WIDETYPE(float,WIDTH) (&mu),                       \
    const WIDETYPE(float,WIDTH) (&mu_stab),                  \
    const WIDETYPE(float,WIDTH) (&kappa),                    \
    const WIDETYPE(float,WIDTH) (&alpha),                    \
    const WIDETYPE(float,WIDTH) (&cutoff),                   \
    const WIDETYPE(float,WIDTH) (&one_over_h),               \
    const WIDETYPE(float,WIDTH) (&cell_volume),              \
    WIDETYPE(float,WIDTH) (&U)[9],                           \
    WIDETYPE(float,WIDTH) (&V)[9],                           \
    WIDETYPE(float,WIDTH) (&Sigma)[3],                       \
    WIDETYPE(float,WIDTH) (&Q_hat)[3],                       \
    WIDETYPE(float,WIDTH) (&dPdF)[12],                       \
    WIDETYPE(float,WIDTH) (&d)[3][8] 



INSTANCE_KERNEL_MATERIAL(Update_Position_Based_State,COROTATED_TAG);
INSTANCE_KERNEL_MATERIAL(Update_Position_Based_State,NEOHOOKEAN_TAG);
#if defined(BIPHASIC_SUPPORT)
INSTANCE_KERNEL_MATERIAL(Update_Position_Based_State,BIPHASIC_TAG);
#endif

#undef INSTANCE_KERNEL_Update_Position_Based_State

#define INSTANCE_KERNEL_Update_Position_Based_State(WIDTH)   \
    const WIDETYPE(float,WIDTH) (&u)[3][8],                  \
    const WIDETYPE(float,WIDTH) (&p),                        \
    const WIDETYPE(float,WIDTH) (&mu),                       \
    const WIDETYPE(float,WIDTH) (&mu_stab),                  \
    const WIDETYPE(float,WIDTH) (&kappa),                    \
    const WIDETYPE(float,WIDTH) (&alpha),                    \
    const WIDETYPE(float,WIDTH) (&cutoff),                   \
    const WIDETYPE(float,WIDTH) (&one_over_h),               \
    const WIDETYPE(float,WIDTH) (&cell_volume),              \
    WIDETYPE(float,WIDTH) (&U)[9],                           \
    WIDETYPE(float,WIDTH) (&V)[9],                           \
    WIDETYPE(float,WIDTH) (&Sigma)[3],                       \
    WIDETYPE(float,WIDTH) (&Q_hat)[3],                       \
    WIDETYPE(float,WIDTH) (&dPdF)[12],                       \
    WIDETYPE(float,WIDTH) (&d)[3][8],                        \
    WIDETYPE(float,WIDTH) (&system_matrix)[300]



INSTANCE_KERNEL_MATERIAL(Update_Position_Based_State,COROTATED_TAG);
INSTANCE_KERNEL_MATERIAL(Update_Position_Based_State,NEOHOOKEAN_TAG);
#if defined(BIPHASIC_SUPPORT)
INSTANCE_KERNEL_MATERIAL(Update_Position_Based_State,BIPHASIC_TAG);
#endif

#undef INSTANCE_KERNEL_Update_Position_Based_State

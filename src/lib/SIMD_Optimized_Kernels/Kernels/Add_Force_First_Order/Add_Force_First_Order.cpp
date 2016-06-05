//#####################################################################
//  Copyright (c) 2011-2013 Nathan Mitchell, Eftychios Sifakis.
//  This file is covered by the FreeBSD license. Please refer to the 
//  license.txt file for more information.
//#####################################################################


#include <assert.h>

#include "KernelCommon.h"


#ifdef FORCE_INLINE
#include "../Unweighted_Accumulation/Unweighted_Accumulation.h"
#include "../Matrix_Times_Transpose/Matrix_Times_Transpose.h"
#include "../Pressure_Force/Pressure_Force.h"
#include "../Piola_Kirchhoff_Stress_Tensor/Piola_Kirchhoff_Stress_Tensor.h"
namespace {
    BUILD_CONSTANT(one,1.f);
}
#else


#define SUBROUTINE_Unweighted_Accumulation
#include "../Unweighted_Accumulation/Unweighted_Accumulation.cpp"
#undef SUBROUTINE_Unweighted_Accumulation

#define SUBROUTINE_Matrix_Times_Transpose
#include "../Matrix_Times_Transpose/Matrix_Times_Transpose.cpp"
#undef SUBROUTINE_Matrix_Times_Transpose

#define SUBROUTINE_Pressure_Force
#include "../Pressure_Force/Pressure_Force.cpp"
#undef SUBROUTINE_Pressure_Force

#define SUBROUTINE_Piola_Kirchhoff_Stress_Tensor
#include "../Piola_Kirchhoff_Stress_Tensor/Piola_Kirchhoff_Stress_Tensor.cpp"
#undef SUBROUTINE_Piola_Kirchhoff_Stress_Tensor

#endif

template<class T_MATERIAL,class Tw,class T_DATA=void,class I_DATA=void> 
struct Add_Force_First_Order
{
    static void Run(const T_DATA (&u)[3][8],
                    const T_DATA (&p),

                    const T_DATA (&mu),
                    const T_DATA (&alpha),
                    const T_DATA (&alpha_sqr_over_kappa),
                    const T_DATA (&kappa),
                    const T_DATA (&one_over_h),
                    const T_DATA (&cell_volume),
                    const T_DATA (&U)[9],
                    const T_DATA (&V)[9],
                    const T_DATA (&Sigma)[3],
                    const T_DATA (&Q_Hat)[3],

                    T_DATA (&P_Hat)[3],
                    T_DATA (&f)[3][8],
                    T_DATA (&q))
    {
        typedef Number<Tw> Tn;
        typedef enum { x11=0,x21,x31,x12,x22,x32,x13,x23,x33} Matrix_Entry;
       
        Piola_Kirchhoff_Stress_Tensor<T_MATERIAL,Tw,T_DATA,I_DATA>::Run(P_Hat, Sigma, 
                                                                        Q_Hat, p, mu, alpha, kappa);
        
        T_DATA P[9];
        Tn rP;
        Store(P[x21], rP);
        Store(P[x31], rP);
        Store(P[x12], rP);
        Store(P[x32], rP);
        Store(P[x13], rP);
        Store(P[x23], rP);
        rP.Load(P_Hat[0]);
        Store(P[x11], rP);
        rP.Load(P_Hat[1]);
        Store(P[x22], rP);
        rP.Load(P_Hat[2]);
        Store(P[x33], rP);

        Matrix_Times_Transpose<Tw,T_DATA,I_DATA>(V,P,P);
        Matrix_Times_Transpose<Tw,T_DATA,I_DATA>(U,P,P);

        T_DATA neg_CellVolume;
        Tn CellVolume;
        CellVolume.Load(cell_volume);

        T_DATA qf;
        Pressure_Force<T_MATERIAL,Tw,T_DATA,I_DATA>::Run(qf, Sigma, p, alpha, alpha_sqr_over_kappa);
        Tn rq;
        Tn rqf;
        rq.Load(q);
        rqf.Load(qf);
        rq = rq + rqf * CellVolume;
        Store(q, rq);

        CellVolume = Tn() - CellVolume;
        Store(neg_CellVolume,CellVolume);
        Unweighted_Accumulation<Tw,T_DATA,I_DATA>(f, P, one_over_h, neg_CellVolume);
    }

};


#define INSTANCE_KERNEL_Add_Force_First_Order(WIDTH) \
                       const WIDETYPE(float,WIDTH) (&u)[3][8], \
                       const WIDETYPE(float,WIDTH) (&p), \
                       const WIDETYPE(float,WIDTH) (&mu), \
                       const WIDETYPE(float,WIDTH) (&alpha), \
                       const WIDETYPE(float,WIDTH) (&alpha_sqr_over_kappa), \
                       const WIDETYPE(float,WIDTH) (&kappa), \
                       const WIDETYPE(float,WIDTH) (&one_over_h), \
                       const WIDETYPE(float,WIDTH) (&cell_volume), \
                       const WIDETYPE(float,WIDTH) (&U)[9], \
                       const WIDETYPE(float,WIDTH) (&V)[9], \
                       const WIDETYPE(float,WIDTH) (&Sigma)[3], \
                       const WIDETYPE(float,WIDTH) (&Q_Hat)[3], \
                       WIDETYPE(float,WIDTH) (&P_Hat)[3], \
                       WIDETYPE(float,WIDTH) (&f)[3][8], \
                       WIDETYPE(float,WIDTH) (&q)
INSTANCE_KERNEL_MATERIAL(Add_Force_First_Order,COROTATED_TAG);
INSTANCE_KERNEL_MATERIAL(Add_Force_First_Order,NEOHOOKEAN_TAG);
INSTANCE_KERNEL_MATERIAL(Add_Force_First_Order,BIPHASIC_TAG);
#undef INSTANCE_KERNEL_Add_Force_First_Order

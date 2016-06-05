//#####################################################################
//  Copyright (c) 2011-2013 Nathan Mitchell, Eftychios Sifakis.
//  This file is covered by the FreeBSD license. Please refer to the 
//  license.txt file for more information.
//#####################################################################


/*

#ifndef SUBROUTINE_Masked_Math
#include <assert.h>
#include "KernelCommon.h"
#else
namespace {
#endif

template<class Tw, class TIw, class T_DATA, class TI_DATA, int width>
#ifndef SUBROUTINE_Masked_Math
inline
#endif
    void Masked_Add(T_DATA (&Result),
                    const T_DATA (&A), const T_DATA (&B),
                    const TI_DATA (&Mask))
{
    typedef Number<Tw> Tn;
    typedef Mask<TIw> Tm;

    Tn rA, rB, rResult;
    Tm rM;

    rA.Load(A);
    rB.Load(B);
    rM.Load(Mask);

    rResult = rA + rB;
    MaskStore( Result, rResult, rM );
}

template<class Tw, class TIw, class T_DATA, class TI_DATA, int width>
#ifndef SUBROUTINE_Masked_Math
inline
#endif
    void Masked_Subtract(T_DATA (&Result),
                         const T_DATA (&A), const T_DATA (&B),
                         const TI_DATA (&Mask))
{
    typedef Number<Tw> Tn;
    typedef Mask<TIw> Tm;


    Tn rA, rB, rResult;
    Tm rM;

    rA.Load(A);
    rB.Load(B);
    rM.Load(Mask);

    rResult = rA - rB;
    MaskStore( Result, rResult, rM );

}

template<class Tw, class TIw, class T_DATA, class TI_DATA, int width>
#ifndef SUBROUTINE_Masked_Math
inline
#endif
    void Masked_Times(T_DATA (&Result),
                      const T_DATA (&A), const T_DATA (&B),
                      const TI_DATA (&Mask))
{
    typedef Number<Tw> Tn;
    typedef Mask<TIw> Tm;


    Tn rA, rB, rResult;
    Tm rM;

    rA.Load(A);
    rB.Load(B);
    rM.Load(Mask);

    rResult = rA * rB;
    MaskStore( Result, rResult, rM );

}

template<class Tw, class TIw, class T_DATA, class TI_DATA, int width>
#ifndef SUBROUTINE_Masked_Math
    inline
#endif
void Masked_SAXPY_A(T_DATA (&Result),
                        const T_DATA (&A), const T_DATA (&constant),
                        const TI_DATA (&Mask))
{
    typedef Number<Tw> Tn;
    typedef Mask<TIw> Tm;


    Tn rA, rconstant, rResult;
    Tm rM;

    rA.Load(A);
    rconstant.Load(constant);
    rM.Load(Mask);

    rResult = rconstant * rA;
    MaskStore( Result, rResult, rM );

}

template<class Tw, class TIw, class T_DATA, class TI_DATA, int width>
#ifndef SUBROUTINE_Masked_Math
inline
#endif
    void Masked_SAXPY_AB(T_DATA (&Result),
                         const T_DATA (&A), const T_DATA (&B), const T_DATA (&constant), 
                         const TI_DATA (&Mask))
{
    typedef Number<Tw> Tn;
    typedef Mask<TIw> Tm;


    Tn rA, rB, rconstant, rResult;
    Tm rM;

    rA.Load(A);
    rB.Load(B);
    rconstant.Load(constant);
    rM.Load(Mask);

    rResult = rconstant * rA + rB;
    MaskStore( Result, rResult, rM );


}


#ifndef SUBROUTINE_Masked_Math

KERNEL_EXPLICIT_MASKEDDEF_SCALAR(Masked_Add, float (&Result), const float (&A), const float (&B), const int (&Mask));
KERNEL_EXPLICIT_MASKEDDEF_SIMD(Masked_Add, float (&Result)[8], const float (&A)[8], const float (&B)[8], const int (&Mask)[8]);

KERNEL_EXPLICIT_MASKEDDEF_SCALAR(Masked_Subtract, float (&Result), const float (&A), const float (&B), const int (&Mask));
KERNEL_EXPLICIT_MASKEDDEF_SIMD(Masked_Subtract, float (&Result)[8], const float (&A)[8], const float (&B)[8], const int (&Mask)[8]);

KERNEL_EXPLICIT_MASKEDDEF_SCALAR(Masked_Times, float (&Result), const float (&A), const float (&B), const int (&Mask));
KERNEL_EXPLICIT_MASKEDDEF_SIMD(Masked_Times, float (&Result)[8], const float (&A)[8], const float (&B)[8], const int (&Mask)[8]);

KERNEL_EXPLICIT_MASKEDDEF_SCALAR(Masked_SAXPY_A, float (&Result), const float (&A), const float (&constant), const int (&Mask));
KERNEL_EXPLICIT_MASKEDDEF_SIMD(Masked_SAXPY_A, float (&Result)[8], const float (&A)[8], const float (&constant)[8], const int (&Mask)[8]);

KERNEL_EXPLICIT_MASKEDDEF_SCALAR(Masked_SAXPY_AB, float (&Result), const float (&A), const float (&B), const float (&constant), const int (&Mask));
KERNEL_EXPLICIT_MASKEDDEF_SIMD(Masked_SAXPY_AB, float (&Result)[8], const float (&A)[8], const float (&B)[8], const float (&constant)[8], const int (&Mask)[8]);

#else
}
#endif


*/

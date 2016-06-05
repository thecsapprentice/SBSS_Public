//#####################################################################
//  Copyright (c) 2011-2013 Nathan Mitchell, Eftychios Sifakis.
//  This file is covered by the FreeBSD license. Please refer to the 
//  license.txt file for more information.
//#####################################################################


#include <assert.h>
#include "KernelCommon.h"



template<class Tw,class T_DATA=void,class I_DATA=void>
void NumberTest(const T_DATA (&A),
                const T_DATA (&B),
                const T_DATA (&C),
                const T_DATA (&D),

                T_DATA (&Add),
                T_DATA (&Multiply),
                T_DATA (&Subtract),
                T_DATA (&Divide),

                T_DATA (&LessThan),
                T_DATA (&GreaterThan),
                T_DATA (&LessEquals),
                T_DATA (&GreaterEquals),
                T_DATA (&Equals),

                T_DATA (&And),
                T_DATA (&Or),
                T_DATA (&Xor),
                T_DATA (&AndNot),
                T_DATA (&Not),

                T_DATA (&Sqrt),
                T_DATA (&RSqrt),
                T_DATA (&Log),
                T_DATA (&Exp),
                T_DATA (&Inverse),
                T_DATA (&AbsoluteValue),
                T_DATA (&Sign),

                T_DATA (&Minimum),
                T_DATA (&Maximum),
                T_DATA (&Blend) )
{
    typedef Number<Tw> Tn;
    typedef typename Number<Tw>::Mask Tm;

    Tn rA, rB, rC, rD;
    Tn Result;
    Tm MResult;


    rA.Load(A);    rB.Load(B);
    rC.Load(C);    rD.Load(D);

    
    // Tests

    // Arithmetic
    Result = rA + rB;
    Store(Add, Result);

    Result = rA * rB;
    Store(Multiply, Result);

    Result = rA - rB;
    Store(Subtract, Result);

    Result = rA / rB;
    Store(Divide, Result);

#if 0
    // Comparisons
    MResult = rA < rB;
    Store(LessThan, MResult);

    MResult = rA > rB;
    Store(GreaterThan, MResult);

    MResult = rA <= rB;
    Store(LessEquals, MResult);

    MResult = rA >= rB;
    Store(GreaterEquals, MResult);

    MResult = rA == rB;
    Store(Equals, MResult);
#endif
    
    // Bitwise

    Result = rA & rB;
    Store(And, Result);

    Result = rA | rB;
    Store(Or, Result);

    Result = rA ^ rB;
    Store(Xor, Result);

    Result = rA.andnot(rB);
    Store(AndNot, Result);

    Result = ~rA;
    Store(Not, Result);


    // Functions

    Result = rA.sqrt();
    Store(Sqrt, Result);

    Result = rA.rsqrt();
    Store(RSqrt, Result);
    
    Result = rA.log();
    Store(Log, Result);
    
    Result = rA.exp();
    Store(Exp, Result);

    Result = rA.inverse();
    Store(Inverse, Result);

    Result = rA.abs();
    Store(AbsoluteValue, Result);
    
    Result = rA.sign();
    Store(Sign, Result);
    
    Result = min(rA,rB);
    Store(Minimum, Result);
    
    Result = max(rA,rC);
    Store(Maximum, Result);

#if 1
    Result = blend( rC < rD, rA, rB);
    Store(Blend, Result);
#endif

}

#define INSTANCE_KERNEL_NumberTest(WIDTH)  const WIDETYPE(float,WIDTH) (&A),  const WIDETYPE(float,WIDTH) (&B), const WIDETYPE(float,WIDTH) (&C),  const WIDETYPE(float,WIDTH) (&D),  WIDETYPE(float,WIDTH) (&Add), WIDETYPE(float,WIDTH) (&Multiply), WIDETYPE(float,WIDTH) (&Subtract), WIDETYPE(float,WIDTH) (&Divide),  WIDETYPE(float,WIDTH) (&LessThan),  WIDETYPE(float,WIDTH) (&GreaterThan),   WIDETYPE(float,WIDTH) (&LessEquals),  WIDETYPE(float,WIDTH) (&GreaterEquals), WIDETYPE(float,WIDTH) (&Equals), WIDETYPE(float,WIDTH) (&And),  WIDETYPE(float,WIDTH) (&Or), WIDETYPE(float,WIDTH) (&Xor), WIDETYPE(float,WIDTH) (&AndNot), WIDETYPE(float,WIDTH) (&Not), WIDETYPE(float,WIDTH) (&Sqrt), WIDETYPE(float,WIDTH) (&RSqrt), WIDETYPE(float,WIDTH) (&Log),  WIDETYPE(float,WIDTH) (&Exp), WIDETYPE(float,WIDTH) (&Inverse),  WIDETYPE(float,WIDTH) (&AbsoluteValue), WIDETYPE(float,WIDTH) (&Sign),  WIDETYPE(float,WIDTH) (&Minimum),  WIDETYPE(float,WIDTH) (&Maximum), WIDETYPE(float,WIDTH) (&Blend) 
INSTANCE_KERNEL(NumberTest);
#undef INSTANCE_KERNEL_NumberTest

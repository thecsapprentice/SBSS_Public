//#####################################################################
//  Copyright (c) 2011-2013 Nathan Mitchell, Eftychios Sifakis.
//  This file is covered by the FreeBSD license. Please refer to the 
//  license.txt file for more information.
//#####################################################################


#ifndef SUBROUTINE_Matrix_Transpose_Times
#include <assert.h>
#include "KernelCommon.h"
#else
namespace {
#endif

template<class Tw,class T_DATA=void,class I_DATA=void>
#ifndef SUBROUTINE_Matrix_Transpose_Times
inline
#endif
void Matrix_Transpose_Times(const T_DATA (&A)[9], const T_DATA (&B)[9], T_DATA (&C)[9])
{
    typedef enum { x11=0,x21,x31,x12,x22,x32,x13,x23,x33} Matrix_Entry;

    typedef Number<Tw> Tn;

    Tn rB;

    Vector3<Tn> vA;
    Vector3<Tn> vC1;
    Vector3<Tn> vC2;
    Vector3<Tn> vC3;

    // Use 1st column of A & B
    
    vA.Load_Aligned(A[x11],A[x12],A[x13]);

    rB.Load(B[x11]);
    vC1=vA*rB;

    rB.Load(B[x12]);
    vC2=vA*rB;

    rB.Load(B[x13]);
    vC3=vA*rB;

    // Use 2nd column of A & B

    vA.Load_Aligned(A[x21],A[x22],A[x23]);

    rB.Load(B[x21]);
    vC1=vC1+vA*rB;

    rB.Load(B[x22]);
    vC2=vC2+vA*rB;

    rB.Load(B[x23]);
    vC3=vC3+vA*rB;

    // Use 3rd column of A & B

    vA.Load_Aligned(A[x31],A[x32],A[x33]);

    rB.Load(B[x31]);
    vC1=vC1+vA*rB;

    rB.Load(B[x32]);
    vC2=vC2+vA*rB;

    rB.Load(B[x33]);
    vC3=vC3+vA*rB;

    // Write result

    vC1.Store( C[x11], C[x21], C[x31] );
    vC2.Store( C[x12], C[x22], C[x32] );
    vC3.Store( C[x13], C[x23], C[x33] );
}

#ifndef SUBROUTINE_Matrix_Transpose_Times
#define INSTANCE_KERNEL_Matrix_Transpose_Times(WIDTH) const WIDETYPE(float,WIDTH) (&A)[9], const WIDETYPE(float,WIDTH) (&B)[9], WIDETYPE(float,WIDTH) (&C)[9]
INSTANCE_KERNEL(Matrix_Transpose_Times);
#undef INSTANCE_KERNEL_Matrix_Transpose_Times
#else
}
#endif

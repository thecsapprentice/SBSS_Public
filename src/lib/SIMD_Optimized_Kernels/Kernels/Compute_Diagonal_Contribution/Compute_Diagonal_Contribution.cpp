//#####################################################################
//  Copyright (c) 2011-2013 Nathan Mitchell, Eftychios Sifakis.
//  This file is covered by the FreeBSD license. Please refer to the
//  license.txt file for more information.
//#####################################################################

#ifndef SUBROUTINE_Compute_Diagonal_Contribution
#include <assert.h>
#include "KernelCommon.h"
#else
namespace {
#endif


#ifdef FORCE_INLINE
#include "../Unweighted_Accumulation/Unweighted_Accumulation.cpp"
#else
#define SUBROUTINE_Unweighted_Accumulation
#include "../Unweighted_Accumulation/Unweighted_Accumulation.cpp"
#undef SUBROUTINE_Unweighted_Accumulation
#endif

namespace CDC{
    BUILD_CONSTANT(seven,7.f);
    BUILD_CONSTANT(one,1.f);

    template<class Tw,class T_DATA=void, class I_DATA=void>
        void Transpose(T_DATA (&At)[9], const T_DATA (&A)[9])
    {
        typedef enum { x11=0,x21,x31,x12,x22,x32,x13,x23,x33} Matrix_Entry;
        typedef enum { t11=0,t21=3,t31=6,t12=1,t22=4,t32=7,t13=2,t23=5,t33=8} TransMatrix_Entry;
        typedef Number<Tw> Tn;
        
        Tn Entry;
        Entry.Load(A[x11]);
        Store(At[t11],Entry);

        Entry.Load(A[x21]);
        Store(At[t21],Entry);

        Entry.Load(A[x31]);
        Store(At[t31],Entry);

        Entry.Load(A[x12]);
        Store(At[t12],Entry);

        Entry.Load(A[x22]);
        Store(At[t22],Entry);

        Entry.Load(A[x32]);
        Store(At[t32],Entry);

        Entry.Load(A[x13]);
        Store(At[t13],Entry);

        Entry.Load(A[x23]);
        Store(At[t23],Entry);

        Entry.Load(A[x33]);
        Store(At[t33],Entry);
    }
    
    template<class Tw,class T_DATA=void, class I_DATA=void>
        void Build_M(T_DATA (&M)[9], const T_DATA (&dPdF)[12], const T_DATA (&W)[9], int row )
    {
        
        typedef enum { x1111=0,x1122,x1133,x2222,x2233,x3333,x1212,x1221,x1313,x1331,x2323,x2332 } RSD_Entry;
        typedef enum { x11=0,x21,x31,x12,x22,x32,x13,x23,x33} Matrix_Entry;
        typedef Number<Tw> Tn;
        
        Tn ra1111, ra1122, ra1133, ra2222, ra2233, ra3333, ra1212, ra1221, ra1313, ra1331, ra2323, ra2332;
        Tn rM[6];
        Tn rW[3];
        
        ra1111.Load(dPdF[x1111]);
        ra1122.Load(dPdF[x1122]);
        ra1133.Load(dPdF[x1133]);
        ra2222.Load(dPdF[x2222]);
        ra2233.Load(dPdF[x2233]);
        ra3333.Load(dPdF[x3333]);
        ra1212.Load(dPdF[x1212]);
        ra1221.Load(dPdF[x1221]);
        ra1313.Load(dPdF[x1313]);
        ra1331.Load(dPdF[x1331]);
        ra2323.Load(dPdF[x2323]);
        ra2332.Load(dPdF[x2332]);
        
        rW[0].Load(W[x11+row]);
        rW[1].Load(W[x12+row]);
        rW[2].Load(W[x13+row]);
        
        rM[0] = rW[0]*rW[0]*ra1111 + rW[1]*rW[1]*ra1212 + rW[2]*rW[2]*ra1313;
        rM[1] = rW[0]*rW[0]*ra1212 + rW[1]*rW[1]*ra2222 + rW[2]*rW[2]*ra2323; 
        rM[2] = rW[0]*rW[0]*ra1313 + rW[1]*rW[1]*ra2323 + rW[2]*rW[2]*ra3333; 
        rM[3] = rW[0]*rW[1]*ra1122 + rW[1]*rW[0]*ra1221;
        rM[4] = rW[0]*rW[2]*ra1133 + rW[2]*rW[0]*ra1331;
        rM[5] = rW[2]*rW[1]*ra2233 + rW[1]*rW[2]*ra2332;
        
        Store(M[x11], rM[0]);
        Store(M[x22], rM[1]);
        Store(M[x33], rM[2]);
        Store(M[x12], rM[3]);
        Store(M[x13], rM[4]);
        Store(M[x32], rM[5]);
        Store(M[x21], rM[3]);
        Store(M[x31], rM[4]);
        Store(M[x23], rM[5]);
    }
    
}

template<class Tw,class T_DATA=void, class I_DATA=void>
#ifdef SUBROUTINE_Compute_Diagonal_Contribution
inline
#endif
    void Compute_Diagonal_Contribution(const T_DATA (&one_over_h),
                                       const T_DATA (&mu_stab),
                                       const T_DATA (&cell_volume),
                                       const T_DATA (&U)[9],
                                       const T_DATA (&V)[9],
                                       const T_DATA (&dPdF)[12],
                                       T_DATA (&d)[3][8])
{
    typedef Number<Tw> Tn;
    typedef enum { x11=0,x21,x31,x12,x22,x32,x13,x23,x33} Matrix_Entry;

    Tn Zero;
    Tn rSev;
    Tn rmustab;
    Tn rcV;
    Tn rone;
    BUILD_TDATA(_one, );

    rSev.Load_Aligned( CDC::seven );
    rmustab.Load(mu_stab);
    rcV.Load(cell_volume);
    rone.Load_Aligned( CDC::one );
    Store(_one, rone);
    
    BUILD_TDATA(H,[3][8]);
    for( int k=0; k<3; k++)
        for( int l=0; l<8; l++)
            Store(H[k][l], Zero);
    
    BUILD_TDATA(Vt,[9]);
    CDC::Transpose<Tw,T_DATA,I_DATA>(Vt, V);
    Unweighted_Accumulation<Tw,T_DATA,I_DATA>(H,Vt,one_over_h,_one);

    Tn h[3];
    Tn result;
    Tn rM[3];
    for( int k=0; k<3; k++){
        CDC::Build_M<Tw,T_DATA,I_DATA>(Vt, dPdF, U, k);
        for( int l=0; l<8; l++){
            h[0].Load(H[0][l]);
            h[1].Load(H[1][l]);
            h[2].Load(H[2][l]);

            result = Tn();
            for( int p=0; p<3; p++){
                rM[0].Load(Vt[0+p]); rM[1].Load(Vt[3+p]); rM[2].Load(Vt[6+p]);
                result = result + h[p]*(rM[0]*h[0] + rM[1]*h[1] + rM[2]*h[2]);
            }
            result = Tn() - (result * rcV + rSev * rmustab);
            Store(d[k][l], result);                   
        }
    }

}

#ifdef SUBROUTINE_Compute_Diagonal_Contribution
}
#else
#define INSTANCE_KERNEL_Compute_Diagonal_Contribution(WIDTH) const WIDETYPE(float,WIDTH) (&one_over_h), const WIDETYPE(float,WIDTH) (&mu_stab), const WIDETYPE(float,WIDTH) (&cell_volume), const WIDETYPE(float,WIDTH) (&U)[9], const WIDETYPE(float,WIDTH) (&V)[9], const WIDETYPE(float,WIDTH) (&dPdF)[12], WIDETYPE(float,WIDTH) (&d)[3][8]
INSTANCE_KERNEL(Compute_Diagonal_Contribution);
#undef INSTANCE_KERNEL_Compute_Diagonal_Contribution
#endif

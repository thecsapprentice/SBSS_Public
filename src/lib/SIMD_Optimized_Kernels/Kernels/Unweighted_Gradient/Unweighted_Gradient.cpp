//#####################################################################
//  Copyright (c) 2011-2013 Nathan Mitchell, Eftychios Sifakis.
//  This file is covered by the FreeBSD license. Please refer to the 
//  license.txt file for more information.
//#####################################################################


#ifndef SUBROUTINE_Unweighted_Gradient
#include <assert.h>
#include "KernelCommon.h"
#else
namespace {
#endif

    namespace Unweighted_Gradient_DEFINES{
        BUILD_CONSTANT(K_ONE_OVER_FOUR,0.25f);
    }

template<class Tw,class T_DATA=void,class I_DATA=void>
#ifdef SUBROUTINE_Unweighted_Gradient
inline
#endif
    void Unweighted_Gradient(const T_DATA (&u)[3][8], T_DATA (&F)[9], const T_DATA (&one_over_h))
{
    #define M_I(x,y) (y * 3 + x )

    typedef Number<Tw> Tn;

    Tn Ru000;  Tn Ru001;  Tn Ru010;  Tn Ru011;  
    Tn Ru100;  Tn Ru101;  Tn Ru110;  Tn Ru111;  

    // Interpolated along a single axes (located on edges)
    Tn RuI00;  Tn RuI01;  Tn RuI10;  Tn RuI11;  
    Tn Ru0I0;  Tn Ru0I1;  Tn Ru1I0;  Tn Ru1I1;  

    // Differentiated along one axis, interpolated along one (located on faces)
    Tn RuID0;  Tn RuID1;  
    Tn RuI0D;  Tn RuI1D;  
    Tn RuDI0;  Tn RuDI1;  

    // Differentiated along one axis, interpolated along two (located at cell interior)
    Tn RuDII;  Tn RuIDI;  Tn RuIID;  
	
    // ONE_OVER_H
    Tn SCALE;  
    Tn ONE_OVER_FOUR;

    SCALE.Load(one_over_h);
    ONE_OVER_FOUR.Load_Aligned(Unweighted_Gradient_DEFINES::K_ONE_OVER_FOUR);
    SCALE = SCALE * ONE_OVER_FOUR;

    // V = 0
    for( int v =0; v < 3; v++)
        {
            Ru000.Load(u[v][0]);
            Ru001.Load(u[v][1]);
            Ru010.Load(u[v][2]);
            Ru011.Load(u[v][3]);
            Ru100.Load(u[v][4]);
            Ru101.Load(u[v][5]);
            Ru110.Load(u[v][6]);
            Ru111.Load(u[v][7]);

            RuI00 = Ru100 + Ru000;
            RuI01 = Ru101 + Ru001;
            RuI10 = Ru110 + Ru010;
            RuI11 = Ru111 + Ru011;
            
            Ru0I0 = Ru010 + Ru000;         
            Ru0I1 = Ru011 + Ru001;           
            Ru1I0 = Ru110 + Ru100;          
            Ru1I1 = Ru111 + Ru101;            
            
            RuID0 = RuI10 - RuI00;
            RuID1 = RuI11 - RuI01;
            
            RuI0D = RuI01 - RuI00;
            RuI1D = RuI11 - RuI10;
            
            RuDI0 = Ru1I0 - Ru0I0;
            RuDI1 = Ru1I1 - Ru0I1;
            
            RuDII = (RuDI1 + RuDI0) * SCALE;
            RuIDI = (RuID1 + RuID0) * SCALE;
            RuIID = (RuI1D + RuI0D) * SCALE;

            Store(F[M_I(v,0)], RuDII);
            Store(F[M_I(v,1)], RuIDI);
            Store(F[M_I(v,2)], RuIID);
        }
}

#ifdef SUBROUTINE_Unweighted_Gradient
}
#else
#define INSTANCE_KERNEL_Unweighted_Gradient(WIDTH) const WIDETYPE(float,WIDTH) (&u)[3][8], WIDETYPE(float,WIDTH) (&F)[9], const WIDETYPE(float,WIDTH) (&one_over_h)
INSTANCE_KERNEL(Unweighted_Gradient);
#undef INSTANCE_KERNEL_Unweighted_Gradient
#endif

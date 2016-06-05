//#####################################################################
//  Copyright (c) 2011-2013 Nathan Mitchell, Eftychios Sifakis.
//  This file is covered by the FreeBSD license. Please refer to the 
//  license.txt file for more information.
//#####################################################################


#ifndef SUBROUTINE_Unweighted_Accumulation
#include <assert.h>
#include "KernelCommon.h"
#else
namespace {
#endif

    namespace Unweighted_Accumulation_DEFINES {
        BUILD_CONSTANT(K_ONE_OVER_FOUR,0.25f);
    };

template<class Tw,class T_DATA=void,class I_DATA=void>
#ifdef SUBROUTINE_Unweighted_Accumulation
inline
#endif
void Unweighted_Accumulation(T_DATA (&u)[3][8], const T_DATA (&F)[9], const T_DATA (&one_over_h), const T_DATA (&scale))
{

#define M_I(x,y) (y * 3 + x )

    typedef Number<Tw> Tn;

    // Nodal values
    Tn Ru000;  Tn Ru001;  Tn Ru010;  Tn Ru011;  
    Tn Ru100;  Tn Ru101;  Tn Ru110;  Tn Ru111;  

    // Differentiated along one axis, interpolated along one (located on faces)
    Tn RuID0;  Tn RuID1;  
    Tn RuI0D;  Tn RuI1D;  
    Tn RuDI0;  Tn RuDI1;  

    // Differentiated along one axis, interpolated along two (located at cell interior)
    Tn RuDII;  Tn RuIDI;  Tn RuIID;  
	
    // ONE_OVER_H
    Tn ONE_OVER_H;  

    // SCALE_OVER_FOUR_H
    Tn ONE_OVER_FOUR;
    Tn SCALE_OVER_FOUR_H;  

    // Load Constants for computation
    ONE_OVER_H.Load(one_over_h);
    SCALE_OVER_FOUR_H.Load(scale);
    ONE_OVER_FOUR.Load_Aligned(Unweighted_Accumulation_DEFINES::K_ONE_OVER_FOUR);
    SCALE_OVER_FOUR_H = SCALE_OVER_FOUR_H * ONE_OVER_H * ONE_OVER_FOUR;

    for( int v = 0; v < 3; v++)
        {
            //----------------------------------------
            // Read F values (resulting from 2 interpolations and 1 differentiation)
            //----------------------------------------
            
            // Load row of F
            RuDII.Load(F[M_I(v,0)]);
            RuIDI.Load(F[M_I(v,1)]);
            RuIID.Load(F[M_I(v,2)]);

            RuDII = RuDII * SCALE_OVER_FOUR_H;
            RuIDI = RuIDI * SCALE_OVER_FOUR_H;
            RuIID = RuIID * SCALE_OVER_FOUR_H;
                       
            //----------------------------------------
            // Undo last interpolation
            // (results contain 1 interpolation and 1 differentiation)
            //----------------------------------------
            
            //RuDII=(1-bcd.w[2][cell_num])*RuDI0+bcd.w[2][cell_num]*RuDI1;
            RuDI1 = RuDII;
            RuDI0 = RuDI1;
            
            //RuIDI=(1-bcd.w[2][cell_num])*RuID0+bcd.w[2][cell_num]*RuID1;
            RuID1 = RuIDI;
            RuID0 = RuID1;
            
            //RuIID=(1-bcd.w[1][cell_num])*RuI0D+bcd.w[1][cell_num]*RuI1D;
            RuI1D = RuIID;
            RuI0D = RuI1D;
            
            //----------------------------------------
            // Undo differentiation
            // (results contain 1 interpolation)
            //----------------------------------------
            
            // Interpolated along a single axes (located on edges)
            Tn RuI00;  Tn RuI01;  Tn RuI10;  Tn RuI11;  
            Tn Ru0I0;  Tn Ru0I1;  Tn Ru1I0;  Tn Ru1I1;  

            //RuDI0=(Ru1I0-Ru0I0)*bcd.one_over_h;
            Ru1I0 = RuDI0;
            Ru0I0 = Ru0I0 - RuDI0;
            
            //RuDI1=(Ru1I1-Ru0I1)*bcd.one_over_h;
            Ru1I1 = RuDI1;
            Ru0I1 = Ru0I1 - RuDI1;
                  
            //RuI0D=(RuI01-RuI00)*bcd.one_over_h;
            RuI01 = RuI0D - RuID1;
            RuI00 = RuI00 -RuID0 - RuI0D;

            //RuI1D=(RuI11-RuI10)*bcd.one_over_h;
            RuI11 = RuID1 + RuI1D;
            RuI10 = RuID0 - RuI1D;
            
            //----------------------------------------
            // Load preexisting u values
            // (results need to be accumulated)
            //----------------------------------------
            
            Ru000.Load(u[v][0]);
            Ru001.Load(u[v][1]);
            Ru010.Load(u[v][2]);
            Ru011.Load(u[v][3]);
            Ru100.Load(u[v][4]);
            Ru101.Load(u[v][5]);
            Ru110.Load(u[v][6]);
            Ru111.Load(u[v][7]);
                      
            //----------------------------------------
            // Undo remaining interpolation
            // (results are nodal values)
            //----------------------------------------
            
            //RuI00=(1-bcd.w[0][cell_num])*Ru000+bcd.w[0][cell_num]*Ru100;
            Ru100 = Ru100 + RuI00;
            Ru000 = Ru000 + RuI00;
            
            //RuI01=(1-bcd.w[0][cell_num])*Ru001+bcd.w[0][cell_num]*Ru101;
            Ru101 = Ru101 + RuI01;
            Ru001 = Ru001 + RuI01;

            //RuI10=(1-bcd.w[0][cell_num])*Ru010+bcd.w[0][cell_num]*Ru110;
            Ru110 = Ru110 + RuI10;
            Ru010 = Ru010 + RuI10;
            
            //RuI11=(1-bcd.w[0][cell_num])*Ru011+bcd.w[0][cell_num]*Ru111;
            Ru111 = Ru111 + RuI11;
            Ru011 = Ru011 + RuI11;
            
            //Ru0I0=(1-bcd.w[1][cell_num])*Ru000+bcd.w[1][cell_num]*Ru010;
            Ru010 = Ru010 + Ru0I0;
            Ru000 = Ru000 + Ru0I0;

            //Ru0I1=(1-bcd.w[1][cell_num])*Ru001+bcd.w[1][cell_num]*Ru011;
            Ru011 = Ru011 + Ru0I1;
            Ru001 = Ru001 + Ru0I1;

            //Ru1I0=(1-bcd.w[1][cell_num])*Ru100+bcd.w[1][cell_num]*Ru110;
            Ru110 = Ru110 + Ru1I0;
            Ru100 = Ru100 + Ru1I0;
            
            //Ru1I1=(1-bcd.w[1][cell_num])*Ru101+bcd.w[1][cell_num]*Ru111; 
            Ru111 = Ru111 + Ru1I1;
            Ru101 = Ru101 + Ru1I1;
            
            // //----------------------------------------
            // // Apply Cell Weight
            // //----------------------------------------
            
            // VECTOR_MUL( Ru000, Ru000, CELL_WEIGHT );
            // VECTOR_MUL( Ru001, Ru001, CELL_WEIGHT );
            // VECTOR_MUL( Ru010, Ru010, CELL_WEIGHT );
            // VECTOR_MUL( Ru011, Ru011, CELL_WEIGHT );
            // VECTOR_MUL( Ru100, Ru100, CELL_WEIGHT );
            // VECTOR_MUL( Ru101, Ru101, CELL_WEIGHT );
            // VECTOR_MUL( Ru110, Ru110, CELL_WEIGHT );
            // VECTOR_MUL( Ru111, Ru111, CELL_WEIGHT );
            
            //----------------------------------------
            // Write accumulated values
            //----------------------------------------
            
            Store(u[v][0],Ru000);
            Store(u[v][1],Ru001);
            Store(u[v][2],Ru010);
            Store(u[v][3],Ru011);
            Store(u[v][4],Ru100);
            Store(u[v][5],Ru101);
            Store(u[v][6],Ru110);
            Store(u[v][7],Ru111);

        }
    
    
    
}


#ifndef SUBROUTINE_Unweighted_Accumulation
#define INSTANCE_KERNEL_Unweighted_Accumulation(WIDTH) WIDETYPE(float,WIDTH) (&u)[3][8], const WIDETYPE(float,WIDTH) (&F)[9], const WIDETYPE(float,WIDTH) (&one_over_h), const WIDETYPE(float,WIDTH) (&scale)
INSTANCE_KERNEL(Unweighted_Accumulation);
#undef INSTANCE_KERNEL_Unweighted_Accumulation
#else
}
#endif

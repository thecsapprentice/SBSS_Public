//#####################################################################
//  Copyright (c) 2011-2013 Nathan Mitchell, Eftychios Sifakis.
//  This file is covered by the FreeBSD license. Please refer to the 
//  license.txt file for more information.
//#####################################################################


#ifndef SUBROUTINE_Collision_Force_Differential
#include <assert.h>
#include "KernelCommon.h"
#else
namespace {
#endif

namespace {
}

template<class T_RAW,class T_DATA=void,class I_DATA=void>
#ifdef SUBROUTINE_Collision_Force_Differential
inline
#endif
void Collision_Force_Differential(T_DATA (&f)[3][8],
                      const T_DATA (&u)[3][8],
                      const T_DATA (&W)[3],
                      const I_DATA (&spring_id),
                      const float* extern_collision_stiffness )
{
    typedef Number<T_RAW> Tn;

    Tn Ru000;  Tn Ru001;  Tn Ru010;  Tn Ru011;  
    Tn Ru100;  Tn Ru101;  Tn Ru110;  Tn Ru111;  

    // Interpolated along a single axes (located on edges)
    Tn RuI00;  Tn RuI01;  Tn RuI10;  Tn RuI11;  

	// Interpolated along 2 axes
    Tn RuI0I;  Tn RuI1I;  

	// Interpolated along 3 axes
	Tn RuIII;  
	
    // Weight Values
    Tn w[3];   

    // External Location
    Tn RScale;

    Tn rf;

    w[0].Load(W[0]);
    w[1].Load(W[1]);
    w[2].Load(W[2]);
        
    RScale.Gather( extern_collision_stiffness, spring_id );

    for( int v = 0; v < 3; v++)
        {
            // Do the Interpolation Step

            Ru000.Load(u[v][0]);
            Ru001.Load(u[v][1]);
            Ru010.Load(u[v][2]);
            Ru011.Load(u[v][3]);
            Ru100.Load(u[v][4]);
            Ru101.Load(u[v][5]);
            Ru110.Load(u[v][6]);
            Ru111.Load(u[v][7]);

            RuI00 = (Ru100 - Ru000) * w[0] + Ru000;
            RuI01 = (Ru101 - Ru001) * w[0] + Ru001;
            RuI10 = (Ru110 - Ru010) * w[0] + Ru010;
            RuI11 = (Ru111 - Ru011) * w[0] + Ru011;
            
            RuI0I = (RuI01 - RuI00) * w[2] + RuI00;
            RuI1I = (RuI11 - RuI10) * w[2] + RuI10;

            RuIII = (RuI1I - RuI0I) * w[1] + RuI0I;     

            // Do the Scaling
            RuIII = RuIII * RScale;

            // Perform Extrapolation Step

            RuI1I = RuIII * w[1];
            RuI0I = RuIII - RuI1I;

            RuI11 = RuI1I * w[2];
            RuI10 = RuI1I - RuI11;           
            RuI01 = RuI0I * w[2];
            RuI00 = RuI0I - RuI01;

            Ru100 = RuI00 * w[0];
            Ru000 = RuI00 - Ru100;           
            Ru101 = RuI01 * w[0];
            Ru001 = RuI01 - Ru101;
            Ru110 = RuI10 * w[0];
            Ru010 = RuI10 - Ru110;           
            Ru111 = RuI11 * w[0];
            Ru011 = RuI11 - Ru111;

            // Accumulate onto f

            rf.Load(f[v][0]);
            rf = rf-Ru000;
            Store(f[v][0],rf);

            rf.Load(f[v][1]);
            rf = rf-Ru001;
            Store(f[v][1],rf);

            rf.Load(f[v][2]);
            rf = rf-Ru010;
            Store(f[v][2],rf);

            rf.Load(f[v][3]);
            rf = rf-Ru011;
            Store(f[v][3],rf);

            rf.Load(f[v][4]);
            rf = rf-Ru100;
            Store(f[v][4],rf);

            rf.Load(f[v][5]);
            rf = rf-Ru101;
            Store(f[v][5],rf);

            rf.Load(f[v][6]);
            rf = rf-Ru110;
            Store(f[v][6],rf);

            rf.Load(f[v][7]);
            rf = rf-Ru111;
            Store(f[v][7],rf);
        }

}


#ifdef SUBROUTINE_Collision_Force_Differential
}
#else
#define INSTANCE_KERNEL_Collision_Force_Differential(WIDTH) WIDETYPE(float,WIDTH) (&f)[3][8], const WIDETYPE(float,WIDTH) (&u)[3][8], const WIDETYPE(float,WIDTH) (&W)[3], const WIDETYPE(int,WIDTH) (&spring_id), const float* extern_collision_stiffness
INSTANCE_KERNEL(Collision_Force_Differential);
#undef INSTANCE_KERNEL_Collision_Force_Differential
#endif


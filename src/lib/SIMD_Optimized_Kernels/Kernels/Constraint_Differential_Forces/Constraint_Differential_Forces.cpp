//#####################################################################
//  Copyright (c) 2011-2013 Nathan Mitchell, Eftychios Sifakis.
//  This file is covered by the FreeBSD license. Please refer to the 
//  license.txt file for more information.
//#####################################################################


#ifndef SUBROUTINE_Constraint_Differential_Forces
#include <assert.h>
#include "KernelCommon.h"
#else
namespace {
#endif

template<class T_RAW,class T_DATA=void,class I_DATA=void>
#ifdef SUBROUTINE_Constraint_Differential_Forces
inline
#endif
void Constraint_Differential_Forces(T_DATA (&df)[3][8], const T_DATA (&u)[3][8],
                                    const T_DATA (&W)[3], const T_DATA (&scale) )
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
    Tn w0;  
    Tn w1;  
    Tn w2;  

    Tn RScale;
    Tn rdf;

    w0.Load(W[0]);
    w1.Load(W[1]);
    w2.Load(W[2]);
    RScale.Load(scale);

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

            RuI00 = (Ru100 - Ru000) * w0 + Ru000;
            RuI01 = (Ru101 - Ru001) * w0 + Ru001;
            RuI10 = (Ru110 - Ru010) * w0 + Ru010;
            RuI11 = (Ru111 - Ru011) * w0 + Ru011;
            
            RuI0I = (RuI01 - RuI00) * w2 + RuI00;
            RuI1I = (RuI11 - RuI10) * w2 + RuI10;

            RuIII = (RuI1I - RuI0I) * w1 + RuI0I;     
            
            // Do the Scaling

            RuIII = RuIII * RScale;

            // Perform Extrapolation Step

            RuI1I = RuIII * w1;
            RuI0I = RuIII - RuI1I;

            RuI11 = RuI1I * w2;
            RuI10 = RuI1I - RuI11;           
            RuI01 = RuI0I * w2;
            RuI00 = RuI0I - RuI01;

            Ru100 = RuI00 * w0;
            Ru000 = RuI00 - Ru100;           
            Ru101 = RuI01 * w0;
            Ru001 = RuI01 - Ru101;
            Ru110 = RuI10 * w0;
            Ru010 = RuI10 - Ru110;           
            Ru111 = RuI11 * w0;
            Ru011 = RuI11 - Ru111;

            // Accumulate onto df

            rdf.Load(df[v][0]);
            rdf = rdf-Ru000;
            Store(df[v][0],rdf);

            rdf.Load(df[v][1]);
            rdf = rdf-Ru001;
            Store(df[v][1],rdf);

            rdf.Load(df[v][2]);
            rdf = rdf-Ru010;
            Store(df[v][2],rdf);

            rdf.Load(df[v][3]);
            rdf = rdf-Ru011;
            Store(df[v][3],rdf);

            rdf.Load(df[v][4]);
            rdf = rdf-Ru100;
            Store(df[v][4],rdf);

            rdf.Load(df[v][5]);
            rdf = rdf-Ru101;
            Store(df[v][5],rdf);

            rdf.Load(df[v][6]);
            rdf = rdf-Ru110;
            Store(df[v][6],rdf);

            rdf.Load(df[v][7]);
            rdf = rdf-Ru111;
            Store(df[v][7],rdf);
        }

}


#ifdef SUBROUTINE_Constraint_Differential_Forces
}
#else
#define INSTANCE_KERNEL_Constraint_Differential_Forces(WIDTH) WIDETYPE(float,WIDTH) (&df)[3][8], const WIDETYPE(float,WIDTH) (&u)[3][8], const WIDETYPE(float,WIDTH) (&W)[3], const WIDETYPE(float,WIDTH) (&scale)
INSTANCE_KERNEL(Constraint_Differential_Forces);
#undef INSTANCE_KERNEL_Constraint_Differential_Forces
#endif

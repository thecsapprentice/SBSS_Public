//#####################################################################
//  Copyright (c) 2011-2013 Nathan Mitchell, Eftychios Sifakis.
//  This file is covered by the FreeBSD license. Please refer to the 
//  license.txt file for more information.
//#####################################################################


#ifndef SUBROUTINE_Expanded_Collision_Forces
#include <assert.h>
#include "KernelCommon.h"
#else
namespace {
#endif

namespace {
    BUILD_CONSTANT(Attach_Index_Scale,3);
    BUILD_CONSTANT(Attach_X_Shift,0);
    BUILD_CONSTANT(Attach_Y_Shift,1);
    BUILD_CONSTANT(Attach_Z_Shift,2);

    BUILD_ICONSTANT(IS_REFERENCE_SPRING_CONSTANT,1);
    BUILD_ICONSTANT(IS_REFERENCE_ATTACHMENT,2);
}

template<class T_RAW,class T_DATA=void,class I_DATA=void>
#ifdef SUBROUTINE_Expanded_Collision_Forces
inline
#endif
void Expanded_Collision_Forces(T_DATA (&f)[3][8],
                      const T_DATA (&u)[3][8],
                      const T_DATA (&N)[3],
                      const T_DATA (&W)[3],
                      const T_DATA (&h),
                      const I_DATA (&flag),
                      const T_DATA (&stiffness_constant),
                      const T_DATA (&attachment_point)[3],
                      const I_DATA (&spring_id),
                      const I_DATA (&spring_id_X),
                      const I_DATA (&spring_id_Y),
                      const I_DATA (&spring_id_Z),
                      const float* extern_collision_attach,
                      const float* extern_collision_stiffness )
{
    typedef Number<T_RAW> Tn;
    typedef typename NumberPolicy< Number<T_RAW> >::DISCRETE_TYPE Dn;
    typedef typename NumberPolicy< Number<T_RAW> >::MASK_TYPE Tm;

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

    // Node Values
    Tn n[3];  

    // Type Flag
    Dn rFlag;

    // Attachment Point
    Tn rAttach[3];  
    Tn rScale;

    // External Location
    Tn rExtAttach[3];
    Tn rExtScale;

    Tn rDx;
    Tn rf;

    w[0].Load(W[0]);
    w[1].Load(W[1]);
    w[2].Load(W[2]);

    n[0].Load(N[0]);
    n[1].Load(N[1]);
    n[2].Load(N[2]);
        


    rDx.Load(h);
  
    // On a per channel basis, decide if we have hardcoded attachment/spring constants, or not
    rFlag.Load(flag);
    Dn Zero;
    Dn TEST_Is_Reference_Spring_Constant;
    Dn TEST_Is_Reference_Attachment;
    TEST_Is_Reference_Spring_Constant.Load( IS_REFERENCE_SPRING_CONSTANT );
    TEST_Is_Reference_Attachment.Load( IS_REFERENCE_ATTACHMENT );

    Tm Is_Reference_Spring_Constant;
    Tm Is_Reference_Attachment;

    Is_Reference_Spring_Constant = (rFlag & TEST_Is_Reference_Spring_Constant) > Zero;   
    Is_Reference_Attachment = (rFlag & TEST_Is_Reference_Attachment) > Zero;

    rScale.Load(stiffness_constant);
    rExtScale.Gather( extern_collision_stiffness, spring_id );
    rAttach[0].Load(attachment_point[0]);
    rAttach[1].Load(attachment_point[1]);
    rAttach[2].Load(attachment_point[2]);
    rExtAttach[0].Gather( extern_collision_attach, spring_id_X );
    rExtAttach[1].Gather( extern_collision_attach, spring_id_Y );
    rExtAttach[2].Gather( extern_collision_attach, spring_id_Z );

    rScale = blend( Is_Reference_Spring_Constant, rScale, rExtScale );
    rAttach[0] = blend( Is_Reference_Attachment, rAttach[0], rExtAttach[0] );
    rAttach[1] = blend( Is_Reference_Attachment, rAttach[1], rExtAttach[1] );
    rAttach[2] = blend( Is_Reference_Attachment, rAttach[2], rExtAttach[2] );

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

            // Compute a Real Deformation
            RuIII = RuIII + n[v] + rDx * w[v];

            // Compute Real Offset
            RuIII = rAttach[v] - RuIII;
            
            // Do the Scaling
            RuIII = RuIII * rScale;

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
            rf = rf+Ru000;
            Store(f[v][0],rf);

            rf.Load(f[v][1]);
            rf = rf+Ru001;
            Store(f[v][1],rf);

            rf.Load(f[v][2]);
            rf = rf+Ru010;
            Store(f[v][2],rf);

            rf.Load(f[v][3]);
            rf = rf+Ru011;
            Store(f[v][3],rf);

            rf.Load(f[v][4]);
            rf = rf+Ru100;
            Store(f[v][4],rf);

            rf.Load(f[v][5]);
            rf = rf+Ru101;
            Store(f[v][5],rf);

            rf.Load(f[v][6]);
            rf = rf+Ru110;
            Store(f[v][6],rf);

            rf.Load(f[v][7]);
            rf = rf+Ru111;
            Store(f[v][7],rf);
        }

}


#ifdef SUBROUTINE_Expanded_Collision_Forces
}
#else
#define INSTANCE_KERNEL_Expanded_Collision_Forces(WIDTH) WIDETYPE(float,WIDTH) (&f)[3][8], \
                                                         const WIDETYPE(float,WIDTH) (&u)[3][8], \
                                                         const WIDETYPE(float,WIDTH) (&N)[3], \ 
                                                         const WIDETYPE(float,WIDTH) (&W)[3], \
                                                         const WIDETYPE(float,WIDTH) (&h), \
                                                         const WIDETYPE(int,WIDTH) (&flag), \
                                                         const WIDETYPE(float,WIDTH) (&stiffness_constant), \
                                                         const WIDETYPE(float,WIDTH) (&attachment_point)[3], \
                                                         const WIDETYPE(int,WIDTH) (&spring_id), \
                                                         const WIDETYPE(int,WIDTH) (&spring_id_X), \
                                                         const WIDETYPE(int,WIDTH) (&spring_id_Y), \ 
                                                         const WIDETYPE(int,WIDTH) (&spring_id_Z), \
                                                         const float* extern_collision_attach, \
                                                         const float* extern_collision_stiffness
INSTANCE_KERNEL(Expanded_Collision_Forces);
#undef INSTANCE_KERNEL_Expanded_Collision_Forces
#endif


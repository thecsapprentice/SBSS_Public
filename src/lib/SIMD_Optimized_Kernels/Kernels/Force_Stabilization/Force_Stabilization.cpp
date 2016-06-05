//#####################################################################
//  Copyright (c) 2011-2013 Nathan Mitchell, Eftychios Sifakis.
//  This file is covered by the FreeBSD license. Please refer to the 
//  license.txt file for more information.
//#####################################################################


#ifndef SUBROUTINE_Force_Stabilization
#include <assert.h>
#include "KernelCommon.h"
#else
namespace {
#endif

    namespace nm_Force_Stabilization{
        BUILD_CONSTANT(minusseven,-7.f);
        BUILD_CONSTANT(three,3.f);
        BUILD_CONSTANT(minusfive,-5.f);
        BUILD_CONSTANT(one,1.f);
    }
    
    using namespace nm_Force_Stabilization;

    
    template<class Tw,class T_DATA=void,class I_DATA=void>
#ifdef SUBROUTINE_Force_Stabilization
    inline
#endif
    void Force_Stabilization(const T_DATA (&Du)[3][8], const T_DATA (&constant), T_DATA (&dH)[3][8] )
    {
        typedef Number<Tw> Tn;

        const int OpOrder[8][8] = {
            {0, 1, 2, 4, 3, 5, 6, 7},
            {1, 0, 3, 5, 2, 4, 7, 6},
            {2, 3, 0, 6, 1, 7, 4, 5},
            {3, 2, 1, 7, 0, 6, 5, 4},
            {4, 5, 6, 0, 7, 1, 2, 3},
            {5, 4, 7, 1, 6, 0, 3, 2},
            {6, 7, 4, 2, 5, 3, 0, 1},
            {7, 6, 5, 3, 4, 2, 1, 0}};

        Tn Ru[8];
        Tn Rf[8];
      
        Tn TEMP;
        Tn Rminusseven;
        Tn Rthree;
        Tn Rone;
        Tn Rminusfive;
        Tn RConstant;
        Tn RWeight;
        
        Rminusseven.Load_Aligned(minusseven);
        Rthree.Load_Aligned(three);
        Rone.Load_Aligned(one);
        Rminusfive.Load_Aligned(minusfive);
        
        RConstant.Load( constant );
        /*RWeight.Load( weight );*/
        RConstant = RConstant;// * RWeight;

        for( int d = 0; d < 3; d++ )
            {
                for( int i = 0; i < 8; i++)
                    Ru[i].Load(Du[d][i]);
                
#if 0
                for( int i = 0; i < 8; i++)
                    {
                        Rf[i] = 
                            (Rminusseven * Ru[OpOrder[i][0]]) + 
                            (Rthree * Ru[OpOrder[i][1]])+
                            (Rthree * Ru[OpOrder[i][2]])+
                            (Rthree * Ru[OpOrder[i][3]])+
                            (Rone * Ru[OpOrder[i][4]])+
                            (Rone * Ru[OpOrder[i][5]])+
                            (Rone * Ru[OpOrder[i][6]])+
                            (Rminusfive * Ru[OpOrder[i][7]]);
                    }
#else
                Rf[0] = 
                    (Rminusseven * Ru[0]) + 
                    (Rthree * Ru[1])+
                    (Rthree * Ru[2])+
                    (Rthree * Ru[4])+
                    (Rone * Ru[3])+
                    (Rone * Ru[5])+
                    (Rone * Ru[6])+
                    (Rminusfive * Ru[7]);

                Rf[1] = 
                    (Rminusseven * Ru[1]) + 
                    (Rthree * Ru[0])+
                    (Rthree * Ru[3])+
                    (Rthree * Ru[5])+
                    (Rone * Ru[2])+
                    (Rone * Ru[4])+
                    (Rone * Ru[7])+
                    (Rminusfive * Ru[6]);

                Rf[2] = 
                    (Rminusseven * Ru[2]) + 
                    (Rthree * Ru[3])+
                    (Rthree * Ru[0])+
                    (Rthree * Ru[6])+
                    (Rone * Ru[1])+
                    (Rone * Ru[7])+
                    (Rone * Ru[4])+
                    (Rminusfive * Ru[5]);

                Rf[3] = 
                    (Rminusseven * Ru[3]) + 
                    (Rthree * Ru[2])+
                    (Rthree * Ru[1])+
                    (Rthree * Ru[7])+
                    (Rone * Ru[0])+
                    (Rone * Ru[6])+
                    (Rone * Ru[5])+
                    (Rminusfive * Ru[4]);

                Rf[4] = 
                    (Rminusseven * Ru[4]) + 
                    (Rthree * Ru[5])+
                    (Rthree * Ru[6])+
                    (Rthree * Ru[0])+
                    (Rone * Ru[7])+
                    (Rone * Ru[1])+
                    (Rone * Ru[2])+
                    (Rminusfive * Ru[3]);

                Rf[5] = 
                    (Rminusseven * Ru[5]) + 
                    (Rthree * Ru[4])+
                    (Rthree * Ru[7])+
                    (Rthree * Ru[1])+
                    (Rone * Ru[6])+
                    (Rone * Ru[0])+
                    (Rone * Ru[3])+
                    (Rminusfive * Ru[2]);

                Rf[6] = 
                    (Rminusseven * Ru[6]) + 
                    (Rthree * Ru[7])+
                    (Rthree * Ru[4])+
                    (Rthree * Ru[2])+
                    (Rone * Ru[5])+
                    (Rone * Ru[3])+
                    (Rone * Ru[0])+
                    (Rminusfive * Ru[1]);

                Rf[7] = 
                    (Rminusseven * Ru[7]) + 
                    (Rthree * Ru[6])+
                    (Rthree * Ru[5])+
                    (Rthree * Ru[3])+
                    (Rone * Ru[4])+
                    (Rone * Ru[2])+
                    (Rone * Ru[1])+
                    (Rminusfive * Ru[0]);
#endif

            

                for( int i = 0; i < 8; i++)
                    {
                        Ru[i].Load(dH[d][i]);
                        Rf[i] = Rf[i] * RConstant;
                        Rf[i] = Rf[i] + Ru[i];
                        Store( dH[d][i], Rf[i]);
                    }
            }

    }

#ifdef SUBROUTINE_Force_Stabilization
}
#else
#define INSTANCE_KERNEL_Force_Stabilization(WIDTH) const WIDETYPE(float,WIDTH) (&Du)[3][8], const WIDETYPE(float,WIDTH) (&constant),  WIDETYPE(float,WIDTH) (&dH)[3][8]
INSTANCE_KERNEL(Force_Stabilization);
#undef INSTANCE_KERNEL_Force_Stabilization
#endif

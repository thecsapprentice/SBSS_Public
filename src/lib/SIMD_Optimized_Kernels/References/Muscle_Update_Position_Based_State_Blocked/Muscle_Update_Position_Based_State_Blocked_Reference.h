//#####################################################################
//  Copyright (c) 2011-2013 Nathan Mitchell, Eftychios Sifakis.
//  This file is covered by the FreeBSD license. Please refer to the 
//  license.txt file for more information.
//#####################################################################


template<class T, class I> 
struct Muscle_Update_Position_Based_State_Blocked_Reference
{
static void Run(const float u[3][27], 
                const I muscle_id[8],
                const T fiber[3][8],
                const T density[8],
                const T one_over_h[8],
                
                T c1[8],
                T c2[8],
                T F_fiber[3][8],
                
                const float *activations,
                const float *fiber_max_stresses);
};


template<class T, class I>
    bool Muscle_Update_Position_Based_State_Blocked_Compare(const T c1[8],
                                                            const T c2[8],
                                                            const T F_fiber[3][8],
                                                            const T c1_reference[8],
                                                            const T c2_reference[8],
                                                            const T F_fiber_reference[3][8]
                                                            );
                                         


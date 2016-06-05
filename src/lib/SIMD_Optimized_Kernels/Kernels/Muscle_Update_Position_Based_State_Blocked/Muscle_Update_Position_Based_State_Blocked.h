//#####################################################################
//  Copyright (c) 2011-2013 Nathan Mitchell, Eftychios Sifakis.
//  This file is covered by the FreeBSD license. Please refer to the 
//  license.txt file for more information.
//#####################################################################


template<class Tw,class T_DATA=void,class I_DATA=void> 
void Muscle_Update_Position_Based_State_Blocked(const float (&u)[3][27], 
                                                const I_DATA (&muscle_id),
                                                const T_DATA (&fiber)[3],
                                                const T_DATA (&density),
                                                const T_DATA (&one_over_h),
                                                
                                                T_DATA (&c1),
                                                T_DATA (&c2),
                                                T_DATA (&F_fiber)[3],
                                                
                                                const float *activations,
                                                const float *fiber_max_stresses
                                                );


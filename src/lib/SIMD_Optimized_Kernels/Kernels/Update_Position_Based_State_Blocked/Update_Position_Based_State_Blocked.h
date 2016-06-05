//#####################################################################
//  Copyright (c) 2011-2013 Nathan Mitchell, Eftychios Sifakis.
//  This file is covered by the FreeBSD license. Please refer to the 
//  license.txt file for more information.
//#####################################################################


template<class T_MATERIAL,class Tw,class T_DATA=void,class I_DATA=void> 
struct Update_Position_Based_State_Blocked
{
    static void Run(const float (&u)[3][27], 
                    const T_DATA (&p),
                    const T_DATA (&mu),
                    const T_DATA (&kappa),
                    const T_DATA (&alpha),
                    const T_DATA (&cutoff),
                    const T_DATA (&one_over_h),
                    
                    T_DATA (&U)[9],
                    T_DATA (&V)[9],
                    T_DATA (&Sigma)[3],
                    T_DATA (&Q_hat)[3],
                    T_DATA (&dPdF)[12]);
};

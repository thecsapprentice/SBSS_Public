//#####################################################################
//  Copyright (c) 2011-2013 Nathan Mitchell, Eftychios Sifakis.
//  This file is covered by the FreeBSD license. Please refer to the 
//  license.txt file for more information.
//#####################################################################


template<class T_RAW,class T_DATA=void,class I_DATA=void>
    void Stress_Tensor_Differential(T_DATA (&dP_hat)[9], const T_DATA (&dPdF)[12], const T_DATA (&dF_hat)[9],
                                    const T_DATA (&Q_hat)[3], const T_DATA (&dp), const T_DATA (&alpha));

//#####################################################################
//  Copyright (c) 2011-2013 Nathan Mitchell, Eftychios Sifakis.
//  This file is covered by the FreeBSD license. Please refer to the 
//  license.txt file for more information.
//#####################################################################


template<class Tw,class T_DATA=void,class I_DATA=void> 
void Add_Force_Differential(const T_DATA (&du)[3][8], 
                            const T_DATA (&dp),
                            const T_DATA (&alpha_squared_over_kappa),
                            const T_DATA (&alpha),
                            const T_DATA (&one_over_h),
                            const T_DATA (&cell_volume),
                            const T_DATA (&Q_hat)[3],
                            const T_DATA (&U)[9],
                            const T_DATA (&V)[9],
                            const T_DATA (&dPdF)[12],
                            
                            T_DATA (&df)[3][8],
                            T_DATA (&dq));


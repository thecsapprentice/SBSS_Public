//#####################################################################
//  Copyright (c) 2011-2013 Nathan Mitchell, Eftychios Sifakis.
//  This file is covered by the FreeBSD license. Please refer to the 
//  license.txt file for more information.
//#####################################################################


template<class Tw,class T_DATA=void, class I_DATA=void> 
    void Compute_Cell_Matrix(const T_DATA (&one_over_h),
                             const T_DATA (&mu_stab),
                             const T_DATA (&cell_volume),                                       
                             const T_DATA (&U)[9],
                             const T_DATA (&V)[9],
                             const T_DATA (&dPdF)[12],
                             T_DATA (&matrix)[300]);

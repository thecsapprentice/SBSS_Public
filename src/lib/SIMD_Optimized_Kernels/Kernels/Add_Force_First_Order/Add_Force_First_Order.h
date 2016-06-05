//#####################################################################
//  Copyright (c) 2011-2013 Nathan Mitchell, Eftychios Sifakis.
//  This file is covered by the FreeBSD license. Please refer to the 
//  license.txt file for more information.
//#####################################################################


template<class T_MATERIAL,class Tw,class T_DATA=void, class I_DATA=void>
struct Add_Force_First_Order
{
    static void Run(const T_DATA (&u)[3][8],
                    const T_DATA (&p),

                    const T_DATA (&mu),
                    const T_DATA (&alpha),
                    const T_DATA (&alpha_sqr_over_kappa),
                    const T_DATA (&kappa), 
                    const T_DATA (&one_over_h),
                    const T_DATA (&cell_volume),
                    const T_DATA (&U)[9],
                    const T_DATA (&V)[9],
                    const T_DATA (&Sigma)[3],
                    const T_DATA (&Q_Hat)[3],

                    T_DATA (&P_Hat)[3],
                    T_DATA (&f)[3][8],
                    T_DATA (&q));

};



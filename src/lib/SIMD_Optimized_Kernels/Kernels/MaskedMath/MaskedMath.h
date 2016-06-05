//#####################################################################
//  Copyright (c) 2011-2013 Nathan Mitchell, Eftychios Sifakis.
//  This file is covered by the FreeBSD license. Please refer to the 
//  license.txt file for more information.
//#####################################################################


/*

template<class T_RAW, class TI_RAW, class T_DATA, class TI_DATA, int width>
    void Masked_Add(T_DATA (&Result),
                    const T_DATA (&A), const T_DATA (&B),
                    const TI_DATA (&Mask));

template<class T_RAW, class TI_RAW, class T_DATA, class TI_DATA, int width>
    void Masked_Subtract(T_DATA (&Result),
                         const T_DATA (&A), const T_DATA (&B),
                         const TI_DATA (&Mask));

template<class T_RAW, class TI_RAW, class T_DATA, class TI_DATA, int width>
    void Masked_Times(T_DATA (&Result),
                      const T_DATA (&A), const T_DATA (&B),
                      const TI_DATA (&Mask));

template<class T_RAW, class TI_RAW, class T_DATA, class TI_DATA, int width>
    void Masked_SAXPY_A(T_DATA (&Result),
                        const T_DATA (&A), const T_DATA (&B), const T_DATA (&constant),
                        const TI_DATA (&Mask));

template<class T_RAW, class TI_RAW, class T_DATA, class TI_DATA, int width>
    void Masked_SAXPY_AB(T_DATA (&Result),
                         const T_DATA (&A), const T_DATA (&B), const T_DATA (&constant), const T_DATA (&B2),
                         const TI_DATA (&Mask));
*/

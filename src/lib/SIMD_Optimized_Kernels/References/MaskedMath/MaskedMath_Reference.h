//#####################################################################
//  Copyright (c) 2011-2013 Nathan Mitchell, Eftychios Sifakis.
//  This file is covered by the FreeBSD license. Please refer to the 
//  license.txt file for more information.
//#####################################################################


/*

template<class T>
    void Masked_Add(T &Result),
                    const T &A, const T &B,
                    const int &Mask);

template<class T>
    void Masked_Subtract(T (&Result),
                         const T (&A), const T (&B),
                         const TI_DATA (&Mask));

template<class T>
    void Masked_Times(T (&Result),
                      const T (&A), const T (&B),
                      const TI_DATA (&Mask));

template<class T>
    void Masked_SAXPY_A(T (&Result),
                        const T (&A), const T (&B), const T (&constant),
                        const TI_DATA (&Mask));

template<class T>
    void Masked_SAXPY_AB(T (&Result),
                         const T (&A), const T (&B), const T (&constant), const T (&B2),
                         const TI_DATA (&Mask));

/*

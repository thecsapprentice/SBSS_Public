//#####################################################################
//  Copyright (c) 2011-2013 Nathan Mitchell, Eftychios Sifakis.
//  This file is covered by the FreeBSD license. Please refer to the 
//  license.txt file for more information.
//#####################################################################

template<class T>
void Add_Force_First_Order_Corotated_Reference(const T u[3][8],
                                               const T p,
                                               
                                               const T mu,
                                               const T alpha,
                                               const T alpha_sqr_over_kappa,
                                               const T one_over_h,
                                               const T cell_volume,
                                               const T U[9],
                                               const T V[9],
                                               const T Sigma[3],
                                               const T Q_Hat[3],
                                               
                                               T f[3][8],
                                               T q);
                                               

template<class T>
void Add_Force_First_Order_Neohookean_Reference(const T u[3][8],
                                                const T p,
                                                
                                                const T mu,
                                                const T alpha,
                                                const T alpha_sqr_over_kappa,
                                                const T one_over_h,
                                                const T cell_volume,
                                                const T U[9],
                                                const T V[9],
                                                const T Sigma[3],
                                                const T Q_Hat[3],
                                                
                                                T f[3][8],
                                                T q);

template<class T>
void Add_Force_First_Order_Corotated_Compare(const T f[3][8],
                                             const T q,
                                             const T f_reference[3][8],
                                             const T q_reference);


template<class T>
void Add_Force_First_Order_Neohookean_Compare(const T f[3][8],
                                              const T q,
                                              const T f_reference[3][8],
                                              const T q_reference);





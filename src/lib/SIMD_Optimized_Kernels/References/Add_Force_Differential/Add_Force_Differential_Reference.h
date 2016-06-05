//#####################################################################
//  Copyright (c) 2011-2013 Nathan Mitchell, Eftychios Sifakis.
//  This file is covered by the FreeBSD license. Please refer to the 
//  license.txt file for more information.
//#####################################################################


template<class T> 
void Add_Force_Differential_Reference(const T du[3][8], 
                                      const T dp,
                                      const T alpha_squared_over_kappa,
                                      const T alpha,
                                      const T one_over_h,
                                      const T cell_volume,
                                      const T Q_hat[3],
                                      const T U[9],
                                      const T V[9],
                                      const T dPdF[12],
                                      
                                      T df[3][8],
                                      T &dq);


template<class T>
bool Add_Force_Differential_Compare( const T df[3][8], const T dq, const T df_reference[3][8], const T dq_reference);



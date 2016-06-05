//#####################################################################
//  Copyright (c) 2011-2013 Nathan Mitchell, Eftychios Sifakis.
//  This file is covered by the FreeBSD license. Please refer to the 
//  license.txt file for more information.
//#####################################################################


template<class T>
void Pressure_Force_Differential_Reference(T &dq, const T Q_hat[3], const T dF_hat[9], const T dp, const T alpha, const T alpha_squared_over_kappa);

template<class T>
bool Pressure_Force_Differential_Compare(const T dq, const T dq_reference);


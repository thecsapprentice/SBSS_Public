//#####################################################################
//  Copyright (c) 2011-2013 Nathan Mitchell, Eftychios Sifakis.
//  This file is covered by the FreeBSD license. Please refer to the 
//  license.txt file for more information.
//#####################################################################


template<class T>
void Penalty_Measure_Gradient_Neohookean_Reference(const T Sigma[3], T Q_hat[3]);

template<class T>
bool Penalty_Measure_Gradient_Neohookean_Compare(const T Q_hat[3], const T Q_hat_reference[3]);

template<class T>
void Penalty_Measure_Gradient_Corotated_Reference(const T Sigma[3], T Q_hat[3]);

template<class T>
bool Penalty_Measure_Gradient_Corotated_Compare(const T Q_hat[3], const T Q_hat_reference[3]);

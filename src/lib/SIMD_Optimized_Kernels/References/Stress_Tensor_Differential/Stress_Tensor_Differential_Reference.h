//#####################################################################
//  Copyright (c) 2011-2013 Nathan Mitchell, Eftychios Sifakis.
//  This file is covered by the FreeBSD license. Please refer to the 
//  license.txt file for more information.
//#####################################################################


template<class T>
    void Stress_Tensor_Differential_Reference(T dP_hat[9], const T dPdF[12], const T dF_hat[9],const T Q_hat[3], const T dp, const T alpha);

template<class T>
    bool Stress_Tensor_Differential_Compare(const T dP_hat[9], const T dP_hat_reference[9]);


//#####################################################################
//  Copyright (c) 2011-2013 Nathan Mitchell, Eftychios Sifakis.
//  This file is covered by the FreeBSD license. Please refer to the 
//  license.txt file for more information.
//#####################################################################


template<class T_MATERIAL,class T> 
struct Update_Position_Based_State_Blocked_Reference
{
static void Run(const T u[3][27], 
                                               const T p[8],
                                               const T mu[8],
                                               const T kappa[8],
                                               const T alpha[8],
                                               const T cutoff[8],
                                               const T one_over_h[8],
                                               T U[9][8],
                                               T V[9][8],
                                               T Sigma[3][8],
                                               T Q_hat[3][8],
                                               T dPdF[12][8]);
};


template<class T>
bool Update_Position_Based_State_Blocked_Compare(const T U[9][8],
                                         const T V[9][8],
                                         const T Sigma[3][8],
                                         const T Q_hat[3][8],
                                         const T dPdF[12][8],
                                         const T U_reference[9][8],
                                         const T V_reference[9][8],
                                         const T Sigma_reference[3][8],
                                         const T Q_hat_reference[3][8],
                                         const T dPdF_reference[12][8]);
                                         

struct NEOHOOKEAN_TAG;
struct COROTATED_TAG;

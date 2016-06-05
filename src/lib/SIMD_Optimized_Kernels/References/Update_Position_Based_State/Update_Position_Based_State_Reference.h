//#####################################################################
//  Copyright (c) 2011-2013 Nathan Mitchell, Eftychios Sifakis.
//  This file is covered by the FreeBSD license. Please refer to the 
//  license.txt file for more information.
//#####################################################################


template<class T> 
void Update_Position_Based_State_Corotated_Reference(
                                               const T u[3][8], 
                                               const T p,
                                               const T mu,
                                               const T mu_stab,
                                               const T kappa,
                                               const T alpha,
                                               const T cutoff,
                                               const T one_over_h,
                                               const T cell_volume,
                                               T U[9],
                                               T V[9],
                                               T Sigma[3],
                                               T Q_hat[3],
                                               T dPdF[12],
                                               T d[3][8],
                                               T system_matrix[300]);


template<class T>
bool Update_Position_Based_State_Corotated_Compare(const T U[9],
                                                   const T V[9],
                                                   const T Sigma[3],
                                                   const T Q_hat[3],
                                                   const T dPdF[12],
                                                   const T d[3][8],
                                                   const T system_matrix[300],
                                                   const T U_reference[9],
                                                   const T V_reference[9],
                                                   const T Sigma_reference[3],
                                                   const T Q_hat_reference[3],
                                                   const T dPdF_reference[12],
                                                   const T d_reference[3][8],
                                                   const T system_matrix_reference[300]);


template<class T> 
void Update_Position_Based_State_Neohookean_Reference(
                                               const T u[3][8], 
                                               const T p,
                                               const T mu,
                                               const T mu_stab,
                                               const T kappa,
                                               const T alpha,
                                               const T cutoff,
                                               const T one_over_h,
                                               const T cell_volume,
                                               T U[9],
                                               T V[9],
                                               T Sigma[3],
                                               T Q_hat[3],
                                               T dPdF[12],
                                               T d[3][8],
                                               T system_matrix[300]);



template<class T>
bool Update_Position_Based_State_Neohookean_Compare(const T U[9],
                                                    const T V[9],
                                                    const T Sigma[3],
                                                    const T Q_hat[3],
                                                    const T dPdF[12],
                                                    const T d[3][8],
                                                    const T system_matrix[300],
                                                    const T U_reference[9],
                                                    const T V_reference[9],
                                                    const T Sigma_reference[3],
                                                    const T Q_hat_reference[3],
                                                    const T dPdF_reference[12],
                                                    const T d_reference[3][8],
                                                    const T system_matrix_reference[300]);                                         

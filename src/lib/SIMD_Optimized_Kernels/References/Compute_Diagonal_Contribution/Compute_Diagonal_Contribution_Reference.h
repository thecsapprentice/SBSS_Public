

template<class T>
void  Compute_Diagonal_Contribution_Reference(const T one_over_h,
                                                       const T mu_stab,
                                                       const T cell_volume,
                                                       const T U[9],
                                                       const T V[9],
                                                       const T dPdF[12],
                                                       T d[3][8]);

template<class T>
bool  Compute_Diagonal_Contribution_Compare(const T d[3][8], const T d_reference[3][8]);
 

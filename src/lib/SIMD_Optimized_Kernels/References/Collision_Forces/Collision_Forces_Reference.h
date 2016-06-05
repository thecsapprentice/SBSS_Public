//#####################################################################
//  Copyright (c) 2011-2013 Nathan Mitchell, Eftychios Sifakis.
//  This file is covered by the FreeBSD license. Please refer to the 
//  license.txt file for more information.
//#####################################################################


template<class T>
void Collision_Forces_Reference(T df[3][8], const T du[3][8],
                                              const T N[3], const T W[3], const T h,
                                              const int spring_id, const int spring_id_X,
                                              const int spring_id_Y, const int spring_id_Z,
                                              const float* extern_collision_attach,
                                              const float* extern_collision_stiffness);

template<class T>
bool Collision_Forces_Compare(const T df[3][8], const T df_reference[3][8]);

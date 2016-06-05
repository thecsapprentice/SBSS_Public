//#####################################################################
//  Copyright (c) 2011-2013 Nathan Mitchell, Eftychios Sifakis.
//  This file is covered by the FreeBSD license. Please refer to the 
//  license.txt file for more information.
//#####################################################################


template<class T>
void Collision_Force_Differential_Reference(T df[3][8], const T du[3][8],
                                            const T W[3], const int spring_id,
                                            const float* extern_collision_stiffness);

template<class T>
bool Collision_Force_Differential_Compare(const T df[3][8], const T df_reference[3][8]);

//#####################################################################
//  Copyright (c) 2011-2013 Nathan Mitchell, Eftychios Sifakis.
//  This file is covered by the FreeBSD license. Please refer to the 
//  license.txt file for more information.
//#####################################################################


template<class T_RAW,class T_DATA=void,class I_DATA=void>
    void Collision_Force_Differential(T_DATA (&f)[3][8],
                                      const T_DATA (&u)[3][8],
                                      const T_DATA (&W)[3],
                                      const I_DATA (&spring_id),                                      
                                      const float* extern_collision_stiffness
                                      );

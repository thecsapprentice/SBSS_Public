//#####################################################################
//  Copyright (c) 2011-2013 Nathan Mitchell, Eftychios Sifakis.
//  This file is covered by the FreeBSD license. Please refer to the 
//  license.txt file for more information.
//#####################################################################


template<class T_RAW,class T_DATA=void,class I_DATA=void>
    void Expanded_Collision_Forces(T_DATA (&f)[3][8],
                                   const T_DATA (&u)[3][8],
                                   const T_DATA (&N)[3],
                                   const T_DATA (&W)[3],
                                   const T_DATA (&h),
                                   const I_DATA (&flag),
                                   const T_DATA (&stiffness_constant),
                                   const T_DATA (&attachment_point)[3],
                                   const I_DATA (&spring_id),
                                   const I_DATA (&spring_id_X),
                                   const I_DATA (&spring_id_Y),
                                   const I_DATA (&spring_id_Z),
                                   const float* extern_collision_attach,
                                   const float* extern_collision_stiffness
                                   );

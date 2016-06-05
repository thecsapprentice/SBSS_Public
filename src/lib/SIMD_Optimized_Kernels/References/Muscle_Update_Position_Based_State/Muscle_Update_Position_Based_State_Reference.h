//#####################################################################
//  Copyright (c) 2011-2013 Nathan Mitchell, Eftychios Sifakis.
//  This file is covered by the FreeBSD license. Please refer to the 
//  license.txt file for more information.
//#####################################################################


template<class T=float, class I=int > 
void Muscle_Update_Position_Based_State_Reference(const T u[3][8], 
                const I muscle_id,
                const T fiber[3],
                const T density,
                const T one_over_h,
                
                T& c1,
                T& c2,
                T F_fiber[3],
                
                const float *activations,
                const float *fiber_max_stresses);



template<class T=float, class I=int >
    bool Muscle_Update_Position_Based_State_Compare(const T c1,
                                                    const T c2,
                                                    const T F_fiber[3],
                                                    const T c1_reference,
                                                    const T c2_reference,
                                                    const T F_fiber_reference[3]
                                                    );
                                         


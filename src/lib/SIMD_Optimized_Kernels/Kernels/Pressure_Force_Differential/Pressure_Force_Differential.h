//#####################################################################
//  Copyright (c) 2011-2013 Nathan Mitchell, Eftychios Sifakis.
//  This file is covered by the FreeBSD license. Please refer to the 
//  license.txt file for more information.
//#####################################################################


template<class T_RAW,class T_DATA=void,class I_DATA=void>
    void Pressure_Force_Differential(T_DATA (&dq), const T_DATA (&Q_hat)[3], const T_DATA (&dF_hat)[9], const T_DATA (&dp), const T_DATA (&alpha), const T_DATA (&alpha_squared_over_kappa));

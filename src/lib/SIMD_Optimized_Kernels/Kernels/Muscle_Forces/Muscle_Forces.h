//#####################################################################
//  Copyright (c) 2011-2013 Nathan Mitchell, Eftychios Sifakis.
//  This file is covered by the FreeBSD license. Please refer to the 
//  license.txt file for more information.
//#####################################################################


template<class T_RAW,class T_DATA=void, class I_DATA=void>
    void Muscle_Forces(T_DATA (&f)[3][8], const T_DATA (&fiber)[3],
                       const T_DATA (&Ffiber)[3], const T_DATA (&c1),
                       const T_DATA (&one_over_h), const T_DATA (&cell_volume));

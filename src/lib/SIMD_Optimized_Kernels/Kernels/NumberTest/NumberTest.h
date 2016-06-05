//#####################################################################
//  Copyright (c) 2011-2013 Nathan Mitchell, Eftychios Sifakis.
//  This file is covered by the FreeBSD license. Please refer to the 
//  license.txt file for more information.
//#####################################################################


template<class T_RAW,class T_DATA=void,class I_DATA=void>
void NumberTest(const T_DATA (&A),
                const T_DATA (&B),
                const T_DATA (&C),
                const T_DATA (&D),

                T_DATA (&Add),
                T_DATA (&Multiply),
                T_DATA (&Subtract),
                T_DATA (&Divide),

                T_DATA (&LessThan),
                T_DATA (&GreaterThan),
                T_DATA (&LessEquals),
                T_DATA (&GreaterEquals),
                T_DATA (&Equals),

                T_DATA (&And),
                T_DATA (&Or),
                T_DATA (&Xor),
                T_DATA (&AndNot),
                T_DATA (&Not),

                T_DATA (&Sqrt),
                T_DATA (&RSqrt),
                T_DATA (&Log),
                T_DATA (&Exp),
                T_DATA (&Inverse),
                T_DATA (&AbsoluteValue),
                T_DATA (&Sign),

                T_DATA (&Minimum),
                T_DATA (&Maximum),
                T_DATA (&Blend) );

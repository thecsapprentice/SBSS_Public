//#####################################################################
//  Copyright (c) 2011-2013 Nathan Mitchell, Eftychios Sifakis.
//  This file is covered by the FreeBSD license. Please refer to the 
//  license.txt file for more information.
//#####################################################################


template<class T_RAW>
void Block_DeDuplicate(const float u_duplicated[3][8][8],float u_compact[3][27]);

template<class T_RAW>
void Block_DeDuplicate(const float u_duplicated[3][8][16],float u_compact[3][27]);

template<class T_RAW>
void Block_DeDuplicate(const float u_duplicated[3][8][16], float u_compact[2][3][27]);

//#####################################################################
//  Copyright (c) 2011-2013 Nathan Mitchell, Eftychios Sifakis.
//  This file is covered by the FreeBSD license. Please refer to the 
//  license.txt file for more information.
//#####################################################################


template <class T>
void Block_Duplicate_Reference(const T u_compact[3][27], T u_duplicated[3][8][8]);

template <class T>
bool Block_Duplicate_Compare(const T u_duplicated[3][8][8], const T u_duplicated_reference[3][8][8]);

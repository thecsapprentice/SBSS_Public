//#####################################################################
//  Copyright (c) 2011-2013 Nathan Mitchell, Eftychios Sifakis.
//  This file is covered by the FreeBSD license. Please refer to the 
//  license.txt file for more information.
//#####################################################################


template <class T>
void Block_DeDuplicate_Reference(const T u_duplicated[3][8][8],T u_compact[3][27]);

template <class T>
bool Block_DeDuplicate_Compare(const T u_compact[3][27],const T u_compact_reference[3][27]);

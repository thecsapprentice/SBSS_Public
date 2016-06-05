//#####################################################################
//  Copyright (c) 2011-2013 Nathan Mitchell, Eftychios Sifakis.
//  This file is covered by the FreeBSD license. Please refer to the 
//  license.txt file for more information.
//#####################################################################


template<class T>
void RSD_Positive_Definite_Part_Reference(const T RSD[12], T RSDpd[12]);

template<class T>
bool RSD_Positive_Definite_Part_Compare(const T RSDpd[12], const T RSDpd_reference[12]);

//#####################################################################
//  Copyright (c) 2011-2013 Nathan Mitchell, Eftychios Sifakis.
//  This file is covered by the FreeBSD license. Please refer to the 
//  license.txt file for more information.
//#####################################################################


#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Arrays_Computations/MAGNITUDE.h>
#include <PhysBAM_Tools/Matrices/DIAGONAL_MATRIX_3X3.h>
#include <PhysBAM_Tools/Matrices/MATRIX_3X3.h>
#include <PhysBAM_Tools/Matrices/SYMMETRIC_MATRIX_3X3.h>

#include "ROTATED_STRESS_DERIVATIVE.h"

using namespace PhysBAM;

template<class T>
void RSD_Positive_Definite_Part_Reference(const T RSD[12], T RSDpd[12])
{
    const ROTATED_STRESS_DERIVATIVE<T,3>& rRSD=*(const ROTATED_STRESS_DERIVATIVE<T,3>*)(RSD);
    ROTATED_STRESS_DERIVATIVE<T,3>& rRSDpd=*(ROTATED_STRESS_DERIVATIVE<T,3>*)(RSDpd);
    rRSDpd=rRSD;
    rRSDpd.Make_Positive_Definite();
}

template<class T>
bool RSD_Positive_Definite_Part_Compare(const T RSDpd[12], const T RSDpd_reference[12])
{
    ARRAY_VIEW<const T> aRSDpd(12,RSDpd);
    ARRAY_VIEW<const T> aRSDpd_reference(12,RSDpd_reference);

    std::cout<<"Computed RSDpd : "<<aRSDpd<<std::endl;
    std::cout<<"Reference RSDpd :"<<aRSDpd_reference<<std::endl;
    ARRAY<T> difference(aRSDpd-aRSDpd_reference);
    std::cout<<"Difference = "<<ARRAYS_COMPUTATIONS::Maximum_Magnitude(difference)<<std::endl;

    if( ARRAYS_COMPUTATIONS::Maximum_Magnitude(difference) < 0.00001 )
        return true;
    else
        return false;
}

template void RSD_Positive_Definite_Part_Reference(const float RSD[12], float RSDpd[12]);
template bool RSD_Positive_Definite_Part_Compare(const float RSDpd[12], const float RSDpd_reference[12]);

 

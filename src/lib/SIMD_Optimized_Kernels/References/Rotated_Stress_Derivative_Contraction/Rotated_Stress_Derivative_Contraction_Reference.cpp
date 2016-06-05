//#####################################################################
//  Copyright (c) 2011-2013 Nathan Mitchell, Eftychios Sifakis.
//  This file is covered by the FreeBSD license. Please refer to the 
//  license.txt file for more information.
//#####################################################################


#include <PhysBAM_Tools/Matrices/MATRIX_3X3.h>

#include "ROTATED_STRESS_DERIVATIVE.h"

using namespace PhysBAM;

namespace{
template<class T_MATRIX>
void Print_Formatted(const T_MATRIX& A,std::ostream& output)
{
    for(int i=1;i<=A.Rows();i++){
        for(int j=1;j<=A.Columns();j++){
            if(A.Valid_Index(i,j))
                output<<std::setw(12)<<A(i,j);
            else
                output<<"            ";
            if(j<A.Columns()) output<<" ";}
        output<<std::endl;}
}
}
 

template<class T>
void Rotated_Stress_Derivative_Contraction_Reference(const T dPdF[12], const T dF_Hat[9], T dP_Hat[9])
{
    const ROTATED_STRESS_DERIVATIVE<T,3>& mdPdF=*(const ROTATED_STRESS_DERIVATIVE<T,3>*)(dPdF);
    const MATRIX<T,3>& mdF_Hat=*(const MATRIX<T,3>*)(dF_Hat);
    MATRIX<T,3>& mdP_Hat=*(MATRIX<T,3>*)(dP_Hat);

    mdP_Hat=mdPdF.dP_hat(mdF_Hat);
}

template<class T>
bool Rotated_Stress_Derivative_Contraction_Compare(const T dP_Hat[9], const T dP_Hat_reference[9])
{
    const MATRIX<T,3>& mdP_Hat=*(const MATRIX<T,3>*)(dP_Hat);
    const MATRIX<T,3>& mdP_Hat_reference=*(const MATRIX<T,3>*)(dP_Hat_reference);

    std::cout<<"Computed matrix dP_Hat :"<<std::endl;Print_Formatted(mdP_Hat,std::cout);
    std::cout<<"Reference matrix dP_Hat :"<<std::endl;Print_Formatted(mdP_Hat_reference,std::cout);
    std::cout<<"Difference = "<<(mdP_Hat-mdP_Hat_reference).Frobenius_Norm()<<std::endl;

    if( (mdP_Hat-mdP_Hat_reference).Frobenius_Norm() < 0.00001 )
        return true;
    else
        return false;
}

template void Rotated_Stress_Derivative_Contraction_Reference(const float dPdF[12], const float dF_Hat[9], float dP_Hat[9]);
template bool Rotated_Stress_Derivative_Contraction_Compare(const float dP_Hat[9], const float dP_Hat_reference[9]);
 

//#####################################################################
//  Copyright (c) 2011-2013 Nathan Mitchell, Eftychios Sifakis.
//  This file is covered by the FreeBSD license. Please refer to the 
//  license.txt file for more information.
//#####################################################################


#include <PhysBAM_Tools/Matrices/DIAGONAL_MATRIX_3X3.h>
#include <PhysBAM_Tools/Matrices/MATRIX_3X3.h>
#include <PhysBAM_Tools/Matrices/SYMMETRIC_MATRIX_3X3.h>
#include <iomanip>

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
void Symmetric_Definite_Projection_Reference(const T A[6], T Apd[6])
{
    const SYMMETRIC_MATRIX<T,3>& mA=*(const SYMMETRIC_MATRIX<T,3>*)(A);
    SYMMETRIC_MATRIX<T,3>& mApd=*(SYMMETRIC_MATRIX<T,3>*)(Apd);

    mApd=mA.Positive_Definite_Part();
}

template<class T>
bool Symmetric_Definite_Projection_Compare(const T Apd[6], const T Apd_reference[6])
{
    const SYMMETRIC_MATRIX<T,3>& mApd=*(const SYMMETRIC_MATRIX<T,3>*)(Apd);
    const SYMMETRIC_MATRIX<T,3>& mApd_reference=*(const SYMMETRIC_MATRIX<T,3>*)(Apd_reference);

    std::cout<<"Computed matrix Apd :"<<std::endl;Print_Formatted(mApd,std::cout);
    std::cout<<"Reference matrix Apd :"<<std::endl;Print_Formatted(mApd_reference,std::cout);
    std::cout<<"Difference = "<<(mApd-mApd_reference).Frobenius_Norm()<<std::endl;

    if( (mApd-mApd_reference).Frobenius_Norm() < 0.00001 )
        return true;
    else
        return false;
}

template void Symmetric_Definite_Projection_Reference(const float A[6], float Apd[6]);
template bool Symmetric_Definite_Projection_Compare(const float Apd[6], const float Apd_reference[6]);
 

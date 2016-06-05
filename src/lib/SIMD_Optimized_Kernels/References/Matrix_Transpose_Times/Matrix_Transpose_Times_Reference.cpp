//#####################################################################
//  Copyright (c) 2011-2013 Nathan Mitchell, Eftychios Sifakis.
//  This file is covered by the FreeBSD license. Please refer to the 
//  license.txt file for more information.
//#####################################################################


#include <PhysBAM_Tools/Matrices/MATRIX_3X3.h>

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
void Matrix_Transpose_Times_Reference(const T A[9], const T B[9], T C[9])
{
    const MATRIX<T,3>& mA=*(const MATRIX<T,3>*)(A);
    const MATRIX<T,3>& mB=*(const MATRIX<T,3>*)(B);
    MATRIX<T,3>& mC=*(MATRIX<T,3>*)(C);

    mC=mA.Transpose_Times(mB);
}

template<class T>
bool Matrix_Transpose_Times_Compare(const T C[9], const T C_reference[9])
{
    const MATRIX<T,3>& mC=*(const MATRIX<T,3>*)(C);
    const MATRIX<T,3>& mC_reference=*(const MATRIX<T,3>*)(C_reference);

    std::cout<<"Computed matrix C :"<<std::endl;Print_Formatted(mC,std::cout);
    std::cout<<"Reference matrix C :"<<std::endl;Print_Formatted(mC_reference,std::cout);
    std::cout<<"Difference = "<<(mC-mC_reference).Frobenius_Norm()<<std::endl;

    if( (mC-mC_reference).Frobenius_Norm() < 0.00001 )
        return true;
    else
        return false;
}

template void Matrix_Transpose_Times_Reference(const float A[9], const float B[9], float C[9]);
template bool Matrix_Transpose_Times_Compare(const float C[9], const float C_reference[9]);
 

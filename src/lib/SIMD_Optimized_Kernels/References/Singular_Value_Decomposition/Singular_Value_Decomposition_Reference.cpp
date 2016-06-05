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
void Singular_Value_Decomposition_Reference(const T A[9], T U[9], T Sigma[3], T V[9])
{
    const MATRIX<T,3>& mA=*(const MATRIX<T,3>*)(A);

    MATRIX<T,3>& mU=*(MATRIX<T,3>*)(U);
    DIAGONAL_MATRIX<T,3>& mSigma=*(DIAGONAL_MATRIX<T,3>*)(Sigma);
    MATRIX<T,3>& mV=*(MATRIX<T,3>*)(V);

    mA.Fast_Singular_Value_Decomposition(mU,mSigma,mV);
}

template<class T>
bool Singular_Value_Decomposition_Compare(const T U[9], const T Sigma[3], const T V[9],
                                          const T U_reference[9], const T Sigma_reference[3], const T V_reference[9])
{

    MATRIX<T,3>& mU=*(MATRIX<T,3>*)(U);
    DIAGONAL_MATRIX<T,3>& mSigma=*(DIAGONAL_MATRIX<T,3>*)(Sigma);
    MATRIX<T,3>& mV=*(MATRIX<T,3>*)(V);

    MATRIX<T,3>& mU_reference=*(MATRIX<T,3>*)(U_reference);
    DIAGONAL_MATRIX<T,3>& mSigma_reference=*(DIAGONAL_MATRIX<T,3>*)(Sigma_reference);
    MATRIX<T,3>& mV_reference=*(MATRIX<T,3>*)(V_reference);

    std::cout<<"Computed matrix U :"<<std::endl;Print_Formatted(mU,std::cout);
    std::cout<<"Reference matrix U :"<<std::endl;Print_Formatted(mU_reference,std::cout);
    std::cout<<"Difference = "<< (mU-mU_reference).Frobenius_Norm()  <<std::endl;

    std::cout<<std::endl;
    std::cout<<"Computed matrix Sigma :"<<std::endl;Print_Formatted(mSigma,std::cout);
    std::cout<<"Reference matrix Sigma :"<<std::endl;Print_Formatted(mSigma_reference,std::cout);
    std::cout<<"Difference = "<<  (mSigma-mSigma_reference).Frobenius_Norm()  <<std::endl;
    
    std::cout<<std::endl;
    std::cout<<"Computed matrix V :"<<std::endl;Print_Formatted(mV,std::cout);
    std::cout<<"Reference matrix V :"<<std::endl;Print_Formatted(mV_reference,std::cout);
    std::cout<<"Difference = "<< (mV-mV_reference).Frobenius_Norm()  <<std::endl;   

    MATRIX<T,3> mA=mU*mSigma.Times_Transpose(mV);
    MATRIX<T,3> mA_reference=mU_reference*mSigma_reference.Times_Transpose(mV_reference);

    std::cout<<std::endl;
    std::cout<<"Computed matrix A=U*Sigma*V^T :"<<std::endl;Print_Formatted(mA,std::cout);
    std::cout<<"Reference matrix A=U*Sigma*V^T :"<<std::endl;Print_Formatted(mA_reference,std::cout);
    std::cout<<"Difference = "<< (mA-mA_reference).Frobenius_Norm()  <<std::endl;


    if(!((mU-mU_reference).Frobenius_Norm() < 0.00001) )
        return false;
    

    if(!((mSigma-mSigma_reference).Frobenius_Norm() < 0.00001) )
        return false;
    

    if(!((mV-mV_reference).Frobenius_Norm() < 0.00001) )
        return false;
    

    






    return true;
}

template void Singular_Value_Decomposition_Reference(const float A[9], float U[9], float Sigma[3], float V[9]);
template bool Singular_Value_Decomposition_Compare(const float U[9], const float Sigma[3], const float V[9],
                                                   const float U_reference[9], const float Sigma_reference[3], const float V_reference[9]);
 

//#####################################################################
//  Copyright (c) 2011-2013 Nathan Mitchell, Eftychios Sifakis.
//  This file is covered by the FreeBSD license. Please refer to the 
//  license.txt file for more information.
//#####################################################################


#include "NEOHOOKEAN.h"
#include "COROTATED.h"

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
void Penalty_Measure_Gradient_Corotated_Reference(const T Sigma[3], T Q_hat[3])
{
    const DIAGONAL_MATRIX<T,3>& mSigma=*(const DIAGONAL_MATRIX<T,3>*)(Sigma);    
    DIAGONAL_MATRIX<T,3>& mQ_hat=*(DIAGONAL_MATRIX<T,3>*)(Q_hat);    
    mQ_hat=COROTATED<T,3>::Q_hat(mSigma);
}

template<class T>
bool Penalty_Measure_Gradient_Corotated_Compare(const T Q_hat[3], const T Q_hat_reference[3])
{
    const DIAGONAL_MATRIX<T,3>& mQ_hat=*(const DIAGONAL_MATRIX<T,3>*)(Q_hat);    
    const DIAGONAL_MATRIX<T,3>& mQ_hat_reference=*(const DIAGONAL_MATRIX<T,3>*)(Q_hat_reference);    

    std::cout<<"Computed matrix Q_hat :"<<std::endl;Print_Formatted(mQ_hat,std::cout);
    std::cout<<"Reference matrix Q_hat :"<<std::endl;Print_Formatted(mQ_hat_reference,std::cout);
    std::cout<<"Difference = "<<(mQ_hat-mQ_hat_reference).Frobenius_Norm()<<std::endl;

    if( (mQ_hat-mQ_hat_reference).Frobenius_Norm() < 0.00001 )
        return true;
    else
        return false;
}

template void Penalty_Measure_Gradient_Corotated_Reference(const float Sigma[3], float Q_hat[3]);
template bool Penalty_Measure_Gradient_Corotated_Compare(const float Q_hat[3], const float Q_hat_reference[3]);



template<class T>
void Penalty_Measure_Gradient_Neohookean_Reference(const T Sigma[3], T Q_hat[3])
{
    const DIAGONAL_MATRIX<T,3>& mSigma=*(const DIAGONAL_MATRIX<T,3>*)(Sigma);    
    DIAGONAL_MATRIX<T,3>& mQ_hat=*(DIAGONAL_MATRIX<T,3>*)(Q_hat);    
    mQ_hat=NEOHOOKEAN<T,3>::Q_hat(mSigma);
}

template<class T>
bool Penalty_Measure_Gradient_Neohookean_Compare(const T Q_hat[3], const T Q_hat_reference[3])
{
    const DIAGONAL_MATRIX<T,3>& mQ_hat=*(const DIAGONAL_MATRIX<T,3>*)(Q_hat);    
    const DIAGONAL_MATRIX<T,3>& mQ_hat_reference=*(const DIAGONAL_MATRIX<T,3>*)(Q_hat_reference);    

    std::cout<<"Computed matrix Q_hat :"<<std::endl;Print_Formatted(mQ_hat,std::cout);
    std::cout<<"Reference matrix Q_hat :"<<std::endl;Print_Formatted(mQ_hat_reference,std::cout);
    std::cout<<"Difference = "<<(mQ_hat-mQ_hat_reference).Frobenius_Norm()<<std::endl;

    if( (mQ_hat-mQ_hat_reference).Frobenius_Norm() < 0.00001 )
        return true;
    else
        return false;
}

template void Penalty_Measure_Gradient_Neohookean_Reference(const float Sigma[3], float Q_hat[3]);
template bool Penalty_Measure_Gradient_Neohookean_Compare(const float Q_hat[3], const float Q_hat_reference[3]);

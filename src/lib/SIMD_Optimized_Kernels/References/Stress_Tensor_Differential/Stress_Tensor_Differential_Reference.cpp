//#####################################################################
//  Copyright (c) 2011-2013 Nathan Mitchell, Eftychios Sifakis.
//  This file is covered by the FreeBSD license. Please refer to the 
//  license.txt file for more information.
//#####################################################################


#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Arrays_Computations/MAGNITUDE.h>
#include <PhysBAM_Tools/Matrices/MATRIX_2X2.h>
#include <PhysBAM_Tools/Matrices/MATRIX_3X3.h>

#include "COROTATED.h"
#include "ROTATED_STRESS_DERIVATIVE.h"
#include "MATERIAL_MODEL.h"

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
void Stress_Tensor_Differential_Reference(T dP_hat[9], const T dPdF[12], const T dF_hat[9],const T Q_hat[3], const T dp, const T alpha)
{
    MATRIX<T,3>& mdP_hat=*(MATRIX<T,3>*)(dP_hat);    
    const ROTATED_STRESS_DERIVATIVE<T,3>& rdPdF=*(const ROTATED_STRESS_DERIVATIVE<T,3>*)(dPdF);
    const MATRIX<T,3>& mdF_hat=*(const MATRIX<T,3>*)(dF_hat);    
    const DIAGONAL_MATRIX<T,3>& mQhat=*(const DIAGONAL_MATRIX<T,3>*)(Q_hat);    
    
    mdP_hat=MATERIAL_MODEL<COROTATED<T,3> >::dP_hat(rdPdF,mQhat,mdF_hat,dp,alpha);
}

template<class T>
bool Stress_Tensor_Differential_Compare(const T dP_hat[9], const T dP_hat_reference[9])
{
    const MATRIX<T,3>& mdP_hat=*(const MATRIX<T,3>*)(dP_hat);
    const MATRIX<T,3>& mdP_hat_reference=*(const MATRIX<T,3>*)(dP_hat_reference);

    std::cout<<"Computed matrix dP_hat :"<<std::endl;Print_Formatted(mdP_hat,std::cout);
    std::cout<<"Reference matrix dP_hat :"<<std::endl;Print_Formatted(mdP_hat_reference,std::cout);
    std::cout<<"Difference = "<<(mdP_hat-mdP_hat_reference).Frobenius_Norm()<<std::endl;

    if( (mdP_hat-mdP_hat_reference).Frobenius_Norm() < 0.00001 )
        return true;
    else
        return false;
}

template void Stress_Tensor_Differential_Reference(float dP_hat[9], const float dPdF[12], const float dF_hat[9],const float Q_hat[3], const float p, const float alpha);
template bool Stress_Tensor_Differential_Compare(const float dP_hat[9], const float dP_hat_reference[9]);

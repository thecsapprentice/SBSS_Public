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
void Pressure_Force_Differential_Reference(T &dq, const T Q_hat[3], const T dF_hat[9], const T dp, const T alpha, const T alpha_squared_over_kappa)
{
    const MATRIX<T,3>& mdF_hat=*(const MATRIX<T,3>*)(dF_hat);    
    const DIAGONAL_MATRIX<T,3>& mQ_hat=*(const DIAGONAL_MATRIX<T,3>*)(Q_hat);    
    dq = -alpha*mQ_hat.Times_Transpose(mdF_hat).Trace()+alpha_squared_over_kappa*dp;
}



template<class T>
bool  Pressure_Force_Differential_Compare(const T dq, const T dq_reference)
{

    std::cout<<"Computed  dq:"<<std::endl << dq << std::endl;
    std::cout<<"Reference dq:"<<std::endl << dq_reference <<std::endl;
    std::cout<<"Difference = "<<(dq-dq_reference)<<std::endl;

    if( (dq-dq_reference) < 0.00001 )
        return true;
    else
        return false;
}

template void Pressure_Force_Differential_Reference(float &dq, const float Q_hat[3], const float dF_hat[9], const float dp, const float alpha, const float alpha_squared_over_kappa);
template bool Pressure_Force_Differential_Compare(const float dq, const float dq_reference);

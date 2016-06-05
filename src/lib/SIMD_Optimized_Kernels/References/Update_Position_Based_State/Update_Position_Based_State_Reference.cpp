//#####################################################################
//  Copyright (c) 2011-2013 Nathan Mitchell, Eftychios Sifakis.
//  This file is covered by the FreeBSD license. Please refer to the 
//  license.txt file for more information.
//#####################################################################


#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Arrays_Computations/MAGNITUDE.h>
#include <PhysBAM_Tools/Matrices/MATRIX_3X3.h>
#include <PhysBAM_Tools/Matrices/SYMMETRIC_MATRIX_3X3.h>
#include <PhysBAM_Tools/Matrices/MATRIX_MXN.h>

#include "COROTATED.h"
#include "NEOHOOKEAN.h"
#include "MATERIAL_MODEL.h"
#include "RANGE_ITERATOR.h"

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
template<class T,int d> void
Gradient_Matrix(MATRIX_MXN<T>& G,const VECTOR<T,d>& weights,const T one_over_h)
{
    typedef VECTOR<int,d> T_INDEX;

    G.Resize(d,(d==2)?4:8);

    int vertex=1;
    for(RANGE_ITERATOR<d> iterator(RANGE<T_INDEX>(T_INDEX(),T_INDEX::All_Ones_Vector()));iterator.Valid();iterator.Next(),vertex++){
        const T_INDEX& index=iterator.Index();
        
        for(int v=1;v<=d;v++){
            T& coefficient=G(v,vertex);
            coefficient=1.;
            
            for(int w=1;w<=d;w++)
                if(w==v)
                    if(index(w)==0) coefficient*=-one_over_h; else coefficient*=one_over_h;
                else
                    if(index(w)==0) coefficient*=(1.-weights(w)); else coefficient*=weights(w);}}
}
}

template<class T> 
void Update_Position_Based_State_Corotated_Reference(
                                                     const T u[3][8], 
                                                     const T p,
                                                     const T mu,
                                                     const T mu_stab,
                                                     const T kappa,
                                                     const T alpha,
                                                     const T cutoff,
                                                     const T one_over_h,
                                                     const T cell_volume,
                                                     T U[9],
                                                     T V[9],
                                                     T Sigma[3],
                                                     T Q_hat[3],
                                                     T dPdF[12],
                                                     T d[3][8],
                                                     T system_matrix[300])
{   
    MATRIX_MXN<float> G;
    MATRIX<float,3> F;
    MATRIX<float,3>& mU=*(MATRIX<float,3>*)(U);
    MATRIX<float,3>& mV=*(MATRIX<float,3>*)(V);
    ROTATED_STRESS_DERIVATIVE<float,3>& rdPdF=*(ROTATED_STRESS_DERIVATIVE<float,3>*)(dPdF);
    DIAGONAL_MATRIX<float,3>& mSigma=*(DIAGONAL_MATRIX<float,3>*)(Sigma);
    DIAGONAL_MATRIX<float,3>& mQ_hat=*(DIAGONAL_MATRIX<float,3>*)(Q_hat);
    MATRIX_MXN<float> Du(3,8);

    const VECTOR<float,3> weights=VECTOR<float,3>::All_Ones_Vector()*(float).5;
    Gradient_Matrix(G,weights,one_over_h);

    for(int i=0;i<3;i++) for(int j=0;j<8;j++) Du(i+1,j+1)=u[i][j];
    MATRIX_MXN<float> mF=Du*G.Transposed();
    for(int i=0;i<9;i++) F.x[i]=mF.x[i];
    F+=1;

    F.Fast_Singular_Value_Decomposition(mU,mSigma,mV);
    mSigma=mSigma.Clamp_Min(cutoff);

    mQ_hat=COROTATED<T,3>::Q_hat(mSigma);    
    rdPdF=MATERIAL_MODEL<COROTATED<T,3> >::Isotropic_Stress_Derivative(mSigma,p,mu,kappa,alpha,true);
}



template<class T>
bool Update_Position_Based_State_Corotated_Compare(const T U[9],
                                                   const T V[9],
                                                   const T Sigma[3],
                                                   const T Q_hat[3],
                                                   const T dPdF[12],
                                                   const T d[3][8],
                                                   const T system_matrix[300],
                                                   const T U_reference[9],
                                                   const T V_reference[9],
                                                   const T Sigma_reference[3],
                                                   const T Q_hat_reference[3],
                                                   const T dPdF_reference[12],
                                                   const T d_reference[3][8],
                                                   const T system_matrix_reference[300]
                                                   )
{
    const DIAGONAL_MATRIX<T,3>& mSigma=*(const DIAGONAL_MATRIX<T,3>*)(Sigma);
    const DIAGONAL_MATRIX<T,3>& mSigma_reference=*(const DIAGONAL_MATRIX<T,3>*)(Sigma_reference);

    std::cout<<"Computed matrix Sigma :"<<std::endl;Print_Formatted(mSigma,std::cout);
    std::cout<<"Reference matrix Sigma :"<<std::endl;Print_Formatted(mSigma_reference,std::cout);
    std::cout<<"Difference = "<<(mSigma-mSigma_reference).Frobenius_Norm()<<std::endl;

    const DIAGONAL_MATRIX<T,3>& mQ_hat=*(const DIAGONAL_MATRIX<T,3>*)(Q_hat);
    const DIAGONAL_MATRIX<T,3>& mQ_hat_reference=*(const DIAGONAL_MATRIX<T,3>*)(Q_hat_reference);

    std::cout<<"Computed matrix Q_hat :"<<std::endl;Print_Formatted(mQ_hat,std::cout);
    std::cout<<"Reference matrix Q_hat :"<<std::endl;Print_Formatted(mQ_hat_reference,std::cout);
    std::cout<<"Difference = "<<(mQ_hat-mQ_hat_reference).Frobenius_Norm()<<std::endl;

    ARRAY_VIEW<const T> adPdF(12,dPdF);
    ARRAY_VIEW<const T> adPdF_reference(12,dPdF_reference);

    std::cout<<"Computed dPdF : "<<adPdF<<std::endl;
    std::cout<<"Reference dPdF :"<<adPdF_reference<<std::endl;
    ARRAY<T> difference(adPdF-adPdF_reference);
    std::cout<<"Difference = "<<ARRAYS_COMPUTATIONS::Maximum_Magnitude(difference)<<std::endl;

    if( (mSigma-mSigma_reference).Frobenius_Norm() < 0.00001 )
        return true;
    else
        return false;

    if( (mQ_hat-mQ_hat_reference).Frobenius_Norm() < 0.00001 )
        return true;
    else
        return false;

    if( ARRAYS_COMPUTATIONS::Maximum_Magnitude(difference) < 0.00001 )
        return true;
    else
        return false;
}



template void Update_Position_Based_State_Corotated_Reference<float>(const float u[3][8], 
                                                                     const float p,
                                                                     const float mu,
                                                                     const float mu_stab,
                                                                     const float kappa,
                                                                     const float alpha,
                                                                     const float cutoff,
                                                                     const float one_over_h,
                                                                     const float cell_volume,
                                                                     float U[9],
                                                                     float V[9],
                                                                     float Sigma[3], 
                                                                     float Q_hat[3],
                                                                     float dPdF[12],
                                                                     float d[3][8],
                                                                     float system_matrix[300] );
template bool Update_Position_Based_State_Corotated_Compare(const float U[9],
                                                            const float V[9],
                                                            const float Sigma[3],
                                                            const float Q_hat[3],
                                                            const float dPdF[12],
                                                            const float d[3][8],
                                                            const float system_matrix[300],
                                                            const float U_reference[9],
                                                            const float V_reference[9],
                                                            const float Sigma_reference[3],
                                                            const float Q_hat_reference[3],
                                                            const float dPdF_reference[12],
                                                            const float d_reference[3][8],
                                                            const float system_matrix_reference[300]
                                                            );


template<class T> 
void Update_Position_Based_State_Neohookean_Reference(
                                                     const T u[3][8], 
                                                     const T p,
                                                     const T mu,
                                                     const T mu_stab,
                                                     const T kappa,
                                                     const T alpha,
                                                     const T cutoff,
                                                     const T one_over_h,
                                                     const T cell_volume,
                                                     T U[9],
                                                     T V[9],
                                                     T Sigma[3],
                                                     T Q_hat[3],
                                                     T dPdF[12],
                                                     T d[3][8],
                                                     T system_matrix[300])
{   
    MATRIX_MXN<float> G;
    MATRIX<float,3> F;
    MATRIX<float,3>& mU=*(MATRIX<float,3>*)(U);
    MATRIX<float,3>& mV=*(MATRIX<float,3>*)(V);
    ROTATED_STRESS_DERIVATIVE<float,3>& rdPdF=*(ROTATED_STRESS_DERIVATIVE<float,3>*)(dPdF);
    DIAGONAL_MATRIX<float,3>& mSigma=*(DIAGONAL_MATRIX<float,3>*)(Sigma);
    DIAGONAL_MATRIX<float,3>& mQ_hat=*(DIAGONAL_MATRIX<float,3>*)(Q_hat);
    MATRIX_MXN<float> Du(3,8);

    const VECTOR<float,3> weights=VECTOR<float,3>::All_Ones_Vector()*(float).5;
    Gradient_Matrix(G,weights,one_over_h);

    for(int i=0;i<3;i++) for(int j=0;j<8;j++) Du(i+1,j+1)=u[i][j];
    MATRIX_MXN<float> mF=Du*G.Transposed();
    for(int i=0;i<9;i++) F.x[i]=mF.x[i];
    F+=1;

    F.Fast_Singular_Value_Decomposition(mU,mSigma,mV);
    mSigma=mSigma.Clamp_Min(cutoff);

    mQ_hat=NEOHOOKEAN<T,3>::Q_hat(mSigma);    
    rdPdF=MATERIAL_MODEL<NEOHOOKEAN<T,3> >::Isotropic_Stress_Derivative(mSigma,p,mu,kappa,alpha,true);
}



template<class T>
bool Update_Position_Based_State_Neohookean_Compare(const T U[9],
                                                    const T V[9],
                                                    const T Sigma[3],
                                                    const T Q_hat[3],
                                                    const T dPdF[12],
                                                    const T d[3][8],
                                                    const T system_matrix[300],
                                                    const T U_reference[9],
                                                    const T V_reference[9],
                                                    const T Sigma_reference[3],
                                                    const T Q_hat_reference[3],
                                                    const T dPdF_reference[12],
                                                    const T d_reference[3][8],
                                                    const T system_matrix_reference[300]
                                                    )
{
    const DIAGONAL_MATRIX<T,3>& mSigma=*(const DIAGONAL_MATRIX<T,3>*)(Sigma);
    const DIAGONAL_MATRIX<T,3>& mSigma_reference=*(const DIAGONAL_MATRIX<T,3>*)(Sigma_reference);

    std::cout<<"Computed matrix Sigma :"<<std::endl;Print_Formatted(mSigma,std::cout);
    std::cout<<"Reference matrix Sigma :"<<std::endl;Print_Formatted(mSigma_reference,std::cout);
    std::cout<<"Difference = "<<(mSigma-mSigma_reference).Frobenius_Norm()<<std::endl;

    const DIAGONAL_MATRIX<T,3>& mQ_hat=*(const DIAGONAL_MATRIX<T,3>*)(Q_hat);
    const DIAGONAL_MATRIX<T,3>& mQ_hat_reference=*(const DIAGONAL_MATRIX<T,3>*)(Q_hat_reference);

    std::cout<<"Computed matrix Q_hat :"<<std::endl;Print_Formatted(mQ_hat,std::cout);
    std::cout<<"Reference matrix Q_hat :"<<std::endl;Print_Formatted(mQ_hat_reference,std::cout);
    std::cout<<"Difference = "<<(mQ_hat-mQ_hat_reference).Frobenius_Norm()<<std::endl;

    ARRAY_VIEW<const T> adPdF(12,dPdF);
    ARRAY_VIEW<const T> adPdF_reference(12,dPdF_reference);

    std::cout<<"Computed dPdF : "<<adPdF<<std::endl;
    std::cout<<"Reference dPdF :"<<adPdF_reference<<std::endl;
    ARRAY<T> difference(adPdF-adPdF_reference);
    std::cout<<"Difference = "<<ARRAYS_COMPUTATIONS::Maximum_Magnitude(difference)<<std::endl;

    if( (mSigma-mSigma_reference).Frobenius_Norm() < 0.00001 )
        return true;
    else
        return false;

    if( (mQ_hat-mQ_hat_reference).Frobenius_Norm() < 0.00001 )
        return true;
    else
        return false;

    if( ARRAYS_COMPUTATIONS::Maximum_Magnitude(difference) < 0.00001 )
        return true;
    else
        return false;
}



template void Update_Position_Based_State_Neohookean_Reference<float>(const float u[3][8], 
                                                                      const float p,
                                                                      const float mu,
                                                                      const float mu_stab,
                                                                      const float kappa,
                                                                      const float alpha,
                                                                      const float cutoff,
                                                                      const float one_over_h,
                                                                      const float cell_volume,
                                                                      float U[9],
                                                                      float V[9],
                                                                      float Sigma[3],
                                                                      float Q_hat[3],
                                                                      float dPdF[12],
                                                                      float d[3][8],
                                                                      float system_matrix[300] );

template bool Update_Position_Based_State_Neohookean_Compare(const float U[9],
                                                            const float V[9],
                                                            const float Sigma[3],
                                                            const float Q_hat[3],
                                                            const float dPdF[12],
                                                            const float d[3][8],
                                                            const float system_matrix[300],
                                                            const float U_reference[9],
                                                            const float V_reference[9],
                                                            const float Sigma_reference[3],
                                                            const float Q_hat_reference[3],
                                                            const float dPdF_reference[12],
                                                            const float d_reference[3][8],
                                                            const float system_matrix_reference[300]
                                                            );

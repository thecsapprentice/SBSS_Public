//#####################################################################
//  Copyright (c) 2011-2013 Nathan Mitchell, Eftychios Sifakis.
//  This file is covered by the FreeBSD license. Please refer to the 
//  license.txt file for more information.
//#####################################################################


struct COROTATED_TAG;
struct NEOHOOKEAN_TAG;

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

template<class T_TAG> struct TAG_TO_MODEL;
template<> struct TAG_TO_MODEL<COROTATED_TAG> {typedef COROTATED<float,3> TYPE;};
template<> struct TAG_TO_MODEL<NEOHOOKEAN_TAG> {typedef NEOHOOKEAN<float,3> TYPE;};

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

template<class T_MATERIAL,class T> 
struct Update_Position_Based_State_Blocked_Reference
{
static void Run(const T u[3][27], 
                                               const T p[8],
                                               const T mu[8],
                                               const T kappa[8],
                                               const T alpha[8],
                                               const T cutoff[8],
                                               const T one_over_h[8],
                                               T U[9][8],
                                               T V[9][8],
                                               T Sigma[3][8],
                                               T Q_hat[3][8],
                                               T dPdF[12][8])
{   



    typedef float T3333[3][3][3][3];

    const T3333& u_compact_3333 = *(reinterpret_cast<const T3333*>(u));
    float u_d[3][8][8];


    for( int v=0; v<3; v++)
        for( int i_base=0; i_base<2; i_base++)
            for( int j_base=0; j_base<2; j_base++)
                for( int k_base=0; k_base<2; k_base++)
                    for( int i=0; i<2; i++)
                        for( int j=0; j<2; j++)
                            for( int k=0; k<2; k++)
                                {
                                    u_d[v][i*4+j*2+k][i_base*4+j_base*2+k_base] =
                                        u_compact_3333[v][i+i_base][j+j_base][k+k_base];
                                }

    for( int i = 0; i < 8; i++)
        {
                     


            MATRIX_MXN<float> G;
            MATRIX<float,3> F;
            MATRIX<float,3> mU;
            for(int k=0;k<3;k++) for(int j=0;j<3;j++) mU(k+1,j+1) = U[k*3+j][i];
            MATRIX<float,3> mV;
            for(int k=0;k<3;k++) for(int j=0;j<3;j++) mV(k+1,j+1) = V[k*3+j][i];
            ROTATED_STRESS_DERIVATIVE<float,3> rdPdF;
            for(int k=0;k<12;k++) *(&(rdPdF.a1111)+k) = dPdF[k][i];
            DIAGONAL_MATRIX<float,3> mSigma;
            for(int k=0;k<3;k++) mSigma(k+1) = Sigma[k][i];
            DIAGONAL_MATRIX<float,3> mQ_hat;
            for(int k=0;k<3;k++) mQ_hat(k+1) = Q_hat[k][i];
            MATRIX_MXN<float> Du(3,8);
            
            const VECTOR<float,3> weights=VECTOR<float,3>::All_Ones_Vector()*(float).5;
            Gradient_Matrix(G,weights,one_over_h[i]);
            
            for(int k=0;k<3;k++) for(int j=0;j<8;j++) Du(k+1,j+1)=u_d[k][j][i];
            MATRIX_MXN<float> mF=Du*G.Transposed();
            for(int k=0;k<9;k++) F.x[k]=mF.x[k];
            F+=1;
            
            F.Fast_Singular_Value_Decomposition(mU,mSigma,mV);
            mSigma=mSigma.Clamp_Min(cutoff[i]);
            
            mQ_hat=TAG_TO_MODEL<T_MATERIAL>::TYPE::Q_hat(mSigma);    
            rdPdF=MATERIAL_MODEL<typename TAG_TO_MODEL<T_MATERIAL>::TYPE>::Isotropic_Stress_Derivative(mSigma,p[i],mu[i],kappa[i],alpha[i],true);

            for(int k=0;k<3;k++)  Q_hat[k][i] = mQ_hat(k+1);
            for(int k=0;k<3;k++) Sigma[k][i] = mSigma(k+1); 
            for(int k=0;k<12;k++)  dPdF[k][i] = *(&(rdPdF.a1111)+k);
            for(int k=0;k<3;k++) for(int j=0;j<3;j++) V[k*3+j][i] = mV(k+1,j+1); 
            for(int k=0;k<3;k++) for(int j=0;j<3;j++) U[k*3+j][i] = mU(k+1,j+1); 
        }
}
};


template<class T>
bool Update_Position_Based_State_Blocked_Compare(const T U[9][8],
                                                 const T V[9][8],
                                                 const T Sigma[3][8],
                                                 const T Q_hat[3][8],
                                                 const T dPdF[12][8],
                                                 const T U_reference[9][8],
                                                 const T V_reference[9][8],
                                                 const T Sigma_reference[3][8],
                                                 const T Q_hat_reference[3][8],
                                                 const T dPdF_reference[12][8])
{


    T Sigma1[8][3];
    T Sigma1_reference[8][3];
    T Q_hat1[8][3];
    T Q_hat1_reference[8][3];
    T dPdF1[8][12];
    T dPdF1_reference[8][12];
            
    for( int j = 0; j < 8; j++)
        {
            for( int k = 0; k < 3; k++)
                {
                    Sigma1[j][k] = Sigma[k][j];
                    Sigma1_reference[j][k] = Sigma_reference[k][j];
                    Q_hat1[j][k] = Q_hat[k][j];
                    Q_hat1_reference[j][k] = Q_hat_reference[k][j];
                }
            for( int k = 0; k < 12; k++)
                {
                    dPdF1[j][k] = dPdF[k][j];
                    dPdF1_reference[j][k] = dPdF_reference[k][j];
                }
        }



    for( int i = 0; i < 8; i++)
        {
            
            
            const DIAGONAL_MATRIX<T,3>& mSigma=*(const DIAGONAL_MATRIX<T,3>*)(Sigma1[i]);
            const DIAGONAL_MATRIX<T,3>& mSigma_reference=*(const DIAGONAL_MATRIX<T,3>*)(Sigma1_reference[i]);
            
            std::cout<<"Computed matrix Sigma :"<<std::endl;Print_Formatted(mSigma,std::cout);
            std::cout<<"Reference matrix Sigma :"<<std::endl;Print_Formatted(mSigma_reference,std::cout);
            std::cout<<"Difference = "<<(mSigma-mSigma_reference).Frobenius_Norm()<<std::endl;
            
            const DIAGONAL_MATRIX<T,3>& mQ_hat=*(const DIAGONAL_MATRIX<T,3>*)(Q_hat1[i]);
            const DIAGONAL_MATRIX<T,3>& mQ_hat_reference=*(const DIAGONAL_MATRIX<T,3>*)(Q_hat1_reference[i]);
            
            std::cout<<"Computed matrix Q_hat :"<<std::endl;Print_Formatted(mQ_hat,std::cout);
            std::cout<<"Reference matrix Q_hat :"<<std::endl;Print_Formatted(mQ_hat_reference,std::cout);
            std::cout<<"Difference = "<<(mQ_hat-mQ_hat_reference).Frobenius_Norm()<<std::endl;
            
            ARRAY_VIEW<const T> adPdF(12,dPdF1[i]);
            ARRAY_VIEW<const T> adPdF_reference(12,dPdF1_reference[i]);
            
            std::cout<<"Computed dPdF : "<<adPdF<<std::endl;
            std::cout<<"Reference dPdF :"<<adPdF_reference<<std::endl;
            ARRAY<T> difference(adPdF-adPdF_reference);
            std::cout<<"Difference = "<<ARRAYS_COMPUTATIONS::Maximum_Magnitude(difference)<<std::endl;
            
            // if( !((mSigma-mSigma_reference).Frobenius_Norm() < 0.00001) )
            //     return false;
            
            // if( !((mQ_hat-mQ_hat_reference).Frobenius_Norm() < 0.00001) )
            //     return false;
            
            // if( !(ARRAYS_COMPUTATIONS::Maximum_Magnitude(difference) < 0.00001) )
            //     return false;
            


        }

    return true;

}

template void Update_Position_Based_State_Blocked_Reference<NEOHOOKEAN_TAG,float>::Run(const float u[3][27], 
    const float p[8], const float mu[8], const float kappa[8], const float alpha[8], const float cutoff[8],
    const float one_over_h[8],float U[9][8], float V[9][8], float Sigma[3][8], float Q_hat[3][8], float dPdF[12][8]);

template bool Update_Position_Based_State_Blocked_Compare(const float U[9][8],
                                                          const float V[9][8],
                                                          const float Sigma[3][8],
                                                          const float Q_hat[3][8],
                                                          const float dPdF[12][8],
                                                          const float U_reference[9][8],
                                                          const float V_reference[9][8],
                                                          const float Sigma_reference[3][8],
                                                          const float Q_hat_reference[3][8],
                                                          const float dPdF_reference[12][8]);

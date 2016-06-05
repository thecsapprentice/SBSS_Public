//#####################################################################
//  Copyright (c) 2011-2013 Nathan Mitchell, Eftychios Sifakis.
//  This file is covered by the FreeBSD license. Please refer to the 
//  license.txt file for more information.
//#####################################################################


struct COROTATED_TAG;
struct NEOHOOKEAN_TAG;

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Arrays_Computations/MAGNITUDE.h>
#include <PhysBAM_Tools/Math_Tools/sqr.h>
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
void Add_Force_Differential_Reference(const T du[3][8], 
                                      const T dp,
                                      const T alpha_squared_over_kappa,
                                      const T alpha,
                                      const T one_over_h,
                                      const T cell_volume,
                                      const T Q_hat[3],
                                      const T U[9],
                                      const T V[9],
                                      const T dPdF[12],
                                      
                                      T df[3][8],
                                      T &dq)
{   
    MATRIX_MXN<float> G;
    MATRIX<float,3> dF;
    const MATRIX<float,3>& mU=*(const MATRIX<float,3>*)(U);
    const MATRIX<float,3>& mV=*(const MATRIX<float,3>*)(V);
    const ROTATED_STRESS_DERIVATIVE<float,3>& rdPdF=*(const ROTATED_STRESS_DERIVATIVE<float,3>*)(dPdF);
    const DIAGONAL_MATRIX<float,3>& mQ_hat=*(const DIAGONAL_MATRIX<float,3>*)(Q_hat);
    MATRIX_MXN<float> Du(3,8);
    MATRIX_MXN<float> Df(3,8);

    const VECTOR<float,3> weights=VECTOR<float,3>::All_Ones_Vector()*(float).5;
    Gradient_Matrix(G,weights,one_over_h);



    for(int i=0;i<3;i++) for(int j=0;j<8;j++) Du(i+1,j+1)=du[i][j];
    for(int i=0;i<3;i++) for(int j=0;j<8;j++) Df(i+1,j+1)=df[i][j];

    MATRIX_MXN<float> mF=Du*G.Transposed();
    for(int i=0;i<9;i++) dF.x[i]=mF.x[i];

    MATRIX<T,3> dF_hat=mU.Transpose_Times(dF)*mV;
    MATRIX<T,3> dP_hat=MATERIAL_MODEL<COROTATED<T,3> >::dP_hat(rdPdF, mQ_hat,dF_hat,dp,alpha);
    MATRIX<T,3> dP=mU*dP_hat.Times_Transpose(mV);

    dq += MATERIAL_MODEL<COROTATED<T,3> >::dq(mQ_hat,dF_hat,dp,alpha,alpha_squared_over_kappa)*cell_volume;

    MATRIX_MXN<T> H=-dP*G*cell_volume;

    Df+=H;

    for(int i=0;i<3;i++) for(int j=0;j<8;j++) df[i][j]=Df(i+1,j+1);
}


template<class T>
bool Add_Force_Differential_Compare( const T df[3][8], const T dq, const T df_reference[3][8], const T dq_reference)
{

    MATRIX_MXN<T> Df(3,8);
    for(int i=0;i<3;i++) for(int j=0;j<8;j++) Df(i+1,j+1)=df[i][j];
    MATRIX_MXN<T> Df_reference(3,8);
    for(int i=0;i<3;i++) for(int j=0;j<8;j++) Df_reference(i+1,j+1)=df_reference[i][j];

    std::cout<<"Computed matrix Df :"<<std::endl;Print_Formatted(Df,std::cout);
    std::cout<<"Reference matrix Df :"<<std::endl;Print_Formatted(Df_reference,std::cout);
    std::cout<<"Difference = "<<(Df-Df_reference).Frobenius_Norm()<<std::endl;

    std::cout<<"Computed dq :"<<std::endl << dq << std::endl;
    std::cout<<"Reference dq :"<<std::endl << dq_reference << std::endl;
    std::cout<<"Difference = "<<(dq-dq_reference)<<std::endl;


    if( !(Df-Df_reference).Frobenius_Norm() < 0.00001 )
        return false;



    if( !(dq-dq_reference) < 0.00001 )
        return false;

    
    return true;

}


template
void Add_Force_Differential_Reference(const float du[3][8], 
                                      const float dp,
                                      const float alpha_squared_over_kappa,
                                      const float alpha,
                                      const float one_over_h,
                                      const float cell_volume,
                                      const float Q_hat[3],
                                      const float U[9],
                                      const float V[9],
                                      const float dPdF[12],
                                      
                                      float df[3][8],
                                      float &dq);


template
bool Add_Force_Differential_Compare( const float df[3][8],
                                     const float dq,
                                     const float df_reference[3][8],
                                     const float dq_reference);


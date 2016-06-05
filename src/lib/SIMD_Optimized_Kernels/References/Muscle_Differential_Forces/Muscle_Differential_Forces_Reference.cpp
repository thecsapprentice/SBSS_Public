//#####################################################################
//  Copyright (c) 2011-2013 Nathan Mitchell, Eftychios Sifakis.
//  This file is covered by the FreeBSD license. Please refer to the 
//  license.txt file for more information.
//#####################################################################


#include <PhysBAM_Tools/Matrices/MATRIX_3X3.h>
#include <PhysBAM_Tools/Math_Tools/RANGE.h>
#include <PhysBAM_Tools/Matrices/MATRIX_MXN.h>
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



template<class T>
void Muscle_Differential_Forces_Reference(T df[3][8], const T du[3][8], const T fiber[3],  const T Ffiber[3], const T c1, const T c2, const T one_over_h, const T cell_volume)
{

    MATRIX_MXN<T> G;
    const VECTOR<T,3> weights(.5,.5,.5);
    Gradient_Matrix(G,weights,one_over_h);
    MATRIX<T,3> mdP_fiber;
    const VECTOR<T,3>& mfiber=*(const VECTOR<T,3>*)(fiber);
    const VECTOR<T,3>& mFfiber=*(const VECTOR<T,3>*)(Ffiber);
    const T mc1 = c1;
    const T mc2 = c2;


    // Gradient Part
    MATRIX<T,3> dF;
    MATRIX_MXN<T> Du(3,8);
    for(int i=0;i<3;i++) for(int j=0;j<8;j++) Du(i+1,j+1)=du[i][j];
    MATRIX_MXN<T> dFp=Du*G.Transposed();
    for(int i=0;i<9;i++) dF.x[i]=dFp.x[i];
 
    // Muscle Differential Part
    MATRIX<T,3> M = MATRIX<T,3>::Outer_Product(mFfiber,mfiber);
    mdP_fiber=mc2*MATRIX<T,3>::Inner_Product(M,dF)*M+mc1*MATRIX<T,3>::Outer_Product(dF*mfiber,mfiber);
    
    //VECTOR<T,3> dw = dF*mfiber;
    //VECTOR<T,3> q = mc1 * dw + mc2 * VECTOR<T,3>::Dot_Product(dw,mFfiber) * mFfiber;
    //mdP_fiber=MATRIX<T,3>::Outer_Product(q,mfiber);   

    // Accumulation Part
    MATRIX_MXN<T> Df(3,8);
    for(int i=0;i<3;i++) for(int j=0;j<8;j++) Df(i+1,j+1)=df[i][j];
    Df+=MATRIX_MXN<T>(mdP_fiber*(-cell_volume))*G;
    for(int i=0;i<3;i++) for(int j=0;j<8;j++) df[i][j]=Df(i+1,j+1);
}

template<class T>
bool Muscle_Differential_Forces_Compare(const T df[3][8], const T df_reference[3][8])
{
    MATRIX_MXN<T> mdf(3,8);
    for(int i=0;i<3;i++) for(int j=0;j<8;j++) mdf(i+1,j+1)=df[i][j];
    MATRIX_MXN<T> mdf_reference(3,8);
    for(int i=0;i<3;i++) for(int j=0;j<8;j++) mdf_reference(i+1,j+1)=df_reference[i][j];

    std::cout<<"Computed matrix mdf :"<<std::endl;Print_Formatted(mdf,std::cout);
    std::cout<<"Reference matrix mdf :"<<std::endl;Print_Formatted(mdf_reference,std::cout);
    std::cout<<"Difference = "<<(mdf-mdf_reference).Frobenius_Norm()<<std::endl;

    if( (mdf-mdf_reference).Frobenius_Norm() < 0.00001 )
        return true;
    else
        return false;

}

template void Muscle_Differential_Forces_Reference(float df[3][8], const float du[3][8], const float fiber[3],  const float Ffiber[3], const float c1, const float c2, const float one_over_h, const float cell_volume);
template bool Muscle_Differential_Forces_Compare(const float df[3][8], const float df_reference[3][8]);


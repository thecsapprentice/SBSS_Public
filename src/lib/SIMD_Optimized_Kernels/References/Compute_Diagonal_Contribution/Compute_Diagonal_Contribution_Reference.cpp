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
//#include <PhysBAM_Tools/Log/LOG.h>
#include "RANGE_ITERATOR.h"
#include "ROTATED_STRESS_DERIVATIVE.h"

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

template<class T>
MATRIX<T,3> Build_M(const ROTATED_STRESS_DERIVATIVE<T,3> &Te, const VECTOR<T,3> &We)
{
    MATRIX<T,3> Matrix_M;
    Matrix_M(1,1) = We(1)*We(1)*Te.a1111 + We(2)*We(2)*Te.a1212 + We(3)*We(3)*Te.a1313;
    Matrix_M(2,2) = We(1)*We(1)*Te.a1212 + We(2)*We(2)*Te.a2222 + We(3)*We(3)*Te.a2323;
    Matrix_M(3,3) = We(1)*We(1)*Te.a1313 + We(2)*We(2)*Te.a2323 + We(3)*We(3)*Te.a3333;
    Matrix_M(2,1) = Matrix_M(1,2) = We(1)*We(2)*Te.a1122 + We(2)*We(1)*Te.a1221;
    Matrix_M(3,1) = Matrix_M(1,3) = We(1)*We(3)*Te.a1133 + We(3)*We(1)*Te.a1331;
    Matrix_M(2,3) = Matrix_M(3,2) = We(3)*We(2)*Te.a2233 + We(2)*We(3)*Te.a2332;
    return Matrix_M;
}
}

template<class T>
void Compute_Diagonal_Contribution_Reference(const T one_over_h,
                                             const T mu_stab,
                                             const T cell_volume,
                                             const T U[9],
                                             const T V[9],
                                             const T dPdF[12],
                                             T d[3][8])
{
    MATRIX<float,3> _U;
    MATRIX<float,3> _V;
    ROTATED_STRESS_DERIVATIVE<float,3> _dPdF;
    MATRIX_MXN<float> diag(3,8);

    
    for( int x=0; x<9; x++){
        _U.x[x] = U[x];
        _V.x[x] = V[x];
    }
    _dPdF.a1111 = dPdF[0];_dPdF.a2233 = dPdF[4];_dPdF.a1313 = dPdF[8];
    _dPdF.a1122 = dPdF[1];_dPdF.a3333 = dPdF[5];_dPdF.a1331 = dPdF[9];
    _dPdF.a1133 = dPdF[2];_dPdF.a1212 = dPdF[6];_dPdF.a2323 = dPdF[10];
    _dPdF.a2222 = dPdF[3];_dPdF.a1221 = dPdF[7];_dPdF.a2332 = dPdF[11];

    for( int x=0; x<3; x++)
        for( int y=0; y<8; y++)
            diag(x+1, y+1) = d[x][y];
    // -------------------------------------------
    
   
    MATRIX_MXN<float> G(3,8); 
    Gradient_Matrix(G, VECTOR<T,3>(.5,.5,.5), one_over_h);
    MATRIX<float,3> Vt(_V.Transposed());
    MATRIX_MXN<float> H = Vt*G;
    
    for(int k=1;k<=3;k++){
        VECTOR<float,3> u(_U(k,1),_U(k,2),_U(k,3));
        MATRIX<float,3> M=Build_M(_dPdF,u);
        for(int l=1;l<=8;l++){
            VECTOR<float,3> h(H(1,l),H(2,l),H(3,l));
            float result=VECTOR<float,3>::Dot_Product(h,M*h)*-cell_volume;                
            diag(k,l)=result-7.f*mu_stab;
        }
    }
    
    // -------------------------------------------
    for( int x=0; x<3; x++)
        for( int y=0; y<8; y++)
            d[x][y] = diag(x+1, y+1);
    
}

template<class T>
bool Compute_Diagonal_Contribution_Compare(const T d[3][8], const T d_reference[3][8])
{
    MATRIX_MXN<T> Dd(3,8);
    for(int i=0;i<3;i++) for(int j=0;j<8;j++) Dd(i+1,j+1)=d[i][j];
    MATRIX_MXN<T> Dd_reference(3,8);
    for(int i=0;i<3;i++) for(int j=0;j<8;j++) Dd_reference(i+1,j+1)=d_reference[i][j];

    std::cout<<"Computed matrix Dd :"<<std::endl;Print_Formatted(Dd,std::cout);
    std::cout<<"Reference matrix Dd :"<<std::endl;Print_Formatted(Dd_reference,std::cout);
    std::cout<<"Difference = "<<(Dd-Dd_reference).Frobenius_Norm()<<std::endl;

    if( (Dd-Dd_reference).Frobenius_Norm() < 0.00001 )
        return true;
    else
        return false;
}

template void  Compute_Diagonal_Contribution_Reference(const float one_over_h,
                                                       const float mu_stab,
                                                       const float cell_volume,
                                                       const float U[9],
                                                       const float V[9],
                                                       const float dPdF[12],
                                                       float d[3][8]);
template bool  Compute_Diagonal_Contribution_Compare(const float d[3][8], const float d_reference[3][8]);
 

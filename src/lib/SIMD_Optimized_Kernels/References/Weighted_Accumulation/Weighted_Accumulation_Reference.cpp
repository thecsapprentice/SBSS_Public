//#####################################################################
//  Copyright (c) 2011-2013 Nathan Mitchell, Eftychios Sifakis.
//  This file is covered by the FreeBSD license. Please refer to the 
//  license.txt file for more information.
//#####################################################################


#include <PhysBAM_Tools/Math_Tools/RANGE.h>
#include <PhysBAM_Tools/Matrices/MATRIX_3X3.h>
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
void Weighted_Accumulation_Reference(T u[3][8], const T F[9],const T W[3], const T one_over_h,const T scale)
{
    MATRIX_MXN<T> G;
    const VECTOR<T,3>& weights=*(const VECTOR<T,3>*)(W);
    Gradient_Matrix(G,weights,one_over_h);
    MATRIX_MXN<T> Du(3,8);
    for(int i=0;i<3;i++) for(int j=0;j<8;j++) Du(i+1,j+1)=u[i][j];
    const MATRIX<T,3>& mF=*(const MATRIX<T,3>*)(F);
    Du+=MATRIX_MXN<T>(mF*scale)*G;
    for(int i=0;i<3;i++) for(int j=0;j<8;j++) u[i][j]=Du(i+1,j+1);
}

template<class T>
bool Weighted_Accumulation_Compare(const T u[3][8], const T u_reference[3][8])
{
    MATRIX_MXN<T> Du(3,8);
    for(int i=0;i<3;i++) for(int j=0;j<8;j++) Du(i+1,j+1)=u[i][j];
    MATRIX_MXN<T> Du_reference(3,8);
    for(int i=0;i<3;i++) for(int j=0;j<8;j++) Du_reference(i+1,j+1)=u_reference[i][j];

    std::cout<<"Computed matrix Du :"<<std::endl;Print_Formatted(Du,std::cout);
    std::cout<<"Reference matrix Du :"<<std::endl;Print_Formatted(Du_reference,std::cout);
    std::cout<<"Difference = "<<(Du-Du_reference).Frobenius_Norm()<<std::endl;

    if( (Du-Du_reference).Frobenius_Norm() < 0.00001 )
        return true;
    else
        return false;
}

template void Weighted_Accumulation_Reference(float u[3][8], const float F[9], const float W[3], const float one_over_h, const float scale);
template bool Weighted_Accumulation_Compare(const float u[3][8], const float u_reference[3][8]);
 

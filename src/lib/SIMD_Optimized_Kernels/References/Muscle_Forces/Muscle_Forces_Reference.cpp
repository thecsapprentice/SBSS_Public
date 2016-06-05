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
void Muscle_Forces_Reference(T f[3][8], const T fiber[3],const T Ffiber[3], const T c1, 
                             const T one_over_h, const T cell_volume)
{
    MATRIX_MXN<T> G;
    const VECTOR<T,3> weights(.5,.5,.5);
    Gradient_Matrix(G,weights,one_over_h);
    MATRIX<T,3> mP_fiber;
    const VECTOR<T,3>& mfiber=*(const VECTOR<T,3>*)(fiber);
    const VECTOR<T,3>& mFfiber=*(const VECTOR<T,3>*)(Ffiber);
    const T mc1 = c1;
    
    mP_fiber = MATRIX<T,3>::Outer_Product(mFfiber, mfiber) * mc1;
    
    // Accumulation Part
    MATRIX_MXN<T> Df(3,8);
    for(int i=0;i<3;i++) for(int j=0;j<8;j++) Df(i+1,j+1)=f[i][j];
    Df+=MATRIX_MXN<T>(mP_fiber*(-cell_volume))*G;
    for(int i=0;i<3;i++) for(int j=0;j<8;j++) f[i][j]=Df(i+1,j+1);
}

template<class T>
bool Muscle_Forces_Compare(const T f[3][8], const T f_reference[3][8])
{
    MATRIX_MXN<T> mf(3,8);
    for(int i=0;i<3;i++) for(int j=0;j<8;j++) mf(i+1,j+1)=f[i][j];
    MATRIX_MXN<T> mf_reference(3,8);
    for(int i=0;i<3;i++) for(int j=0;j<8;j++) mf_reference(i+1,j+1)=f_reference[i][j];

    std::cout<<"Computed matrix mf :"<<std::endl;Print_Formatted(mf,std::cout);
    std::cout<<"Reference matrix mf :"<<std::endl;Print_Formatted(mf_reference,std::cout);
    std::cout<<"Difference = "<<(mf-mf_reference).Frobenius_Norm()<<std::endl;

    if( (mf-mf_reference).Frobenius_Norm() < 0.00001 )
        return true;
    else
        return false;
}

template void Muscle_Forces_Reference(float f[3][8], const float fiber[3],  const float Ffiber[3], const float c1, const float one_over_h, const float cell_volume);
template bool Muscle_Forces_Compare(const float f[3][8], const float f_reference[3][8]);

//#####################################################################
//  Copyright (c) 2011-2013 Nathan Mitchell, Eftychios Sifakis.
//  This file is covered by the FreeBSD license. Please refer to the 
//  license.txt file for more information.
//#####################################################################


#include <PhysBAM_Tools/Matrices/MATRIX_3X3.h>

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
void Muscle_Differential_Reference(T dP_fiber[9], const T dF[9], const T fiber[3],  const T Ffiber[3], const T c1, const T c2 )
{
    MATRIX<T,3>& mdP_fiber=*(MATRIX<T,3>*)(dP_fiber);
    const MATRIX<T,3>& mdF=*(const MATRIX<T,3>*)(dF);
    const VECTOR<T,3>& mfiber=*(const VECTOR<T,3>*)(fiber);
    const VECTOR<T,3>& mFfiber=*(const VECTOR<T,3>*)(Ffiber);
    const T mc1 = c1;
    const T mc2 = c2;
    
    MATRIX<T,3> M = MATRIX<T,3>::Outer_Product(mFfiber,mfiber);
    mdP_fiber=mc2*MATRIX<T,3>::Inner_Product(M,mdF)*M+mc1*MATRIX<T,3>::Outer_Product(mdF*mfiber,mfiber);

    //VECTOR<T,3> dw = mdF*mfiber;
    //VECTOR<T,3> q = mc1 * dw + mc2 * VECTOR<T,3>::Dot_Product(dw,mFfiber) * mFfiber;
    //mdP_fiber=MATRIX<T,3>::Outer_Product(q,mfiber);   

}

template<class T>
bool Muscle_Differential_Compare(const T dP_fiber[9], const T dP_fiber_reference[9])
{
    const MATRIX<T,3>& mdP_fiber=*(const MATRIX<T,3>*)(dP_fiber);
    const MATRIX<T,3>& mdP_fiber_reference=*(const MATRIX<T,3>*)(dP_fiber_reference);

    std::cout<<"Computed matrix dP_fiber :"<<std::endl;Print_Formatted(mdP_fiber,std::cout);
    std::cout<<"Reference matrix dP_fiber :"<<std::endl;Print_Formatted(mdP_fiber_reference,std::cout);
    std::cout<<"Difference = "<<(mdP_fiber-mdP_fiber_reference).Frobenius_Norm()<<std::endl;

    if( (mdP_fiber-mdP_fiber_reference).Frobenius_Norm() < 0.00001 )
        return true;
    else
        return false;
}

template void Muscle_Differential_Reference(float dP_fiber[9], const float dF[9], const float fiber[3],  const float Ffiber[3], const float c1, const float c2);
template bool Muscle_Differential_Compare(const float dP_fiber[9], const float dP_fiber_reference[9]);
 

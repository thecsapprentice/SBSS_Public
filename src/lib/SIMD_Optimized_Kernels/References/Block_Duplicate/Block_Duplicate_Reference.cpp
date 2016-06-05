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
}

template <class T>
void Block_Duplicate_Reference(const T u_compact[3][27], T u_duplicated[3][8][8])
{
    typedef T T3333[3][3][3][3];

    const T3333& u_compact_3333 = *(reinterpret_cast<const T3333*>(u_compact));

    for( int v=0; v<3; v++)
        for( int i_base=0; i_base<2; i_base++)
            for( int j_base=0; j_base<2; j_base++)
                for( int k_base=0; k_base<2; k_base++)
                    for( int i=0; i<2; i++)
                        for( int j=0; j<2; j++)
                            for( int k=0; k<2; k++)
                                {
                                    u_duplicated[v][i*4+j*2+k][i_base*4+j_base*2+k_base] =
                                        u_compact_3333[v][i+i_base][j+j_base][k+k_base];
                                }
}

template <class T>
bool Block_Duplicate_Compare(const T u_duplicated[3][8][8], const T u_duplicated_reference[3][8][8])
{
    
    MATRIX_MXN<T> Du(3,8);
    MATRIX_MXN<T> Du_reference(3,8);

    bool state = true;

    for( int base = 0; base < 8; base ++)
        {
            for(int i=0;i<3;i++) for(int j=0;j<8;j++) Du(i+1,j+1)=u_duplicated[i][j][base];
            for(int i=0;i<3;i++) for(int j=0;j<8;j++) Du_reference(i+1,j+1)=u_duplicated_reference[i][j][base];
            
            std::cout<<"Computed matrix u_duplicated :"<<std::endl;Print_Formatted(Du,std::cout);
            std::cout<<"Reference matrix u_duplicated :"<<std::endl;Print_Formatted(Du_reference,std::cout);
            std::cout<<"Difference = "<<(Du-Du_reference).Frobenius_Norm()<<std::endl;
            
            if( !(Du-Du_reference).Frobenius_Norm() < 0.00001 )
                state = state && false;
        }

    return state;

}



template
void Block_Duplicate_Reference(const float u_compact[3][27], float u_duplicated[3][8][8]);
template
bool Block_Duplicate_Compare(const float u_duplicated[3][8][8], const float u_duplicated_reference[3][8][8]);

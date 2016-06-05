//#####################################################################
//  Copyright (c) 2011-2013 Nathan Mitchell, Eftychios Sifakis.
//  This file is covered by the FreeBSD license. Please refer to the 
//  license.txt file for more information.
//#####################################################################


#include <cmath>
#include <iostream>

template<class T>
void Constraint_Differential_Forces_Reference(T df[3][8], const T du[3][8],
    const T W[3], const T scale )
{

    T w[8];
    
    w[0]=((T)1-W[0])*((T)1-W[1])*((T)1-W[2]);
    w[1]=((T)1-W[0])*((T)1-W[1])*(     W[2]);
    w[2]=((T)1-W[0])*(     W[1])*((T)1-W[2]);
    w[3]=((T)1-W[0])*(     W[1])*(     W[2]);
    w[4]=(     W[0])*((T)1-W[1])*((T)1-W[2]);
    w[5]=(     W[0])*((T)1-W[1])*(     W[2]);
    w[6]=(     W[0])*(     W[1])*((T)1-W[2]);
    w[7]=(     W[0])*(     W[1])*(     W[2]);

    for(int v=0;v<3;v++)
    for(int i=0;i<8;i++)
    for(int j=0;j<8;j++)
        df[v][i]-=scale*w[i]*w[j]*du[v][j];

}

template<class T>
bool Constraint_Differential_Forces_Compare(const T df[3][8], const T df_reference[3][8])
{

    bool is_same = true;

    for(int v=0;v<3;v++)
        for(int i=0;i<8;i++){
            is_same = is_same && (df[v][i] == df_reference[v][i]);
            if( fabs(df[v][i] - df_reference[v][i]) > fabs(1e-6*df_reference[v][i]) )
                {
                    std::cout << "Difference detected at df[" << v << "]["<<i<<"]"<< std::endl;
                    std::cout << "Kernel    : "<< df[v][i] << std::endl;
                    std::cout << "Reference : "<< df_reference[v][i] << std::endl;
                }
        }

    return is_same;

}

template
void Constraint_Differential_Forces_Reference(float df[3][8], const float du[3][8],
    const float W[3], const float scale );

template
bool Constraint_Differential_Forces_Compare(const float df[3][8], const float df_reference[3][8]);


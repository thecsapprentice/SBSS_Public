//#####################################################################
//  Copyright (c) 2011-2013 Nathan Mitchell, Eftychios Sifakis.
//  This file is covered by the FreeBSD license. Please refer to the 
//  license.txt file for more information.
//#####################################################################


#include <cmath>
#include <iostream>

template<class T>
void Collision_Force_Differential_Reference(T df[3][8], const T du[3][8],
                                            const T W[3],  const int spring_id, 
                                            const float* extern_collision_stiffness)
{

    T w[8];
    T force[3];
    T spring_anchor[3];    
    T spring_constant = extern_collision_stiffness[spring_id];
    T spring_attach[3];

    w[0]=((T)1-W[0])*((T)1-W[1])*((T)1-W[2]);
    w[1]=((T)1-W[0])*((T)1-W[1])*(     W[2]);
    w[2]=((T)1-W[0])*(     W[1])*((T)1-W[2]);
    w[3]=((T)1-W[0])*(     W[1])*(     W[2]);
    w[4]=(     W[0])*((T)1-W[1])*((T)1-W[2]);
    w[5]=(     W[0])*((T)1-W[1])*(     W[2]);
    w[6]=(     W[0])*(     W[1])*((T)1-W[2]);
    w[7]=(     W[0])*(     W[1])*(     W[2]);

    for( int i=0; i<3; i++){
        spring_anchor[i] = 0;
        for( int j=0; j<8; j++)
            spring_anchor[i] += w[j] * du[i][j];
        force[i] = spring_constant * (spring_anchor[i]);
    }

    for(int v=0;v<3;v++)
    for(int i=0;i<8;i++)
        df[v][i]-= w[i]*force[v];
}

template<class T>
bool Collision_Force_Differential_Compare(const T df[3][8], const T df_reference[3][8])
{

    bool is_same = true;

    for(int v=0;v<3;v++)
        for(int i=0;i<8;i++){            
            if( fabs(df[v][i] - df_reference[v][i]) > fabs(1e-6*df_reference[v][i]) )
                {
                    is_same == false;
                    std::cout << "Difference detected at df[" << v << "]["<<i<<"]"<< std::endl;
                    std::cout << "Kernel    : "<< df[v][i] << std::endl;
                    std::cout << "Reference : "<< df_reference[v][i] << std::endl;
                }
        }

    return is_same;

}

template
void Collision_Force_Differential_Reference(float df[3][8], const float du[3][8],
                                const float W[3], const int spring_id,
                                const float* extern_collision_stiffness);

template
bool Collision_Force_Differential_Compare(const float df[3][8], const float df_reference[3][8]);


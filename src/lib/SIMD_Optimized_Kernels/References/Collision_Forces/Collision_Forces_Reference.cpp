//#####################################################################
//  Copyright (c) 2011-2013 Nathan Mitchell, Eftychios Sifakis.
//  This file is covered by the FreeBSD license. Please refer to the 
//  license.txt file for more information.
//#####################################################################


#include <cmath>
#include <iostream>

template<class T>
void Collision_Forces_Reference(T df[3][8], const T du[3][8],
                                const T N[3], const T W[3], const T h,
                                const int spring_id, const int spring_id_X,
                                const int spring_id_Y, const int spring_id_Z,
                                const float* extern_collision_attach,
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

    spring_attach[0] = extern_collision_attach[ spring_id_X ];
    spring_attach[1] = extern_collision_attach[ spring_id_Y ];
    spring_attach[2] = extern_collision_attach[ spring_id_Z ];

    for( int i=0; i<3; i++){
        spring_anchor[i] = 0;
        for( int j=0; j<8; j++)
            spring_anchor[i] += w[j] * du[i][j];

        spring_anchor[i] += N[i] + h*W[i];
        force[i] = spring_constant * (spring_attach[i] - spring_anchor[i]);
    }

    for(int v=0;v<3;v++)
    for(int i=0;i<8;i++)
        df[v][i]-= w[i]*force[v];
}

template<class T>
bool Collision_Forces_Compare(const T df[3][8], const T df_reference[3][8])
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
void Collision_Forces_Reference(float df[3][8], const float du[3][8],
                                const float N[3], const float W[3], const float h,
                                const int spring_id, const int spring_id_X,
                                const int spring_id_Y, const int spring_id_Z,
                                const float* extern_collision_attach,
                                const float* extern_collision_stiffness);

template
bool Collision_Forces_Compare(const float df[3][8], const float df_reference[3][8]);


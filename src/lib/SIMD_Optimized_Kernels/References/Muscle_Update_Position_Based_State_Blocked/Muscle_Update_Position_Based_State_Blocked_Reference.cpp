//#####################################################################
//  Copyright (c) 2011-2013 Nathan Mitchell, Eftychios Sifakis.
//  This file is covered by the FreeBSD license. Please refer to the 
//  license.txt file for more information.
//#####################################################################


struct COROTATED_TAG;
struct NEOHOOKEAN_TAG;

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Arrays_Computations/MAGNITUDE.h>
#include <PhysBAM_Tools/Matrices/MATRIX_3X3.h>
#include <PhysBAM_Tools/Matrices/SYMMETRIC_MATRIX_3X3.h>
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


template<class T> T Tension(const T stretch, const T activation, const T density, const T fiber_max_stress)
{
    T fiber_p1=(T).05;
    T fiber_p2=(T)6.6;
    T fiber_cutoff=(T)1.4;
    T cutoff_scaled=fiber_p2*(fiber_cutoff-1);
    T fiber_p3=fiber_p1*fiber_p2*(exp(cutoff_scaled)-1);
    T fiber_p4=fiber_p1*(exp(cutoff_scaled)*(1-fiber_p2*fiber_cutoff)+fiber_p2-1);

    T strain=stretch-1,strain_abs=abs(strain),active_tension=0,passive_tension=0,scale=(T)25/(T)6;
    if(stretch>fiber_cutoff)passive_tension=fiber_p3*stretch+fiber_p4;else if(stretch>1)passive_tension=fiber_p1*(exp(fiber_p2*strain)-fiber_p2*strain-1);
    if(strain_abs<.4)active_tension=activation*density*(1-scale*sqr(strain));else if(strain_abs<.6)active_tension=2*scale*activation*density*sqr(strain_abs-(T).6);
    return fiber_max_stress*(active_tension+passive_tension);
}      
template<class T> T Tension_Derivative(const T stretch, const T activation, const T density, const T fiber_max_stress)
{
    T fiber_p1=(T).05;
    T fiber_p2=(T)6.6;
    T fiber_cutoff=(T)1.4;
    T cutoff_scaled=fiber_p2*(fiber_cutoff-1);
    T fiber_p3=fiber_p1*fiber_p2*(exp(cutoff_scaled)-1);

    T strain=stretch-1,strain_abs=abs(strain),active_tension_derivative=0,passive_tension_derivative=0,scale=(T)25/(T)6;

    if(stretch>fiber_cutoff)passive_tension_derivative=fiber_p3;
    else if(stretch>1)passive_tension_derivative=fiber_p1*fiber_p2*(exp(fiber_p2*strain)-1);
    
    if(strain_abs<.4)active_tension_derivative=-2*scale*activation*density*strain;
    else if(strain_abs<.6)active_tension_derivative=4*scale*activation*density*(strain-sign(strain)*(T).6);

    return fiber_max_stress*(active_tension_derivative+passive_tension_derivative);
}


template<class T, class I> 
struct Muscle_Update_Position_Based_State_Blocked_Reference
{
static void Run(const float u[3][27], 
                const I muscle_id[8],
                const T fiber[3][8],
                const T density[8],
                const T one_over_h[8],
                
                T c1[8],
                T c2[8],
                T F_fiber[3][8],
                
                const float *activations,
                const float *fiber_max_stresses)
{   



    typedef float T3333[3][3][3][3];

    const T3333& u_compact_3333 = *(reinterpret_cast<const T3333*>(u));
    float u_d[3][8][8];


    for( int v=0; v<3; v++)
        for( int i_base=0; i_base<2; i_base++)
            for( int j_base=0; j_base<2; j_base++)
                for( int k_base=0; k_base<2; k_base++)
                    for( int i=0; i<2; i++)
                        for( int j=0; j<2; j++)
                            for( int k=0; k<2; k++)
                                {
                                    u_d[v][i*4+j*2+k][i_base*4+j_base*2+k_base] =
                                        u_compact_3333[v][i+i_base][j+j_base][k+k_base];
                                }

    for( int i = 0; i < 8; i++)
        {       

            MATRIX_MXN<float> G;
            MATRIX<float,3> F;
            MATRIX_MXN<float> Du(3,8);
  
            const VECTOR<float,3> weights=VECTOR<float,3>::All_Ones_Vector()*(float).5;
            Gradient_Matrix(G,weights,one_over_h[i]);
            
            for(int k=0;k<3;k++) for(int j=0;j<8;j++) Du(k+1,j+1)=u_d[k][j][i];
            MATRIX_MXN<float> mF=Du*G.Transposed();
            for(int k=0;k<9;k++) F.x[k]=mF.x[k];
            F+=1;
            
            const VECTOR<float,3> mfiber(fiber[0][i],fiber[1][i],fiber[2][i]);
            VECTOR<float,3> mF_fiber =  F * mfiber;
            F_fiber[0][i] = mF_fiber(1);
            F_fiber[1][i] = mF_fiber(2);
            F_fiber[2][i] = mF_fiber(3);

            T activation=activations[muscle_id[i]];
            T fiber_max_stress=fiber_max_stresses[muscle_id[i]];

            T stretch_squared=mF_fiber.Magnitude_Squared(),stretch=sqrt(stretch_squared);
            stretch = min( 3.0f, stretch );
            T tension=Tension(stretch,activation,density[i],fiber_max_stress);
            T tension_derivative=max((T)0,Tension_Derivative(stretch,activation,density[i],fiber_max_stress));

            c1[i]=tension/stretch;
            c2[i]=(tension_derivative-c1[i])/stretch_squared;

            if( density[i] == 0 )
                {
                    F_fiber[0][i] = 0;
                    F_fiber[1][i] = 0;
                    F_fiber[2][i] = 0;
                    c1[i]=0;
                    c2[i]=0;
                }


        }
}
};


template<class T,class I>
bool Muscle_Update_Position_Based_State_Blocked_Compare(const T c1[8],
                                                            const T c2[8],
                                                            const T F_fiber[3][8],
                                                            const T c1_reference[8],
                                                            const T c2_reference[8],
                                                            const T F_fiber_reference[3][8])
{
            
    ARRAY_VIEW<const T> ac1(8,c1);
    ARRAY_VIEW<const T> ac1_reference(8,c1_reference);           
    std::cout<<"Computed c1 : "<<ac1<<std::endl;
    std::cout<<"Reference c1 :"<<ac1_reference<<std::endl;
    ARRAY<T> c1d(ac1-ac1_reference);
    std::cout<<"Difference = "<<ARRAYS_COMPUTATIONS::Maximum_Magnitude(c1d)<<std::endl;


    ARRAY_VIEW<const T> ac2(8,c2);
    ARRAY_VIEW<const T> ac2_reference(8,c2_reference);           
    std::cout<<"Computed c2 : "<<ac2<<std::endl;
    std::cout<<"Reference c2 :"<<ac2_reference<<std::endl;
    ARRAY<T> c2d(ac2-ac2_reference);
    std::cout<<"Difference = "<<ARRAYS_COMPUTATIONS::Maximum_Magnitude(c2d)<<std::endl;
            
    
    for( int i=0 ; i <3; i++)
        {
            ARRAY_VIEW<const T> aFfiber(8,F_fiber[i]);
            ARRAY_VIEW<const T> aFfiber_reference(8,F_fiber_reference[i]);           
            std::cout<<"Computed F_Fiber ("<<i<<") : "<<aFfiber<<std::endl;
            std::cout<<"Reference F_Fiber ("<<i<<") :"<<aFfiber_reference<<std::endl;
            ARRAY<T> Ffd(aFfiber-aFfiber_reference);
            std::cout<<"Difference = "<<ARRAYS_COMPUTATIONS::Maximum_Magnitude(Ffd)<<std::endl;
        }

    return true;

}

template void Muscle_Update_Position_Based_State_Blocked_Reference<float,int>::Run(const float u[3][27], 
                                                                                   const int muscle_id[8],
                                                                                   const float fiber[3][8],
                                                                                   const float density[8],
                                                                                   const float one_over_h[8],
                                                                                   float c1[8],
                                                                                   float c2[8],
                                                                                   float F_fiber[3][8],
                                                                                   const float *activations,
                                                                                   const float *fiber_max_stresses);



template bool Muscle_Update_Position_Based_State_Blocked_Compare<float,int>(const float c1[8],
                                                                 const float c2[8],
                                                                 const float F_fiber[3][8],
                                                                 const float c1_reference[8],
                                                                 const float c2_reference[8],
                                                                 const float F_fiber_reference[3][8]
                                                                 );
                                         

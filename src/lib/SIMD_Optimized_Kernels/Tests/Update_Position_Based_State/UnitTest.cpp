
#include <cstdlib>
#include <iostream>

struct NEOHOOKEAN_TAG; struct COROTATED_TAG;

#include "Update_Position_Based_State.h"
#include "Update_Position_Based_State_Reference.h"

template<class T>
T Get_Random(const T a=(T)-1.,const T b=(T)1.)
{
    return ((b-a)*(T)rand())/(T)RAND_MAX+a;
}

int main(int argc,char* argv[])
{
    typedef float T;

    int seed=1;
    if(argc==2) seed=atoi(argv[1]);
    srand(seed);

    

{
T u[3][8] __attribute__((aligned(4)));
T p __attribute__((aligned(4)));
T mu __attribute__((aligned(4)));
T mu_stab __attribute__((aligned(4)));
T kappa __attribute__((aligned(4)));
T alpha __attribute__((aligned(4)));
T cutoff __attribute__((aligned(4)));
T one_over_h __attribute__((aligned(4)));
T cell_volume __attribute__((aligned(4)));
T U[9] __attribute__((aligned(4)));
T U_reference[9] __attribute__((aligned(4)));
T U_original[9] __attribute__((aligned(4)));
T V[9] __attribute__((aligned(4)));
T V_reference[9] __attribute__((aligned(4)));
T V_original[9] __attribute__((aligned(4)));
T Sigma[3] __attribute__((aligned(4)));
T Sigma_reference[3] __attribute__((aligned(4)));
T Sigma_original[3] __attribute__((aligned(4)));
T Q_hat[3] __attribute__((aligned(4)));
T Q_hat_reference[3] __attribute__((aligned(4)));
T Q_hat_original[3] __attribute__((aligned(4)));
T dPdF[12] __attribute__((aligned(4)));
T dPdF_reference[12] __attribute__((aligned(4)));
T dPdF_original[12] __attribute__((aligned(4)));
T d[3][8] __attribute__((aligned(4)));
T d_reference[3][8] __attribute__((aligned(4)));
T d_original[3][8] __attribute__((aligned(4)));
T system_matrix[300] __attribute__((aligned(4)));
T system_matrix_reference[300] __attribute__((aligned(4)));
T system_matrix_original[300] __attribute__((aligned(4)));
  

for(int __a=0;__a<3;__a++) for(int __b=0;__b<8;__b++) 
u[__a][__b]=Get_Random<float>();
p=Get_Random<float>();
mu=Get_Random<float>();
mu_stab=Get_Random<float>();
kappa=Get_Random<float>();
alpha=Get_Random<float>();
cutoff=Get_Random<float>();
one_over_h=Get_Random<float>();
cell_volume=Get_Random<float>();for(int __a=0;__a<9;__a++) {
U_original[__a]=Get_Random<float>();
U[__a]=U_original[__a];
U_reference[__a]=U_original[__a];}
for(int __a=0;__a<9;__a++) {
V_original[__a]=Get_Random<float>();
V[__a]=V_original[__a];
V_reference[__a]=V_original[__a];}
for(int __a=0;__a<3;__a++) {
Sigma_original[__a]=Get_Random<float>();
Sigma[__a]=Sigma_original[__a];
Sigma_reference[__a]=Sigma_original[__a];}
for(int __a=0;__a<3;__a++) {
Q_hat_original[__a]=Get_Random<float>();
Q_hat[__a]=Q_hat_original[__a];
Q_hat_reference[__a]=Q_hat_original[__a];}
for(int __a=0;__a<12;__a++) {
dPdF_original[__a]=Get_Random<float>();
dPdF[__a]=dPdF_original[__a];
dPdF_reference[__a]=dPdF_original[__a];}
for(int __a=0;__a<3;__a++) for(int __b=0;__b<8;__b++) {
d_original[__a][__b]=Get_Random<float>();
d[__a][__b]=d_original[__a][__b];
d_reference[__a][__b]=d_original[__a][__b];}
for(int __a=0;__a<300;__a++) {
system_matrix_original[__a]=Get_Random<float>();
system_matrix[__a]=system_matrix_original[__a];
system_matrix_reference[__a]=system_matrix_original[__a];}


for(int i=0; i<1; i+=1) {Update_Position_Based_State<NEOHOOKEAN_TAG,float,float,int>::Run(u,p,mu,mu_stab,kappa,alpha,cutoff,one_over_h,cell_volume,U,V,Sigma,Q_hat,dPdF,d,system_matrix);}

Update_Position_Based_State_Neohookean_Reference<float>(u,p,mu,mu_stab,kappa,alpha,cutoff,one_over_h,cell_volume,U_reference,V_reference,Sigma_reference,Q_hat_reference,dPdF_reference,d_reference,system_matrix_reference);
if( !(Update_Position_Based_State_Neohookean_Compare<float>(U,V,Sigma,Q_hat,dPdF,d,system_matrix,U_reference,V_reference,Sigma_reference,Q_hat_reference,dPdF_reference,d_reference,system_matrix_reference)) ){
   std::cout << "Failed to confirm unit test for Update_Position_Based_State with material NEOHOOKEAN" << std::endl;
   return 1;
}

}



{
T u[3][8] __attribute__((aligned(4)));
T p __attribute__((aligned(4)));
T mu __attribute__((aligned(4)));
T mu_stab __attribute__((aligned(4)));
T kappa __attribute__((aligned(4)));
T alpha __attribute__((aligned(4)));
T cutoff __attribute__((aligned(4)));
T one_over_h __attribute__((aligned(4)));
T cell_volume __attribute__((aligned(4)));
T U[9] __attribute__((aligned(4)));
T U_reference[9] __attribute__((aligned(4)));
T U_original[9] __attribute__((aligned(4)));
T V[9] __attribute__((aligned(4)));
T V_reference[9] __attribute__((aligned(4)));
T V_original[9] __attribute__((aligned(4)));
T Sigma[3] __attribute__((aligned(4)));
T Sigma_reference[3] __attribute__((aligned(4)));
T Sigma_original[3] __attribute__((aligned(4)));
T Q_hat[3] __attribute__((aligned(4)));
T Q_hat_reference[3] __attribute__((aligned(4)));
T Q_hat_original[3] __attribute__((aligned(4)));
T dPdF[12] __attribute__((aligned(4)));
T dPdF_reference[12] __attribute__((aligned(4)));
T dPdF_original[12] __attribute__((aligned(4)));
T d[3][8] __attribute__((aligned(4)));
T d_reference[3][8] __attribute__((aligned(4)));
T d_original[3][8] __attribute__((aligned(4)));
T system_matrix[300] __attribute__((aligned(4)));
T system_matrix_reference[300] __attribute__((aligned(4)));
T system_matrix_original[300] __attribute__((aligned(4)));
  

for(int __a=0;__a<3;__a++) for(int __b=0;__b<8;__b++) 
u[__a][__b]=Get_Random<float>();
p=Get_Random<float>();
mu=Get_Random<float>();
mu_stab=Get_Random<float>();
kappa=Get_Random<float>();
alpha=Get_Random<float>();
cutoff=Get_Random<float>();
one_over_h=Get_Random<float>();
cell_volume=Get_Random<float>();for(int __a=0;__a<9;__a++) {
U_original[__a]=Get_Random<float>();
U[__a]=U_original[__a];
U_reference[__a]=U_original[__a];}
for(int __a=0;__a<9;__a++) {
V_original[__a]=Get_Random<float>();
V[__a]=V_original[__a];
V_reference[__a]=V_original[__a];}
for(int __a=0;__a<3;__a++) {
Sigma_original[__a]=Get_Random<float>();
Sigma[__a]=Sigma_original[__a];
Sigma_reference[__a]=Sigma_original[__a];}
for(int __a=0;__a<3;__a++) {
Q_hat_original[__a]=Get_Random<float>();
Q_hat[__a]=Q_hat_original[__a];
Q_hat_reference[__a]=Q_hat_original[__a];}
for(int __a=0;__a<12;__a++) {
dPdF_original[__a]=Get_Random<float>();
dPdF[__a]=dPdF_original[__a];
dPdF_reference[__a]=dPdF_original[__a];}
for(int __a=0;__a<3;__a++) for(int __b=0;__b<8;__b++) {
d_original[__a][__b]=Get_Random<float>();
d[__a][__b]=d_original[__a][__b];
d_reference[__a][__b]=d_original[__a][__b];}
for(int __a=0;__a<300;__a++) {
system_matrix_original[__a]=Get_Random<float>();
system_matrix[__a]=system_matrix_original[__a];
system_matrix_reference[__a]=system_matrix_original[__a];}


for(int i=0; i<1; i+=1) {Update_Position_Based_State<COROTATED_TAG,float,float,int>::Run(u,p,mu,mu_stab,kappa,alpha,cutoff,one_over_h,cell_volume,U,V,Sigma,Q_hat,dPdF,d,system_matrix);}

Update_Position_Based_State_Corotated_Reference<float>(u,p,mu,mu_stab,kappa,alpha,cutoff,one_over_h,cell_volume,U_reference,V_reference,Sigma_reference,Q_hat_reference,dPdF_reference,d_reference,system_matrix_reference);
if( !(Update_Position_Based_State_Corotated_Compare<float>(U,V,Sigma,Q_hat,dPdF,d,system_matrix,U_reference,V_reference,Sigma_reference,Q_hat_reference,dPdF_reference,d_reference,system_matrix_reference)) ){
   std::cout << "Failed to confirm unit test for Update_Position_Based_State with material COROTATED" << std::endl;
   return 1;
}

}



    return 0;

}



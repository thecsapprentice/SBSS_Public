
#include <cstdlib>
#include <iostream>
#include "KernelCommon.h"

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

    std::cout.precision(10);
    std::cout.setf(std::ios::fixed,std::ios::floatfield);

    

{
std::cout << "Running SIMD Test for Update_Position_Based_State " << std::endl;


//=======================================================
//
//        DEFINE ALL VARIABLES USED BY KERNEL
//
//=======================================================

T u[3][8][16] __attribute__((aligned(64)));
T p[16] __attribute__((aligned(64)));
T mu[16] __attribute__((aligned(64)));
T mu_stab[16] __attribute__((aligned(64)));
T kappa[16] __attribute__((aligned(64)));
T alpha[16] __attribute__((aligned(64)));
T cutoff[16] __attribute__((aligned(64)));
T one_over_h[16] __attribute__((aligned(64)));
T cell_volume[16] __attribute__((aligned(64)));
T U[9][16] __attribute__((aligned(64)));
T U_reference[9][16] __attribute__((aligned(64)));
T U_original[9][16] __attribute__((aligned(64)));
T V[9][16] __attribute__((aligned(64)));
T V_reference[9][16] __attribute__((aligned(64)));
T V_original[9][16] __attribute__((aligned(64)));
T Sigma[3][16] __attribute__((aligned(64)));
T Sigma_reference[3][16] __attribute__((aligned(64)));
T Sigma_original[3][16] __attribute__((aligned(64)));
T Q_hat[3][16] __attribute__((aligned(64)));
T Q_hat_reference[3][16] __attribute__((aligned(64)));
T Q_hat_original[3][16] __attribute__((aligned(64)));
T dPdF[12][16] __attribute__((aligned(64)));
T dPdF_reference[12][16] __attribute__((aligned(64)));
T dPdF_original[12][16] __attribute__((aligned(64)));
T d[3][8][16] __attribute__((aligned(64)));
T d_reference[3][8][16] __attribute__((aligned(64)));
T d_original[3][8][16] __attribute__((aligned(64)));
T system_matrix[300][16] __attribute__((aligned(64)));
T system_matrix_reference[300][16] __attribute__((aligned(64)));
T system_matrix_original[300][16] __attribute__((aligned(64)));
  

for(int __a=0;__a<3;__a++) for(int __b=0;__b<8;__b++) for(int __c=0;__c<16;__c++) 
u[__a][__b][__c]=Get_Random<float>();for(int __a=0;__a<16;__a++) 
p[__a]=Get_Random<float>();for(int __a=0;__a<16;__a++) 
mu[__a]=Get_Random<float>();for(int __a=0;__a<16;__a++) 
mu_stab[__a]=Get_Random<float>();for(int __a=0;__a<16;__a++) 
kappa[__a]=Get_Random<float>();for(int __a=0;__a<16;__a++) 
alpha[__a]=Get_Random<float>();for(int __a=0;__a<16;__a++) 
cutoff[__a]=Get_Random<float>();for(int __a=0;__a<16;__a++) 
one_over_h[__a]=Get_Random<float>();for(int __a=0;__a<16;__a++) 
cell_volume[__a]=Get_Random<float>();for(int __a=0;__a<9;__a++) for(int __b=0;__b<16;__b++) {
U_original[__a][__b]=Get_Random<float>();
U[__a][__b]=U_original[__a][__b];
U_reference[__a][__b]=U_original[__a][__b];}
for(int __a=0;__a<9;__a++) for(int __b=0;__b<16;__b++) {
V_original[__a][__b]=Get_Random<float>();
V[__a][__b]=V_original[__a][__b];
V_reference[__a][__b]=V_original[__a][__b];}
for(int __a=0;__a<3;__a++) for(int __b=0;__b<16;__b++) {
Sigma_original[__a][__b]=Get_Random<float>();
Sigma[__a][__b]=Sigma_original[__a][__b];
Sigma_reference[__a][__b]=Sigma_original[__a][__b];}
for(int __a=0;__a<3;__a++) for(int __b=0;__b<16;__b++) {
Q_hat_original[__a][__b]=Get_Random<float>();
Q_hat[__a][__b]=Q_hat_original[__a][__b];
Q_hat_reference[__a][__b]=Q_hat_original[__a][__b];}
for(int __a=0;__a<12;__a++) for(int __b=0;__b<16;__b++) {
dPdF_original[__a][__b]=Get_Random<float>();
dPdF[__a][__b]=dPdF_original[__a][__b];
dPdF_reference[__a][__b]=dPdF_original[__a][__b];}
for(int __a=0;__a<3;__a++) for(int __b=0;__b<8;__b++) for(int __c=0;__c<16;__c++) {
d_original[__a][__b][__c]=Get_Random<float>();
d[__a][__b][__c]=d_original[__a][__b][__c];
d_reference[__a][__b][__c]=d_original[__a][__b][__c];}
for(int __a=0;__a<300;__a++) for(int __b=0;__b<16;__b++) {
system_matrix_original[__a][__b]=Get_Random<float>();
system_matrix[__a][__b]=system_matrix_original[__a][__b];
system_matrix_reference[__a][__b]=system_matrix_original[__a][__b];}


//=======================================================
//
//             COMPUTE REFERENCE RESULTS
//
//=======================================================
 
T __mu[3][8] __attribute__((aligned(4)));
T __mp __attribute__((aligned(4)));
T __mmu __attribute__((aligned(4)));
T __mmu_stab __attribute__((aligned(4)));
T __mkappa __attribute__((aligned(4)));
T __malpha __attribute__((aligned(4)));
T __mcutoff __attribute__((aligned(4)));
T __mone_over_h __attribute__((aligned(4)));
T __mcell_volume __attribute__((aligned(4)));
T __mU[9] __attribute__((aligned(4)));
T __mU_reference[9] __attribute__((aligned(4)));
T __mU_original[9] __attribute__((aligned(4)));
T __mV[9] __attribute__((aligned(4)));
T __mV_reference[9] __attribute__((aligned(4)));
T __mV_original[9] __attribute__((aligned(4)));
T __mSigma[3] __attribute__((aligned(4)));
T __mSigma_reference[3] __attribute__((aligned(4)));
T __mSigma_original[3] __attribute__((aligned(4)));
T __mQ_hat[3] __attribute__((aligned(4)));
T __mQ_hat_reference[3] __attribute__((aligned(4)));
T __mQ_hat_original[3] __attribute__((aligned(4)));
T __mdPdF[12] __attribute__((aligned(4)));
T __mdPdF_reference[12] __attribute__((aligned(4)));
T __mdPdF_original[12] __attribute__((aligned(4)));
T __md[3][8] __attribute__((aligned(4)));
T __md_reference[3][8] __attribute__((aligned(4)));
T __md_original[3][8] __attribute__((aligned(4)));
T __msystem_matrix[300] __attribute__((aligned(4)));
T __msystem_matrix_reference[300] __attribute__((aligned(4)));
T __msystem_matrix_original[300] __attribute__((aligned(4)));
for( int k=0;k<16;k++){for(int __a=0;__a<3;__a++) for(int __b=0;__b<8;__b++) 
__mu[__a][__b]=u[__a][__b][k];
__mp=p[k];
__mmu=mu[k];
__mmu_stab=mu_stab[k];
__mkappa=kappa[k];
__malpha=alpha[k];
__mcutoff=cutoff[k];
__mone_over_h=one_over_h[k];
__mcell_volume=cell_volume[k];for(int __a=0;__a<9;__a++) 
__mU_reference[__a]=U_reference[__a][k];
for(int __a=0;__a<9;__a++) 
__mV_reference[__a]=V_reference[__a][k];
for(int __a=0;__a<3;__a++) 
__mSigma_reference[__a]=Sigma_reference[__a][k];
for(int __a=0;__a<3;__a++) 
__mQ_hat_reference[__a]=Q_hat_reference[__a][k];
for(int __a=0;__a<12;__a++) 
__mdPdF_reference[__a]=dPdF_reference[__a][k];
for(int __a=0;__a<3;__a++) for(int __b=0;__b<8;__b++) 
__md_reference[__a][__b]=d_reference[__a][__b][k];
for(int __a=0;__a<300;__a++) 
__msystem_matrix_reference[__a]=system_matrix_reference[__a][k];
Update_Position_Based_State<NEOHOOKEAN_TAG,float,float,int>::Run(__mu,__mp,__mmu,__mmu_stab,__mkappa,__malpha,__mcutoff,__mone_over_h,__mcell_volume,__mU_reference,__mV_reference,__mSigma_reference,__mQ_hat_reference,__mdPdF_reference,__md_reference,__msystem_matrix_reference);for(int __a=0;__a<9;__a++) U_reference[__a][k]=
__mU_reference[__a];
for(int __a=0;__a<9;__a++) V_reference[__a][k]=
__mV_reference[__a];
for(int __a=0;__a<3;__a++) Sigma_reference[__a][k]=
__mSigma_reference[__a];
for(int __a=0;__a<3;__a++) Q_hat_reference[__a][k]=
__mQ_hat_reference[__a];
for(int __a=0;__a<12;__a++) dPdF_reference[__a][k]=
__mdPdF_reference[__a];
for(int __a=0;__a<3;__a++) for(int __b=0;__b<8;__b++) d_reference[__a][__b][k]=
__md_reference[__a][__b];
for(int __a=0;__a<300;__a++) system_matrix_reference[__a][k]=
__msystem_matrix_reference[__a];
}

//=======================================================
//
//               COMPUTE SCALAR RESULTS
//
//=======================================================

{
typedef T (&refArray1)[3][8][16];typedef T (&refArray2)[16];typedef T (&refArray3)[16];typedef T (&refArray4)[16];typedef T (&refArray5)[16];typedef T (&refArray6)[16];typedef T (&refArray7)[16];typedef T (&refArray8)[16];typedef T (&refArray9)[16];typedef T (&refArray10)[9][16];typedef T (&refArray11)[9][16];typedef T (&refArray12)[3][16];typedef T (&refArray13)[3][16];typedef T (&refArray14)[12][16];typedef T (&refArray15)[3][8][16];typedef T (&refArray16)[300][16];for(int __a=0;__a<9;__a++) for(int __b=0;__b<16;__b++) U[__a][__b]=U_original[__a][__b];for(int __a=0;__a<9;__a++) for(int __b=0;__b<16;__b++) V[__a][__b]=V_original[__a][__b];for(int __a=0;__a<3;__a++) for(int __b=0;__b<16;__b++) Sigma[__a][__b]=Sigma_original[__a][__b];for(int __a=0;__a<3;__a++) for(int __b=0;__b<16;__b++) Q_hat[__a][__b]=Q_hat_original[__a][__b];for(int __a=0;__a<12;__a++) for(int __b=0;__b<16;__b++) dPdF[__a][__b]=dPdF_original[__a][__b];for(int __a=0;__a<3;__a++) for(int __b=0;__b<8;__b++) for(int __c=0;__c<16;__c++) d[__a][__b][__c]=d_original[__a][__b][__c];for(int __a=0;__a<300;__a++) for(int __b=0;__b<16;__b++) system_matrix[__a][__b]=system_matrix_original[__a][__b];for(int i=0; i<16; i+=1) {refArray1 uk=reinterpret_cast<refArray1>(u[0][0][i]);refArray2 pk=reinterpret_cast<refArray2>(p[i]);refArray3 muk=reinterpret_cast<refArray3>(mu[i]);refArray4 mu_stabk=reinterpret_cast<refArray4>(mu_stab[i]);refArray5 kappak=reinterpret_cast<refArray5>(kappa[i]);refArray6 alphak=reinterpret_cast<refArray6>(alpha[i]);refArray7 cutoffk=reinterpret_cast<refArray7>(cutoff[i]);refArray8 one_over_hk=reinterpret_cast<refArray8>(one_over_h[i]);refArray9 cell_volumek=reinterpret_cast<refArray9>(cell_volume[i]);refArray10 Uk=reinterpret_cast<refArray10>(U[0][i]);refArray11 Vk=reinterpret_cast<refArray11>(V[0][i]);refArray12 Sigmak=reinterpret_cast<refArray12>(Sigma[0][i]);refArray13 Q_hatk=reinterpret_cast<refArray13>(Q_hat[0][i]);refArray14 dPdFk=reinterpret_cast<refArray14>(dPdF[0][i]);refArray15 dk=reinterpret_cast<refArray15>(d[0][0][i]);refArray16 system_matrixk=reinterpret_cast<refArray16>(system_matrix[0][i]);Update_Position_Based_State<NEOHOOKEAN_TAG,float,float[16],int[16]>::Run(uk,pk,muk,mu_stabk,kappak,alphak,cutoffk,one_over_hk,cell_volumek,Uk,Vk,Sigmak,Q_hatk,dPdFk,dk,system_matrixk);}
for(int __a=0;__a<9;__a++) for(int __b=0;__b<16;__b++) if(std::abs((U[__a][__b] - U_reference[__a][__b]) / (U_reference[__a][__b])) > 1 ){std::cerr<<"Mismatch detected in SCALAR implementation"<<std::endl;std::cerr<<"Variable U:"<<std::endl;std::cerr<<"seed="<<seed<<", __a="<<__a<<", __b="<<__b<<std::endl;std::cerr<<"U SCALAR=  "<<U[__a][__b]<<std::endl;std::cerr<<"U Reference=  "<<U_reference[__a][__b]<<std::endl;std::cerr<<"U Rel Difference=  "<< std::abs((U[__a][__b] - U_reference[__a][__b]) / (U_reference[__a][__b])) << std::endl;std::cerr<<"U Abs Difference=  "<< std::abs(U[__a][__b] - U_reference[__a][__b]) << std::endl;return 1;}
for(int __a=0;__a<9;__a++) for(int __b=0;__b<16;__b++) if(std::abs((V[__a][__b] - V_reference[__a][__b]) / (V_reference[__a][__b])) > 1 ){std::cerr<<"Mismatch detected in SCALAR implementation"<<std::endl;std::cerr<<"Variable V:"<<std::endl;std::cerr<<"seed="<<seed<<", __a="<<__a<<", __b="<<__b<<std::endl;std::cerr<<"V SCALAR=  "<<V[__a][__b]<<std::endl;std::cerr<<"V Reference=  "<<V_reference[__a][__b]<<std::endl;std::cerr<<"V Rel Difference=  "<< std::abs((V[__a][__b] - V_reference[__a][__b]) / (V_reference[__a][__b])) << std::endl;std::cerr<<"V Abs Difference=  "<< std::abs(V[__a][__b] - V_reference[__a][__b]) << std::endl;return 1;}
for(int __a=0;__a<3;__a++) for(int __b=0;__b<16;__b++) if(std::abs((Sigma[__a][__b] - Sigma_reference[__a][__b]) / (Sigma_reference[__a][__b])) > 1 ){std::cerr<<"Mismatch detected in SCALAR implementation"<<std::endl;std::cerr<<"Variable Sigma:"<<std::endl;std::cerr<<"seed="<<seed<<", __a="<<__a<<", __b="<<__b<<std::endl;std::cerr<<"Sigma SCALAR=  "<<Sigma[__a][__b]<<std::endl;std::cerr<<"Sigma Reference=  "<<Sigma_reference[__a][__b]<<std::endl;std::cerr<<"Sigma Rel Difference=  "<< std::abs((Sigma[__a][__b] - Sigma_reference[__a][__b]) / (Sigma_reference[__a][__b])) << std::endl;std::cerr<<"Sigma Abs Difference=  "<< std::abs(Sigma[__a][__b] - Sigma_reference[__a][__b]) << std::endl;return 1;}
for(int __a=0;__a<3;__a++) for(int __b=0;__b<16;__b++) if(std::abs((Q_hat[__a][__b] - Q_hat_reference[__a][__b]) / (Q_hat_reference[__a][__b])) > 1 ){std::cerr<<"Mismatch detected in SCALAR implementation"<<std::endl;std::cerr<<"Variable Q_hat:"<<std::endl;std::cerr<<"seed="<<seed<<", __a="<<__a<<", __b="<<__b<<std::endl;std::cerr<<"Q_hat SCALAR=  "<<Q_hat[__a][__b]<<std::endl;std::cerr<<"Q_hat Reference=  "<<Q_hat_reference[__a][__b]<<std::endl;std::cerr<<"Q_hat Rel Difference=  "<< std::abs((Q_hat[__a][__b] - Q_hat_reference[__a][__b]) / (Q_hat_reference[__a][__b])) << std::endl;std::cerr<<"Q_hat Abs Difference=  "<< std::abs(Q_hat[__a][__b] - Q_hat_reference[__a][__b]) << std::endl;return 1;}
for(int __a=0;__a<12;__a++) for(int __b=0;__b<16;__b++) if(std::abs((dPdF[__a][__b] - dPdF_reference[__a][__b]) / (dPdF_reference[__a][__b])) > 1 ){std::cerr<<"Mismatch detected in SCALAR implementation"<<std::endl;std::cerr<<"Variable dPdF:"<<std::endl;std::cerr<<"seed="<<seed<<", __a="<<__a<<", __b="<<__b<<std::endl;std::cerr<<"dPdF SCALAR=  "<<dPdF[__a][__b]<<std::endl;std::cerr<<"dPdF Reference=  "<<dPdF_reference[__a][__b]<<std::endl;std::cerr<<"dPdF Rel Difference=  "<< std::abs((dPdF[__a][__b] - dPdF_reference[__a][__b]) / (dPdF_reference[__a][__b])) << std::endl;std::cerr<<"dPdF Abs Difference=  "<< std::abs(dPdF[__a][__b] - dPdF_reference[__a][__b]) << std::endl;return 1;}
for(int __a=0;__a<3;__a++) for(int __b=0;__b<8;__b++) for(int __c=0;__c<16;__c++) if(std::abs((d[__a][__b][__c] - d_reference[__a][__b][__c]) / (d_reference[__a][__b][__c])) > 1 ){std::cerr<<"Mismatch detected in SCALAR implementation"<<std::endl;std::cerr<<"Variable d:"<<std::endl;std::cerr<<"seed="<<seed<<", __a="<<__a<<", __b="<<__b<<", __c="<<__c<<std::endl;std::cerr<<"d SCALAR=  "<<d[__a][__b][__c]<<std::endl;std::cerr<<"d Reference=  "<<d_reference[__a][__b][__c]<<std::endl;std::cerr<<"d Rel Difference=  "<< std::abs((d[__a][__b][__c] - d_reference[__a][__b][__c]) / (d_reference[__a][__b][__c])) << std::endl;std::cerr<<"d Abs Difference=  "<< std::abs(d[__a][__b][__c] - d_reference[__a][__b][__c]) << std::endl;return 1;}
for(int __a=0;__a<300;__a++) for(int __b=0;__b<16;__b++) if(std::abs((system_matrix[__a][__b] - system_matrix_reference[__a][__b]) / (system_matrix_reference[__a][__b])) > 1 ){std::cerr<<"Mismatch detected in SCALAR implementation"<<std::endl;std::cerr<<"Variable system_matrix:"<<std::endl;std::cerr<<"seed="<<seed<<", __a="<<__a<<", __b="<<__b<<std::endl;std::cerr<<"system_matrix SCALAR=  "<<system_matrix[__a][__b]<<std::endl;std::cerr<<"system_matrix Reference=  "<<system_matrix_reference[__a][__b]<<std::endl;std::cerr<<"system_matrix Rel Difference=  "<< std::abs((system_matrix[__a][__b] - system_matrix_reference[__a][__b]) / (system_matrix_reference[__a][__b])) << std::endl;std::cerr<<"system_matrix Abs Difference=  "<< std::abs(system_matrix[__a][__b] - system_matrix_reference[__a][__b]) << std::endl;return 1;}

}

//=======================================================
//
//               COMPUTE SSE RESULTS
//
//=======================================================

#ifdef ENABLE_SSE_INSTRUCTION_SET
{
typedef T (&refArray1)[3][8][16];typedef T (&refArray2)[16];typedef T (&refArray3)[16];typedef T (&refArray4)[16];typedef T (&refArray5)[16];typedef T (&refArray6)[16];typedef T (&refArray7)[16];typedef T (&refArray8)[16];typedef T (&refArray9)[16];typedef T (&refArray10)[9][16];typedef T (&refArray11)[9][16];typedef T (&refArray12)[3][16];typedef T (&refArray13)[3][16];typedef T (&refArray14)[12][16];typedef T (&refArray15)[3][8][16];typedef T (&refArray16)[300][16];for(int __a=0;__a<9;__a++) for(int __b=0;__b<16;__b++) U[__a][__b]=U_original[__a][__b];for(int __a=0;__a<9;__a++) for(int __b=0;__b<16;__b++) V[__a][__b]=V_original[__a][__b];for(int __a=0;__a<3;__a++) for(int __b=0;__b<16;__b++) Sigma[__a][__b]=Sigma_original[__a][__b];for(int __a=0;__a<3;__a++) for(int __b=0;__b<16;__b++) Q_hat[__a][__b]=Q_hat_original[__a][__b];for(int __a=0;__a<12;__a++) for(int __b=0;__b<16;__b++) dPdF[__a][__b]=dPdF_original[__a][__b];for(int __a=0;__a<3;__a++) for(int __b=0;__b<8;__b++) for(int __c=0;__c<16;__c++) d[__a][__b][__c]=d_original[__a][__b][__c];for(int __a=0;__a<300;__a++) for(int __b=0;__b<16;__b++) system_matrix[__a][__b]=system_matrix_original[__a][__b];for(int i=0; i<16; i+=4) {refArray1 uk=reinterpret_cast<refArray1>(u[0][0][i]);refArray2 pk=reinterpret_cast<refArray2>(p[i]);refArray3 muk=reinterpret_cast<refArray3>(mu[i]);refArray4 mu_stabk=reinterpret_cast<refArray4>(mu_stab[i]);refArray5 kappak=reinterpret_cast<refArray5>(kappa[i]);refArray6 alphak=reinterpret_cast<refArray6>(alpha[i]);refArray7 cutoffk=reinterpret_cast<refArray7>(cutoff[i]);refArray8 one_over_hk=reinterpret_cast<refArray8>(one_over_h[i]);refArray9 cell_volumek=reinterpret_cast<refArray9>(cell_volume[i]);refArray10 Uk=reinterpret_cast<refArray10>(U[0][i]);refArray11 Vk=reinterpret_cast<refArray11>(V[0][i]);refArray12 Sigmak=reinterpret_cast<refArray12>(Sigma[0][i]);refArray13 Q_hatk=reinterpret_cast<refArray13>(Q_hat[0][i]);refArray14 dPdFk=reinterpret_cast<refArray14>(dPdF[0][i]);refArray15 dk=reinterpret_cast<refArray15>(d[0][0][i]);refArray16 system_matrixk=reinterpret_cast<refArray16>(system_matrix[0][i]);Update_Position_Based_State<NEOHOOKEAN_TAG,__m128,float[16],int[16]>::Run(uk,pk,muk,mu_stabk,kappak,alphak,cutoffk,one_over_hk,cell_volumek,Uk,Vk,Sigmak,Q_hatk,dPdFk,dk,system_matrixk);}
for(int __a=0;__a<9;__a++) for(int __b=0;__b<16;__b++) if(std::abs((U[__a][__b] - U_reference[__a][__b]) / (U_reference[__a][__b])) > 1 ){std::cerr<<"Mismatch detected in SSE implementation"<<std::endl;std::cerr<<"Variable U:"<<std::endl;std::cerr<<"seed="<<seed<<", __a="<<__a<<", __b="<<__b<<std::endl;std::cerr<<"U SSE=  "<<U[__a][__b]<<std::endl;std::cerr<<"U Reference=  "<<U_reference[__a][__b]<<std::endl;std::cerr<<"U Rel Difference=  "<< std::abs((U[__a][__b] - U_reference[__a][__b]) / (U_reference[__a][__b])) << std::endl;std::cerr<<"U Abs Difference=  "<< std::abs(U[__a][__b] - U_reference[__a][__b]) << std::endl;return 1;}
for(int __a=0;__a<9;__a++) for(int __b=0;__b<16;__b++) if(std::abs((V[__a][__b] - V_reference[__a][__b]) / (V_reference[__a][__b])) > 1 ){std::cerr<<"Mismatch detected in SSE implementation"<<std::endl;std::cerr<<"Variable V:"<<std::endl;std::cerr<<"seed="<<seed<<", __a="<<__a<<", __b="<<__b<<std::endl;std::cerr<<"V SSE=  "<<V[__a][__b]<<std::endl;std::cerr<<"V Reference=  "<<V_reference[__a][__b]<<std::endl;std::cerr<<"V Rel Difference=  "<< std::abs((V[__a][__b] - V_reference[__a][__b]) / (V_reference[__a][__b])) << std::endl;std::cerr<<"V Abs Difference=  "<< std::abs(V[__a][__b] - V_reference[__a][__b]) << std::endl;return 1;}
for(int __a=0;__a<3;__a++) for(int __b=0;__b<16;__b++) if(std::abs((Sigma[__a][__b] - Sigma_reference[__a][__b]) / (Sigma_reference[__a][__b])) > 1 ){std::cerr<<"Mismatch detected in SSE implementation"<<std::endl;std::cerr<<"Variable Sigma:"<<std::endl;std::cerr<<"seed="<<seed<<", __a="<<__a<<", __b="<<__b<<std::endl;std::cerr<<"Sigma SSE=  "<<Sigma[__a][__b]<<std::endl;std::cerr<<"Sigma Reference=  "<<Sigma_reference[__a][__b]<<std::endl;std::cerr<<"Sigma Rel Difference=  "<< std::abs((Sigma[__a][__b] - Sigma_reference[__a][__b]) / (Sigma_reference[__a][__b])) << std::endl;std::cerr<<"Sigma Abs Difference=  "<< std::abs(Sigma[__a][__b] - Sigma_reference[__a][__b]) << std::endl;return 1;}
for(int __a=0;__a<3;__a++) for(int __b=0;__b<16;__b++) if(std::abs((Q_hat[__a][__b] - Q_hat_reference[__a][__b]) / (Q_hat_reference[__a][__b])) > 1 ){std::cerr<<"Mismatch detected in SSE implementation"<<std::endl;std::cerr<<"Variable Q_hat:"<<std::endl;std::cerr<<"seed="<<seed<<", __a="<<__a<<", __b="<<__b<<std::endl;std::cerr<<"Q_hat SSE=  "<<Q_hat[__a][__b]<<std::endl;std::cerr<<"Q_hat Reference=  "<<Q_hat_reference[__a][__b]<<std::endl;std::cerr<<"Q_hat Rel Difference=  "<< std::abs((Q_hat[__a][__b] - Q_hat_reference[__a][__b]) / (Q_hat_reference[__a][__b])) << std::endl;std::cerr<<"Q_hat Abs Difference=  "<< std::abs(Q_hat[__a][__b] - Q_hat_reference[__a][__b]) << std::endl;return 1;}
for(int __a=0;__a<12;__a++) for(int __b=0;__b<16;__b++) if(std::abs((dPdF[__a][__b] - dPdF_reference[__a][__b]) / (dPdF_reference[__a][__b])) > 1 ){std::cerr<<"Mismatch detected in SSE implementation"<<std::endl;std::cerr<<"Variable dPdF:"<<std::endl;std::cerr<<"seed="<<seed<<", __a="<<__a<<", __b="<<__b<<std::endl;std::cerr<<"dPdF SSE=  "<<dPdF[__a][__b]<<std::endl;std::cerr<<"dPdF Reference=  "<<dPdF_reference[__a][__b]<<std::endl;std::cerr<<"dPdF Rel Difference=  "<< std::abs((dPdF[__a][__b] - dPdF_reference[__a][__b]) / (dPdF_reference[__a][__b])) << std::endl;std::cerr<<"dPdF Abs Difference=  "<< std::abs(dPdF[__a][__b] - dPdF_reference[__a][__b]) << std::endl;return 1;}
for(int __a=0;__a<3;__a++) for(int __b=0;__b<8;__b++) for(int __c=0;__c<16;__c++) if(std::abs((d[__a][__b][__c] - d_reference[__a][__b][__c]) / (d_reference[__a][__b][__c])) > 1 ){std::cerr<<"Mismatch detected in SSE implementation"<<std::endl;std::cerr<<"Variable d:"<<std::endl;std::cerr<<"seed="<<seed<<", __a="<<__a<<", __b="<<__b<<", __c="<<__c<<std::endl;std::cerr<<"d SSE=  "<<d[__a][__b][__c]<<std::endl;std::cerr<<"d Reference=  "<<d_reference[__a][__b][__c]<<std::endl;std::cerr<<"d Rel Difference=  "<< std::abs((d[__a][__b][__c] - d_reference[__a][__b][__c]) / (d_reference[__a][__b][__c])) << std::endl;std::cerr<<"d Abs Difference=  "<< std::abs(d[__a][__b][__c] - d_reference[__a][__b][__c]) << std::endl;return 1;}
for(int __a=0;__a<300;__a++) for(int __b=0;__b<16;__b++) if(std::abs((system_matrix[__a][__b] - system_matrix_reference[__a][__b]) / (system_matrix_reference[__a][__b])) > 1 ){std::cerr<<"Mismatch detected in SSE implementation"<<std::endl;std::cerr<<"Variable system_matrix:"<<std::endl;std::cerr<<"seed="<<seed<<", __a="<<__a<<", __b="<<__b<<std::endl;std::cerr<<"system_matrix SSE=  "<<system_matrix[__a][__b]<<std::endl;std::cerr<<"system_matrix Reference=  "<<system_matrix_reference[__a][__b]<<std::endl;std::cerr<<"system_matrix Rel Difference=  "<< std::abs((system_matrix[__a][__b] - system_matrix_reference[__a][__b]) / (system_matrix_reference[__a][__b])) << std::endl;std::cerr<<"system_matrix Abs Difference=  "<< std::abs(system_matrix[__a][__b] - system_matrix_reference[__a][__b]) << std::endl;return 1;}

}
#endif

//=======================================================
//
//               COMPUTE AVX RESULTS
//
//=======================================================

#ifdef ENABLE_AVX_INSTRUCTION_SET
{
typedef T (&refArray1)[3][8][16];typedef T (&refArray2)[16];typedef T (&refArray3)[16];typedef T (&refArray4)[16];typedef T (&refArray5)[16];typedef T (&refArray6)[16];typedef T (&refArray7)[16];typedef T (&refArray8)[16];typedef T (&refArray9)[16];typedef T (&refArray10)[9][16];typedef T (&refArray11)[9][16];typedef T (&refArray12)[3][16];typedef T (&refArray13)[3][16];typedef T (&refArray14)[12][16];typedef T (&refArray15)[3][8][16];typedef T (&refArray16)[300][16];for(int __a=0;__a<9;__a++) for(int __b=0;__b<16;__b++) U[__a][__b]=U_original[__a][__b];for(int __a=0;__a<9;__a++) for(int __b=0;__b<16;__b++) V[__a][__b]=V_original[__a][__b];for(int __a=0;__a<3;__a++) for(int __b=0;__b<16;__b++) Sigma[__a][__b]=Sigma_original[__a][__b];for(int __a=0;__a<3;__a++) for(int __b=0;__b<16;__b++) Q_hat[__a][__b]=Q_hat_original[__a][__b];for(int __a=0;__a<12;__a++) for(int __b=0;__b<16;__b++) dPdF[__a][__b]=dPdF_original[__a][__b];for(int __a=0;__a<3;__a++) for(int __b=0;__b<8;__b++) for(int __c=0;__c<16;__c++) d[__a][__b][__c]=d_original[__a][__b][__c];for(int __a=0;__a<300;__a++) for(int __b=0;__b<16;__b++) system_matrix[__a][__b]=system_matrix_original[__a][__b];for(int i=0; i<16; i+=8) {refArray1 uk=reinterpret_cast<refArray1>(u[0][0][i]);refArray2 pk=reinterpret_cast<refArray2>(p[i]);refArray3 muk=reinterpret_cast<refArray3>(mu[i]);refArray4 mu_stabk=reinterpret_cast<refArray4>(mu_stab[i]);refArray5 kappak=reinterpret_cast<refArray5>(kappa[i]);refArray6 alphak=reinterpret_cast<refArray6>(alpha[i]);refArray7 cutoffk=reinterpret_cast<refArray7>(cutoff[i]);refArray8 one_over_hk=reinterpret_cast<refArray8>(one_over_h[i]);refArray9 cell_volumek=reinterpret_cast<refArray9>(cell_volume[i]);refArray10 Uk=reinterpret_cast<refArray10>(U[0][i]);refArray11 Vk=reinterpret_cast<refArray11>(V[0][i]);refArray12 Sigmak=reinterpret_cast<refArray12>(Sigma[0][i]);refArray13 Q_hatk=reinterpret_cast<refArray13>(Q_hat[0][i]);refArray14 dPdFk=reinterpret_cast<refArray14>(dPdF[0][i]);refArray15 dk=reinterpret_cast<refArray15>(d[0][0][i]);refArray16 system_matrixk=reinterpret_cast<refArray16>(system_matrix[0][i]);Update_Position_Based_State<NEOHOOKEAN_TAG,__m256,float[16],int[16]>::Run(uk,pk,muk,mu_stabk,kappak,alphak,cutoffk,one_over_hk,cell_volumek,Uk,Vk,Sigmak,Q_hatk,dPdFk,dk,system_matrixk);}
for(int __a=0;__a<9;__a++) for(int __b=0;__b<16;__b++) if(std::abs((U[__a][__b] - U_reference[__a][__b]) / (U_reference[__a][__b])) > 1 ){std::cerr<<"Mismatch detected in AVX implementation"<<std::endl;std::cerr<<"Variable U:"<<std::endl;std::cerr<<"seed="<<seed<<", __a="<<__a<<", __b="<<__b<<std::endl;std::cerr<<"U AVX=  "<<U[__a][__b]<<std::endl;std::cerr<<"U Reference=  "<<U_reference[__a][__b]<<std::endl;std::cerr<<"U Rel Difference=  "<< std::abs((U[__a][__b] - U_reference[__a][__b]) / (U_reference[__a][__b])) << std::endl;std::cerr<<"U Abs Difference=  "<< std::abs(U[__a][__b] - U_reference[__a][__b]) << std::endl;return 1;}
for(int __a=0;__a<9;__a++) for(int __b=0;__b<16;__b++) if(std::abs((V[__a][__b] - V_reference[__a][__b]) / (V_reference[__a][__b])) > 1 ){std::cerr<<"Mismatch detected in AVX implementation"<<std::endl;std::cerr<<"Variable V:"<<std::endl;std::cerr<<"seed="<<seed<<", __a="<<__a<<", __b="<<__b<<std::endl;std::cerr<<"V AVX=  "<<V[__a][__b]<<std::endl;std::cerr<<"V Reference=  "<<V_reference[__a][__b]<<std::endl;std::cerr<<"V Rel Difference=  "<< std::abs((V[__a][__b] - V_reference[__a][__b]) / (V_reference[__a][__b])) << std::endl;std::cerr<<"V Abs Difference=  "<< std::abs(V[__a][__b] - V_reference[__a][__b]) << std::endl;return 1;}
for(int __a=0;__a<3;__a++) for(int __b=0;__b<16;__b++) if(std::abs((Sigma[__a][__b] - Sigma_reference[__a][__b]) / (Sigma_reference[__a][__b])) > 1 ){std::cerr<<"Mismatch detected in AVX implementation"<<std::endl;std::cerr<<"Variable Sigma:"<<std::endl;std::cerr<<"seed="<<seed<<", __a="<<__a<<", __b="<<__b<<std::endl;std::cerr<<"Sigma AVX=  "<<Sigma[__a][__b]<<std::endl;std::cerr<<"Sigma Reference=  "<<Sigma_reference[__a][__b]<<std::endl;std::cerr<<"Sigma Rel Difference=  "<< std::abs((Sigma[__a][__b] - Sigma_reference[__a][__b]) / (Sigma_reference[__a][__b])) << std::endl;std::cerr<<"Sigma Abs Difference=  "<< std::abs(Sigma[__a][__b] - Sigma_reference[__a][__b]) << std::endl;return 1;}
for(int __a=0;__a<3;__a++) for(int __b=0;__b<16;__b++) if(std::abs((Q_hat[__a][__b] - Q_hat_reference[__a][__b]) / (Q_hat_reference[__a][__b])) > 1 ){std::cerr<<"Mismatch detected in AVX implementation"<<std::endl;std::cerr<<"Variable Q_hat:"<<std::endl;std::cerr<<"seed="<<seed<<", __a="<<__a<<", __b="<<__b<<std::endl;std::cerr<<"Q_hat AVX=  "<<Q_hat[__a][__b]<<std::endl;std::cerr<<"Q_hat Reference=  "<<Q_hat_reference[__a][__b]<<std::endl;std::cerr<<"Q_hat Rel Difference=  "<< std::abs((Q_hat[__a][__b] - Q_hat_reference[__a][__b]) / (Q_hat_reference[__a][__b])) << std::endl;std::cerr<<"Q_hat Abs Difference=  "<< std::abs(Q_hat[__a][__b] - Q_hat_reference[__a][__b]) << std::endl;return 1;}
for(int __a=0;__a<12;__a++) for(int __b=0;__b<16;__b++) if(std::abs((dPdF[__a][__b] - dPdF_reference[__a][__b]) / (dPdF_reference[__a][__b])) > 1 ){std::cerr<<"Mismatch detected in AVX implementation"<<std::endl;std::cerr<<"Variable dPdF:"<<std::endl;std::cerr<<"seed="<<seed<<", __a="<<__a<<", __b="<<__b<<std::endl;std::cerr<<"dPdF AVX=  "<<dPdF[__a][__b]<<std::endl;std::cerr<<"dPdF Reference=  "<<dPdF_reference[__a][__b]<<std::endl;std::cerr<<"dPdF Rel Difference=  "<< std::abs((dPdF[__a][__b] - dPdF_reference[__a][__b]) / (dPdF_reference[__a][__b])) << std::endl;std::cerr<<"dPdF Abs Difference=  "<< std::abs(dPdF[__a][__b] - dPdF_reference[__a][__b]) << std::endl;return 1;}
for(int __a=0;__a<3;__a++) for(int __b=0;__b<8;__b++) for(int __c=0;__c<16;__c++) if(std::abs((d[__a][__b][__c] - d_reference[__a][__b][__c]) / (d_reference[__a][__b][__c])) > 1 ){std::cerr<<"Mismatch detected in AVX implementation"<<std::endl;std::cerr<<"Variable d:"<<std::endl;std::cerr<<"seed="<<seed<<", __a="<<__a<<", __b="<<__b<<", __c="<<__c<<std::endl;std::cerr<<"d AVX=  "<<d[__a][__b][__c]<<std::endl;std::cerr<<"d Reference=  "<<d_reference[__a][__b][__c]<<std::endl;std::cerr<<"d Rel Difference=  "<< std::abs((d[__a][__b][__c] - d_reference[__a][__b][__c]) / (d_reference[__a][__b][__c])) << std::endl;std::cerr<<"d Abs Difference=  "<< std::abs(d[__a][__b][__c] - d_reference[__a][__b][__c]) << std::endl;return 1;}
for(int __a=0;__a<300;__a++) for(int __b=0;__b<16;__b++) if(std::abs((system_matrix[__a][__b] - system_matrix_reference[__a][__b]) / (system_matrix_reference[__a][__b])) > 1 ){std::cerr<<"Mismatch detected in AVX implementation"<<std::endl;std::cerr<<"Variable system_matrix:"<<std::endl;std::cerr<<"seed="<<seed<<", __a="<<__a<<", __b="<<__b<<std::endl;std::cerr<<"system_matrix AVX=  "<<system_matrix[__a][__b]<<std::endl;std::cerr<<"system_matrix Reference=  "<<system_matrix_reference[__a][__b]<<std::endl;std::cerr<<"system_matrix Rel Difference=  "<< std::abs((system_matrix[__a][__b] - system_matrix_reference[__a][__b]) / (system_matrix_reference[__a][__b])) << std::endl;std::cerr<<"system_matrix Abs Difference=  "<< std::abs(system_matrix[__a][__b] - system_matrix_reference[__a][__b]) << std::endl;return 1;}

}
#endif

//=======================================================
//
//               COMPUTE NEON RESULTS
//
//=======================================================

#ifdef ENABLE_NEON_INSTRUCTION_SET
{
typedef T (&refArray1)[3][8][16];typedef T (&refArray2)[16];typedef T (&refArray3)[16];typedef T (&refArray4)[16];typedef T (&refArray5)[16];typedef T (&refArray6)[16];typedef T (&refArray7)[16];typedef T (&refArray8)[16];typedef T (&refArray9)[16];typedef T (&refArray10)[9][16];typedef T (&refArray11)[9][16];typedef T (&refArray12)[3][16];typedef T (&refArray13)[3][16];typedef T (&refArray14)[12][16];typedef T (&refArray15)[3][8][16];typedef T (&refArray16)[300][16];for(int __a=0;__a<9;__a++) for(int __b=0;__b<16;__b++) U[__a][__b]=U_original[__a][__b];for(int __a=0;__a<9;__a++) for(int __b=0;__b<16;__b++) V[__a][__b]=V_original[__a][__b];for(int __a=0;__a<3;__a++) for(int __b=0;__b<16;__b++) Sigma[__a][__b]=Sigma_original[__a][__b];for(int __a=0;__a<3;__a++) for(int __b=0;__b<16;__b++) Q_hat[__a][__b]=Q_hat_original[__a][__b];for(int __a=0;__a<12;__a++) for(int __b=0;__b<16;__b++) dPdF[__a][__b]=dPdF_original[__a][__b];for(int __a=0;__a<3;__a++) for(int __b=0;__b<8;__b++) for(int __c=0;__c<16;__c++) d[__a][__b][__c]=d_original[__a][__b][__c];for(int __a=0;__a<300;__a++) for(int __b=0;__b<16;__b++) system_matrix[__a][__b]=system_matrix_original[__a][__b];for(int i=0; i<16; i+=4) {refArray1 uk=reinterpret_cast<refArray1>(u[0][0][i]);refArray2 pk=reinterpret_cast<refArray2>(p[i]);refArray3 muk=reinterpret_cast<refArray3>(mu[i]);refArray4 mu_stabk=reinterpret_cast<refArray4>(mu_stab[i]);refArray5 kappak=reinterpret_cast<refArray5>(kappa[i]);refArray6 alphak=reinterpret_cast<refArray6>(alpha[i]);refArray7 cutoffk=reinterpret_cast<refArray7>(cutoff[i]);refArray8 one_over_hk=reinterpret_cast<refArray8>(one_over_h[i]);refArray9 cell_volumek=reinterpret_cast<refArray9>(cell_volume[i]);refArray10 Uk=reinterpret_cast<refArray10>(U[0][i]);refArray11 Vk=reinterpret_cast<refArray11>(V[0][i]);refArray12 Sigmak=reinterpret_cast<refArray12>(Sigma[0][i]);refArray13 Q_hatk=reinterpret_cast<refArray13>(Q_hat[0][i]);refArray14 dPdFk=reinterpret_cast<refArray14>(dPdF[0][i]);refArray15 dk=reinterpret_cast<refArray15>(d[0][0][i]);refArray16 system_matrixk=reinterpret_cast<refArray16>(system_matrix[0][i]);Update_Position_Based_State<NEOHOOKEAN_TAG,float32x4_t,float[16],int[16]>::Run(uk,pk,muk,mu_stabk,kappak,alphak,cutoffk,one_over_hk,cell_volumek,Uk,Vk,Sigmak,Q_hatk,dPdFk,dk,system_matrixk);}
for(int __a=0;__a<9;__a++) for(int __b=0;__b<16;__b++) if(std::abs((U[__a][__b] - U_reference[__a][__b]) / (U_reference[__a][__b])) > 1 ){std::cerr<<"Mismatch detected in NEON implementation"<<std::endl;std::cerr<<"Variable U:"<<std::endl;std::cerr<<"seed="<<seed<<", __a="<<__a<<", __b="<<__b<<std::endl;std::cerr<<"U NEON=  "<<U[__a][__b]<<std::endl;std::cerr<<"U Reference=  "<<U_reference[__a][__b]<<std::endl;std::cerr<<"U Rel Difference=  "<< std::abs((U[__a][__b] - U_reference[__a][__b]) / (U_reference[__a][__b])) << std::endl;std::cerr<<"U Abs Difference=  "<< std::abs(U[__a][__b] - U_reference[__a][__b]) << std::endl;return 1;}
for(int __a=0;__a<9;__a++) for(int __b=0;__b<16;__b++) if(std::abs((V[__a][__b] - V_reference[__a][__b]) / (V_reference[__a][__b])) > 1 ){std::cerr<<"Mismatch detected in NEON implementation"<<std::endl;std::cerr<<"Variable V:"<<std::endl;std::cerr<<"seed="<<seed<<", __a="<<__a<<", __b="<<__b<<std::endl;std::cerr<<"V NEON=  "<<V[__a][__b]<<std::endl;std::cerr<<"V Reference=  "<<V_reference[__a][__b]<<std::endl;std::cerr<<"V Rel Difference=  "<< std::abs((V[__a][__b] - V_reference[__a][__b]) / (V_reference[__a][__b])) << std::endl;std::cerr<<"V Abs Difference=  "<< std::abs(V[__a][__b] - V_reference[__a][__b]) << std::endl;return 1;}
for(int __a=0;__a<3;__a++) for(int __b=0;__b<16;__b++) if(std::abs((Sigma[__a][__b] - Sigma_reference[__a][__b]) / (Sigma_reference[__a][__b])) > 1 ){std::cerr<<"Mismatch detected in NEON implementation"<<std::endl;std::cerr<<"Variable Sigma:"<<std::endl;std::cerr<<"seed="<<seed<<", __a="<<__a<<", __b="<<__b<<std::endl;std::cerr<<"Sigma NEON=  "<<Sigma[__a][__b]<<std::endl;std::cerr<<"Sigma Reference=  "<<Sigma_reference[__a][__b]<<std::endl;std::cerr<<"Sigma Rel Difference=  "<< std::abs((Sigma[__a][__b] - Sigma_reference[__a][__b]) / (Sigma_reference[__a][__b])) << std::endl;std::cerr<<"Sigma Abs Difference=  "<< std::abs(Sigma[__a][__b] - Sigma_reference[__a][__b]) << std::endl;return 1;}
for(int __a=0;__a<3;__a++) for(int __b=0;__b<16;__b++) if(std::abs((Q_hat[__a][__b] - Q_hat_reference[__a][__b]) / (Q_hat_reference[__a][__b])) > 1 ){std::cerr<<"Mismatch detected in NEON implementation"<<std::endl;std::cerr<<"Variable Q_hat:"<<std::endl;std::cerr<<"seed="<<seed<<", __a="<<__a<<", __b="<<__b<<std::endl;std::cerr<<"Q_hat NEON=  "<<Q_hat[__a][__b]<<std::endl;std::cerr<<"Q_hat Reference=  "<<Q_hat_reference[__a][__b]<<std::endl;std::cerr<<"Q_hat Rel Difference=  "<< std::abs((Q_hat[__a][__b] - Q_hat_reference[__a][__b]) / (Q_hat_reference[__a][__b])) << std::endl;std::cerr<<"Q_hat Abs Difference=  "<< std::abs(Q_hat[__a][__b] - Q_hat_reference[__a][__b]) << std::endl;return 1;}
for(int __a=0;__a<12;__a++) for(int __b=0;__b<16;__b++) if(std::abs((dPdF[__a][__b] - dPdF_reference[__a][__b]) / (dPdF_reference[__a][__b])) > 1 ){std::cerr<<"Mismatch detected in NEON implementation"<<std::endl;std::cerr<<"Variable dPdF:"<<std::endl;std::cerr<<"seed="<<seed<<", __a="<<__a<<", __b="<<__b<<std::endl;std::cerr<<"dPdF NEON=  "<<dPdF[__a][__b]<<std::endl;std::cerr<<"dPdF Reference=  "<<dPdF_reference[__a][__b]<<std::endl;std::cerr<<"dPdF Rel Difference=  "<< std::abs((dPdF[__a][__b] - dPdF_reference[__a][__b]) / (dPdF_reference[__a][__b])) << std::endl;std::cerr<<"dPdF Abs Difference=  "<< std::abs(dPdF[__a][__b] - dPdF_reference[__a][__b]) << std::endl;return 1;}
for(int __a=0;__a<3;__a++) for(int __b=0;__b<8;__b++) for(int __c=0;__c<16;__c++) if(std::abs((d[__a][__b][__c] - d_reference[__a][__b][__c]) / (d_reference[__a][__b][__c])) > 1 ){std::cerr<<"Mismatch detected in NEON implementation"<<std::endl;std::cerr<<"Variable d:"<<std::endl;std::cerr<<"seed="<<seed<<", __a="<<__a<<", __b="<<__b<<", __c="<<__c<<std::endl;std::cerr<<"d NEON=  "<<d[__a][__b][__c]<<std::endl;std::cerr<<"d Reference=  "<<d_reference[__a][__b][__c]<<std::endl;std::cerr<<"d Rel Difference=  "<< std::abs((d[__a][__b][__c] - d_reference[__a][__b][__c]) / (d_reference[__a][__b][__c])) << std::endl;std::cerr<<"d Abs Difference=  "<< std::abs(d[__a][__b][__c] - d_reference[__a][__b][__c]) << std::endl;return 1;}
for(int __a=0;__a<300;__a++) for(int __b=0;__b<16;__b++) if(std::abs((system_matrix[__a][__b] - system_matrix_reference[__a][__b]) / (system_matrix_reference[__a][__b])) > 1 ){std::cerr<<"Mismatch detected in NEON implementation"<<std::endl;std::cerr<<"Variable system_matrix:"<<std::endl;std::cerr<<"seed="<<seed<<", __a="<<__a<<", __b="<<__b<<std::endl;std::cerr<<"system_matrix NEON=  "<<system_matrix[__a][__b]<<std::endl;std::cerr<<"system_matrix Reference=  "<<system_matrix_reference[__a][__b]<<std::endl;std::cerr<<"system_matrix Rel Difference=  "<< std::abs((system_matrix[__a][__b] - system_matrix_reference[__a][__b]) / (system_matrix_reference[__a][__b])) << std::endl;std::cerr<<"system_matrix Abs Difference=  "<< std::abs(system_matrix[__a][__b] - system_matrix_reference[__a][__b]) << std::endl;return 1;}

}
#endif

//=======================================================
//
//               COMPUTE MIC RESULTS
//
//=======================================================

#ifdef ENABLE_MIC_INSTRUCTION_SET
{
typedef T (&refArray1)[3][8][16];typedef T (&refArray2)[16];typedef T (&refArray3)[16];typedef T (&refArray4)[16];typedef T (&refArray5)[16];typedef T (&refArray6)[16];typedef T (&refArray7)[16];typedef T (&refArray8)[16];typedef T (&refArray9)[16];typedef T (&refArray10)[9][16];typedef T (&refArray11)[9][16];typedef T (&refArray12)[3][16];typedef T (&refArray13)[3][16];typedef T (&refArray14)[12][16];typedef T (&refArray15)[3][8][16];typedef T (&refArray16)[300][16];for(int __a=0;__a<9;__a++) for(int __b=0;__b<16;__b++) U[__a][__b]=U_original[__a][__b];for(int __a=0;__a<9;__a++) for(int __b=0;__b<16;__b++) V[__a][__b]=V_original[__a][__b];for(int __a=0;__a<3;__a++) for(int __b=0;__b<16;__b++) Sigma[__a][__b]=Sigma_original[__a][__b];for(int __a=0;__a<3;__a++) for(int __b=0;__b<16;__b++) Q_hat[__a][__b]=Q_hat_original[__a][__b];for(int __a=0;__a<12;__a++) for(int __b=0;__b<16;__b++) dPdF[__a][__b]=dPdF_original[__a][__b];for(int __a=0;__a<3;__a++) for(int __b=0;__b<8;__b++) for(int __c=0;__c<16;__c++) d[__a][__b][__c]=d_original[__a][__b][__c];for(int __a=0;__a<300;__a++) for(int __b=0;__b<16;__b++) system_matrix[__a][__b]=system_matrix_original[__a][__b];for(int i=0; i<16; i+=16) {refArray1 uk=reinterpret_cast<refArray1>(u[0][0][i]);refArray2 pk=reinterpret_cast<refArray2>(p[i]);refArray3 muk=reinterpret_cast<refArray3>(mu[i]);refArray4 mu_stabk=reinterpret_cast<refArray4>(mu_stab[i]);refArray5 kappak=reinterpret_cast<refArray5>(kappa[i]);refArray6 alphak=reinterpret_cast<refArray6>(alpha[i]);refArray7 cutoffk=reinterpret_cast<refArray7>(cutoff[i]);refArray8 one_over_hk=reinterpret_cast<refArray8>(one_over_h[i]);refArray9 cell_volumek=reinterpret_cast<refArray9>(cell_volume[i]);refArray10 Uk=reinterpret_cast<refArray10>(U[0][i]);refArray11 Vk=reinterpret_cast<refArray11>(V[0][i]);refArray12 Sigmak=reinterpret_cast<refArray12>(Sigma[0][i]);refArray13 Q_hatk=reinterpret_cast<refArray13>(Q_hat[0][i]);refArray14 dPdFk=reinterpret_cast<refArray14>(dPdF[0][i]);refArray15 dk=reinterpret_cast<refArray15>(d[0][0][i]);refArray16 system_matrixk=reinterpret_cast<refArray16>(system_matrix[0][i]);Update_Position_Based_State<NEOHOOKEAN_TAG,__m512,float[16],int[16]>::Run(uk,pk,muk,mu_stabk,kappak,alphak,cutoffk,one_over_hk,cell_volumek,Uk,Vk,Sigmak,Q_hatk,dPdFk,dk,system_matrixk);}
for(int __a=0;__a<9;__a++) for(int __b=0;__b<16;__b++) if(std::abs((U[__a][__b] - U_reference[__a][__b]) / (U_reference[__a][__b])) > 1 ){std::cerr<<"Mismatch detected in MIC implementation"<<std::endl;std::cerr<<"Variable U:"<<std::endl;std::cerr<<"seed="<<seed<<", __a="<<__a<<", __b="<<__b<<std::endl;std::cerr<<"U MIC=  "<<U[__a][__b]<<std::endl;std::cerr<<"U Reference=  "<<U_reference[__a][__b]<<std::endl;std::cerr<<"U Rel Difference=  "<< std::abs((U[__a][__b] - U_reference[__a][__b]) / (U_reference[__a][__b])) << std::endl;std::cerr<<"U Abs Difference=  "<< std::abs(U[__a][__b] - U_reference[__a][__b]) << std::endl;return 1;}
for(int __a=0;__a<9;__a++) for(int __b=0;__b<16;__b++) if(std::abs((V[__a][__b] - V_reference[__a][__b]) / (V_reference[__a][__b])) > 1 ){std::cerr<<"Mismatch detected in MIC implementation"<<std::endl;std::cerr<<"Variable V:"<<std::endl;std::cerr<<"seed="<<seed<<", __a="<<__a<<", __b="<<__b<<std::endl;std::cerr<<"V MIC=  "<<V[__a][__b]<<std::endl;std::cerr<<"V Reference=  "<<V_reference[__a][__b]<<std::endl;std::cerr<<"V Rel Difference=  "<< std::abs((V[__a][__b] - V_reference[__a][__b]) / (V_reference[__a][__b])) << std::endl;std::cerr<<"V Abs Difference=  "<< std::abs(V[__a][__b] - V_reference[__a][__b]) << std::endl;return 1;}
for(int __a=0;__a<3;__a++) for(int __b=0;__b<16;__b++) if(std::abs((Sigma[__a][__b] - Sigma_reference[__a][__b]) / (Sigma_reference[__a][__b])) > 1 ){std::cerr<<"Mismatch detected in MIC implementation"<<std::endl;std::cerr<<"Variable Sigma:"<<std::endl;std::cerr<<"seed="<<seed<<", __a="<<__a<<", __b="<<__b<<std::endl;std::cerr<<"Sigma MIC=  "<<Sigma[__a][__b]<<std::endl;std::cerr<<"Sigma Reference=  "<<Sigma_reference[__a][__b]<<std::endl;std::cerr<<"Sigma Rel Difference=  "<< std::abs((Sigma[__a][__b] - Sigma_reference[__a][__b]) / (Sigma_reference[__a][__b])) << std::endl;std::cerr<<"Sigma Abs Difference=  "<< std::abs(Sigma[__a][__b] - Sigma_reference[__a][__b]) << std::endl;return 1;}
for(int __a=0;__a<3;__a++) for(int __b=0;__b<16;__b++) if(std::abs((Q_hat[__a][__b] - Q_hat_reference[__a][__b]) / (Q_hat_reference[__a][__b])) > 1 ){std::cerr<<"Mismatch detected in MIC implementation"<<std::endl;std::cerr<<"Variable Q_hat:"<<std::endl;std::cerr<<"seed="<<seed<<", __a="<<__a<<", __b="<<__b<<std::endl;std::cerr<<"Q_hat MIC=  "<<Q_hat[__a][__b]<<std::endl;std::cerr<<"Q_hat Reference=  "<<Q_hat_reference[__a][__b]<<std::endl;std::cerr<<"Q_hat Rel Difference=  "<< std::abs((Q_hat[__a][__b] - Q_hat_reference[__a][__b]) / (Q_hat_reference[__a][__b])) << std::endl;std::cerr<<"Q_hat Abs Difference=  "<< std::abs(Q_hat[__a][__b] - Q_hat_reference[__a][__b]) << std::endl;return 1;}
for(int __a=0;__a<12;__a++) for(int __b=0;__b<16;__b++) if(std::abs((dPdF[__a][__b] - dPdF_reference[__a][__b]) / (dPdF_reference[__a][__b])) > 1 ){std::cerr<<"Mismatch detected in MIC implementation"<<std::endl;std::cerr<<"Variable dPdF:"<<std::endl;std::cerr<<"seed="<<seed<<", __a="<<__a<<", __b="<<__b<<std::endl;std::cerr<<"dPdF MIC=  "<<dPdF[__a][__b]<<std::endl;std::cerr<<"dPdF Reference=  "<<dPdF_reference[__a][__b]<<std::endl;std::cerr<<"dPdF Rel Difference=  "<< std::abs((dPdF[__a][__b] - dPdF_reference[__a][__b]) / (dPdF_reference[__a][__b])) << std::endl;std::cerr<<"dPdF Abs Difference=  "<< std::abs(dPdF[__a][__b] - dPdF_reference[__a][__b]) << std::endl;return 1;}
for(int __a=0;__a<3;__a++) for(int __b=0;__b<8;__b++) for(int __c=0;__c<16;__c++) if(std::abs((d[__a][__b][__c] - d_reference[__a][__b][__c]) / (d_reference[__a][__b][__c])) > 1 ){std::cerr<<"Mismatch detected in MIC implementation"<<std::endl;std::cerr<<"Variable d:"<<std::endl;std::cerr<<"seed="<<seed<<", __a="<<__a<<", __b="<<__b<<", __c="<<__c<<std::endl;std::cerr<<"d MIC=  "<<d[__a][__b][__c]<<std::endl;std::cerr<<"d Reference=  "<<d_reference[__a][__b][__c]<<std::endl;std::cerr<<"d Rel Difference=  "<< std::abs((d[__a][__b][__c] - d_reference[__a][__b][__c]) / (d_reference[__a][__b][__c])) << std::endl;std::cerr<<"d Abs Difference=  "<< std::abs(d[__a][__b][__c] - d_reference[__a][__b][__c]) << std::endl;return 1;}
for(int __a=0;__a<300;__a++) for(int __b=0;__b<16;__b++) if(std::abs((system_matrix[__a][__b] - system_matrix_reference[__a][__b]) / (system_matrix_reference[__a][__b])) > 1 ){std::cerr<<"Mismatch detected in MIC implementation"<<std::endl;std::cerr<<"Variable system_matrix:"<<std::endl;std::cerr<<"seed="<<seed<<", __a="<<__a<<", __b="<<__b<<std::endl;std::cerr<<"system_matrix MIC=  "<<system_matrix[__a][__b]<<std::endl;std::cerr<<"system_matrix Reference=  "<<system_matrix_reference[__a][__b]<<std::endl;std::cerr<<"system_matrix Rel Difference=  "<< std::abs((system_matrix[__a][__b] - system_matrix_reference[__a][__b]) / (system_matrix_reference[__a][__b])) << std::endl;std::cerr<<"system_matrix Abs Difference=  "<< std::abs(system_matrix[__a][__b] - system_matrix_reference[__a][__b]) << std::endl;return 1;}

}
#endif

}



{
std::cout << "Running SIMD Test for Update_Position_Based_State " << std::endl;


//=======================================================
//
//        DEFINE ALL VARIABLES USED BY KERNEL
//
//=======================================================

T u[3][8][16] __attribute__((aligned(64)));
T p[16] __attribute__((aligned(64)));
T mu[16] __attribute__((aligned(64)));
T mu_stab[16] __attribute__((aligned(64)));
T kappa[16] __attribute__((aligned(64)));
T alpha[16] __attribute__((aligned(64)));
T cutoff[16] __attribute__((aligned(64)));
T one_over_h[16] __attribute__((aligned(64)));
T cell_volume[16] __attribute__((aligned(64)));
T U[9][16] __attribute__((aligned(64)));
T U_reference[9][16] __attribute__((aligned(64)));
T U_original[9][16] __attribute__((aligned(64)));
T V[9][16] __attribute__((aligned(64)));
T V_reference[9][16] __attribute__((aligned(64)));
T V_original[9][16] __attribute__((aligned(64)));
T Sigma[3][16] __attribute__((aligned(64)));
T Sigma_reference[3][16] __attribute__((aligned(64)));
T Sigma_original[3][16] __attribute__((aligned(64)));
T Q_hat[3][16] __attribute__((aligned(64)));
T Q_hat_reference[3][16] __attribute__((aligned(64)));
T Q_hat_original[3][16] __attribute__((aligned(64)));
T dPdF[12][16] __attribute__((aligned(64)));
T dPdF_reference[12][16] __attribute__((aligned(64)));
T dPdF_original[12][16] __attribute__((aligned(64)));
T d[3][8][16] __attribute__((aligned(64)));
T d_reference[3][8][16] __attribute__((aligned(64)));
T d_original[3][8][16] __attribute__((aligned(64)));
T system_matrix[300][16] __attribute__((aligned(64)));
T system_matrix_reference[300][16] __attribute__((aligned(64)));
T system_matrix_original[300][16] __attribute__((aligned(64)));
  

for(int __a=0;__a<3;__a++) for(int __b=0;__b<8;__b++) for(int __c=0;__c<16;__c++) 
u[__a][__b][__c]=Get_Random<float>();for(int __a=0;__a<16;__a++) 
p[__a]=Get_Random<float>();for(int __a=0;__a<16;__a++) 
mu[__a]=Get_Random<float>();for(int __a=0;__a<16;__a++) 
mu_stab[__a]=Get_Random<float>();for(int __a=0;__a<16;__a++) 
kappa[__a]=Get_Random<float>();for(int __a=0;__a<16;__a++) 
alpha[__a]=Get_Random<float>();for(int __a=0;__a<16;__a++) 
cutoff[__a]=Get_Random<float>();for(int __a=0;__a<16;__a++) 
one_over_h[__a]=Get_Random<float>();for(int __a=0;__a<16;__a++) 
cell_volume[__a]=Get_Random<float>();for(int __a=0;__a<9;__a++) for(int __b=0;__b<16;__b++) {
U_original[__a][__b]=Get_Random<float>();
U[__a][__b]=U_original[__a][__b];
U_reference[__a][__b]=U_original[__a][__b];}
for(int __a=0;__a<9;__a++) for(int __b=0;__b<16;__b++) {
V_original[__a][__b]=Get_Random<float>();
V[__a][__b]=V_original[__a][__b];
V_reference[__a][__b]=V_original[__a][__b];}
for(int __a=0;__a<3;__a++) for(int __b=0;__b<16;__b++) {
Sigma_original[__a][__b]=Get_Random<float>();
Sigma[__a][__b]=Sigma_original[__a][__b];
Sigma_reference[__a][__b]=Sigma_original[__a][__b];}
for(int __a=0;__a<3;__a++) for(int __b=0;__b<16;__b++) {
Q_hat_original[__a][__b]=Get_Random<float>();
Q_hat[__a][__b]=Q_hat_original[__a][__b];
Q_hat_reference[__a][__b]=Q_hat_original[__a][__b];}
for(int __a=0;__a<12;__a++) for(int __b=0;__b<16;__b++) {
dPdF_original[__a][__b]=Get_Random<float>();
dPdF[__a][__b]=dPdF_original[__a][__b];
dPdF_reference[__a][__b]=dPdF_original[__a][__b];}
for(int __a=0;__a<3;__a++) for(int __b=0;__b<8;__b++) for(int __c=0;__c<16;__c++) {
d_original[__a][__b][__c]=Get_Random<float>();
d[__a][__b][__c]=d_original[__a][__b][__c];
d_reference[__a][__b][__c]=d_original[__a][__b][__c];}
for(int __a=0;__a<300;__a++) for(int __b=0;__b<16;__b++) {
system_matrix_original[__a][__b]=Get_Random<float>();
system_matrix[__a][__b]=system_matrix_original[__a][__b];
system_matrix_reference[__a][__b]=system_matrix_original[__a][__b];}


//=======================================================
//
//             COMPUTE REFERENCE RESULTS
//
//=======================================================
 
T __mu[3][8] __attribute__((aligned(4)));
T __mp __attribute__((aligned(4)));
T __mmu __attribute__((aligned(4)));
T __mmu_stab __attribute__((aligned(4)));
T __mkappa __attribute__((aligned(4)));
T __malpha __attribute__((aligned(4)));
T __mcutoff __attribute__((aligned(4)));
T __mone_over_h __attribute__((aligned(4)));
T __mcell_volume __attribute__((aligned(4)));
T __mU[9] __attribute__((aligned(4)));
T __mU_reference[9] __attribute__((aligned(4)));
T __mU_original[9] __attribute__((aligned(4)));
T __mV[9] __attribute__((aligned(4)));
T __mV_reference[9] __attribute__((aligned(4)));
T __mV_original[9] __attribute__((aligned(4)));
T __mSigma[3] __attribute__((aligned(4)));
T __mSigma_reference[3] __attribute__((aligned(4)));
T __mSigma_original[3] __attribute__((aligned(4)));
T __mQ_hat[3] __attribute__((aligned(4)));
T __mQ_hat_reference[3] __attribute__((aligned(4)));
T __mQ_hat_original[3] __attribute__((aligned(4)));
T __mdPdF[12] __attribute__((aligned(4)));
T __mdPdF_reference[12] __attribute__((aligned(4)));
T __mdPdF_original[12] __attribute__((aligned(4)));
T __md[3][8] __attribute__((aligned(4)));
T __md_reference[3][8] __attribute__((aligned(4)));
T __md_original[3][8] __attribute__((aligned(4)));
T __msystem_matrix[300] __attribute__((aligned(4)));
T __msystem_matrix_reference[300] __attribute__((aligned(4)));
T __msystem_matrix_original[300] __attribute__((aligned(4)));
for( int k=0;k<16;k++){for(int __a=0;__a<3;__a++) for(int __b=0;__b<8;__b++) 
__mu[__a][__b]=u[__a][__b][k];
__mp=p[k];
__mmu=mu[k];
__mmu_stab=mu_stab[k];
__mkappa=kappa[k];
__malpha=alpha[k];
__mcutoff=cutoff[k];
__mone_over_h=one_over_h[k];
__mcell_volume=cell_volume[k];for(int __a=0;__a<9;__a++) 
__mU_reference[__a]=U_reference[__a][k];
for(int __a=0;__a<9;__a++) 
__mV_reference[__a]=V_reference[__a][k];
for(int __a=0;__a<3;__a++) 
__mSigma_reference[__a]=Sigma_reference[__a][k];
for(int __a=0;__a<3;__a++) 
__mQ_hat_reference[__a]=Q_hat_reference[__a][k];
for(int __a=0;__a<12;__a++) 
__mdPdF_reference[__a]=dPdF_reference[__a][k];
for(int __a=0;__a<3;__a++) for(int __b=0;__b<8;__b++) 
__md_reference[__a][__b]=d_reference[__a][__b][k];
for(int __a=0;__a<300;__a++) 
__msystem_matrix_reference[__a]=system_matrix_reference[__a][k];
Update_Position_Based_State<COROTATED_TAG,float,float,int>::Run(__mu,__mp,__mmu,__mmu_stab,__mkappa,__malpha,__mcutoff,__mone_over_h,__mcell_volume,__mU_reference,__mV_reference,__mSigma_reference,__mQ_hat_reference,__mdPdF_reference,__md_reference,__msystem_matrix_reference);for(int __a=0;__a<9;__a++) U_reference[__a][k]=
__mU_reference[__a];
for(int __a=0;__a<9;__a++) V_reference[__a][k]=
__mV_reference[__a];
for(int __a=0;__a<3;__a++) Sigma_reference[__a][k]=
__mSigma_reference[__a];
for(int __a=0;__a<3;__a++) Q_hat_reference[__a][k]=
__mQ_hat_reference[__a];
for(int __a=0;__a<12;__a++) dPdF_reference[__a][k]=
__mdPdF_reference[__a];
for(int __a=0;__a<3;__a++) for(int __b=0;__b<8;__b++) d_reference[__a][__b][k]=
__md_reference[__a][__b];
for(int __a=0;__a<300;__a++) system_matrix_reference[__a][k]=
__msystem_matrix_reference[__a];
}

//=======================================================
//
//               COMPUTE SCALAR RESULTS
//
//=======================================================

{
typedef T (&refArray1)[3][8][16];typedef T (&refArray2)[16];typedef T (&refArray3)[16];typedef T (&refArray4)[16];typedef T (&refArray5)[16];typedef T (&refArray6)[16];typedef T (&refArray7)[16];typedef T (&refArray8)[16];typedef T (&refArray9)[16];typedef T (&refArray10)[9][16];typedef T (&refArray11)[9][16];typedef T (&refArray12)[3][16];typedef T (&refArray13)[3][16];typedef T (&refArray14)[12][16];typedef T (&refArray15)[3][8][16];typedef T (&refArray16)[300][16];for(int __a=0;__a<9;__a++) for(int __b=0;__b<16;__b++) U[__a][__b]=U_original[__a][__b];for(int __a=0;__a<9;__a++) for(int __b=0;__b<16;__b++) V[__a][__b]=V_original[__a][__b];for(int __a=0;__a<3;__a++) for(int __b=0;__b<16;__b++) Sigma[__a][__b]=Sigma_original[__a][__b];for(int __a=0;__a<3;__a++) for(int __b=0;__b<16;__b++) Q_hat[__a][__b]=Q_hat_original[__a][__b];for(int __a=0;__a<12;__a++) for(int __b=0;__b<16;__b++) dPdF[__a][__b]=dPdF_original[__a][__b];for(int __a=0;__a<3;__a++) for(int __b=0;__b<8;__b++) for(int __c=0;__c<16;__c++) d[__a][__b][__c]=d_original[__a][__b][__c];for(int __a=0;__a<300;__a++) for(int __b=0;__b<16;__b++) system_matrix[__a][__b]=system_matrix_original[__a][__b];for(int i=0; i<16; i+=1) {refArray1 uk=reinterpret_cast<refArray1>(u[0][0][i]);refArray2 pk=reinterpret_cast<refArray2>(p[i]);refArray3 muk=reinterpret_cast<refArray3>(mu[i]);refArray4 mu_stabk=reinterpret_cast<refArray4>(mu_stab[i]);refArray5 kappak=reinterpret_cast<refArray5>(kappa[i]);refArray6 alphak=reinterpret_cast<refArray6>(alpha[i]);refArray7 cutoffk=reinterpret_cast<refArray7>(cutoff[i]);refArray8 one_over_hk=reinterpret_cast<refArray8>(one_over_h[i]);refArray9 cell_volumek=reinterpret_cast<refArray9>(cell_volume[i]);refArray10 Uk=reinterpret_cast<refArray10>(U[0][i]);refArray11 Vk=reinterpret_cast<refArray11>(V[0][i]);refArray12 Sigmak=reinterpret_cast<refArray12>(Sigma[0][i]);refArray13 Q_hatk=reinterpret_cast<refArray13>(Q_hat[0][i]);refArray14 dPdFk=reinterpret_cast<refArray14>(dPdF[0][i]);refArray15 dk=reinterpret_cast<refArray15>(d[0][0][i]);refArray16 system_matrixk=reinterpret_cast<refArray16>(system_matrix[0][i]);Update_Position_Based_State<COROTATED_TAG,float,float[16],int[16]>::Run(uk,pk,muk,mu_stabk,kappak,alphak,cutoffk,one_over_hk,cell_volumek,Uk,Vk,Sigmak,Q_hatk,dPdFk,dk,system_matrixk);}
for(int __a=0;__a<9;__a++) for(int __b=0;__b<16;__b++) if(std::abs((U[__a][__b] - U_reference[__a][__b]) / (U_reference[__a][__b])) > 1 ){std::cerr<<"Mismatch detected in SCALAR implementation"<<std::endl;std::cerr<<"Variable U:"<<std::endl;std::cerr<<"seed="<<seed<<", __a="<<__a<<", __b="<<__b<<std::endl;std::cerr<<"U SCALAR=  "<<U[__a][__b]<<std::endl;std::cerr<<"U Reference=  "<<U_reference[__a][__b]<<std::endl;std::cerr<<"U Rel Difference=  "<< std::abs((U[__a][__b] - U_reference[__a][__b]) / (U_reference[__a][__b])) << std::endl;std::cerr<<"U Abs Difference=  "<< std::abs(U[__a][__b] - U_reference[__a][__b]) << std::endl;return 1;}
for(int __a=0;__a<9;__a++) for(int __b=0;__b<16;__b++) if(std::abs((V[__a][__b] - V_reference[__a][__b]) / (V_reference[__a][__b])) > 1 ){std::cerr<<"Mismatch detected in SCALAR implementation"<<std::endl;std::cerr<<"Variable V:"<<std::endl;std::cerr<<"seed="<<seed<<", __a="<<__a<<", __b="<<__b<<std::endl;std::cerr<<"V SCALAR=  "<<V[__a][__b]<<std::endl;std::cerr<<"V Reference=  "<<V_reference[__a][__b]<<std::endl;std::cerr<<"V Rel Difference=  "<< std::abs((V[__a][__b] - V_reference[__a][__b]) / (V_reference[__a][__b])) << std::endl;std::cerr<<"V Abs Difference=  "<< std::abs(V[__a][__b] - V_reference[__a][__b]) << std::endl;return 1;}
for(int __a=0;__a<3;__a++) for(int __b=0;__b<16;__b++) if(std::abs((Sigma[__a][__b] - Sigma_reference[__a][__b]) / (Sigma_reference[__a][__b])) > 1 ){std::cerr<<"Mismatch detected in SCALAR implementation"<<std::endl;std::cerr<<"Variable Sigma:"<<std::endl;std::cerr<<"seed="<<seed<<", __a="<<__a<<", __b="<<__b<<std::endl;std::cerr<<"Sigma SCALAR=  "<<Sigma[__a][__b]<<std::endl;std::cerr<<"Sigma Reference=  "<<Sigma_reference[__a][__b]<<std::endl;std::cerr<<"Sigma Rel Difference=  "<< std::abs((Sigma[__a][__b] - Sigma_reference[__a][__b]) / (Sigma_reference[__a][__b])) << std::endl;std::cerr<<"Sigma Abs Difference=  "<< std::abs(Sigma[__a][__b] - Sigma_reference[__a][__b]) << std::endl;return 1;}
for(int __a=0;__a<3;__a++) for(int __b=0;__b<16;__b++) if(std::abs((Q_hat[__a][__b] - Q_hat_reference[__a][__b]) / (Q_hat_reference[__a][__b])) > 1 ){std::cerr<<"Mismatch detected in SCALAR implementation"<<std::endl;std::cerr<<"Variable Q_hat:"<<std::endl;std::cerr<<"seed="<<seed<<", __a="<<__a<<", __b="<<__b<<std::endl;std::cerr<<"Q_hat SCALAR=  "<<Q_hat[__a][__b]<<std::endl;std::cerr<<"Q_hat Reference=  "<<Q_hat_reference[__a][__b]<<std::endl;std::cerr<<"Q_hat Rel Difference=  "<< std::abs((Q_hat[__a][__b] - Q_hat_reference[__a][__b]) / (Q_hat_reference[__a][__b])) << std::endl;std::cerr<<"Q_hat Abs Difference=  "<< std::abs(Q_hat[__a][__b] - Q_hat_reference[__a][__b]) << std::endl;return 1;}
for(int __a=0;__a<12;__a++) for(int __b=0;__b<16;__b++) if(std::abs((dPdF[__a][__b] - dPdF_reference[__a][__b]) / (dPdF_reference[__a][__b])) > 1 ){std::cerr<<"Mismatch detected in SCALAR implementation"<<std::endl;std::cerr<<"Variable dPdF:"<<std::endl;std::cerr<<"seed="<<seed<<", __a="<<__a<<", __b="<<__b<<std::endl;std::cerr<<"dPdF SCALAR=  "<<dPdF[__a][__b]<<std::endl;std::cerr<<"dPdF Reference=  "<<dPdF_reference[__a][__b]<<std::endl;std::cerr<<"dPdF Rel Difference=  "<< std::abs((dPdF[__a][__b] - dPdF_reference[__a][__b]) / (dPdF_reference[__a][__b])) << std::endl;std::cerr<<"dPdF Abs Difference=  "<< std::abs(dPdF[__a][__b] - dPdF_reference[__a][__b]) << std::endl;return 1;}
for(int __a=0;__a<3;__a++) for(int __b=0;__b<8;__b++) for(int __c=0;__c<16;__c++) if(std::abs((d[__a][__b][__c] - d_reference[__a][__b][__c]) / (d_reference[__a][__b][__c])) > 1 ){std::cerr<<"Mismatch detected in SCALAR implementation"<<std::endl;std::cerr<<"Variable d:"<<std::endl;std::cerr<<"seed="<<seed<<", __a="<<__a<<", __b="<<__b<<", __c="<<__c<<std::endl;std::cerr<<"d SCALAR=  "<<d[__a][__b][__c]<<std::endl;std::cerr<<"d Reference=  "<<d_reference[__a][__b][__c]<<std::endl;std::cerr<<"d Rel Difference=  "<< std::abs((d[__a][__b][__c] - d_reference[__a][__b][__c]) / (d_reference[__a][__b][__c])) << std::endl;std::cerr<<"d Abs Difference=  "<< std::abs(d[__a][__b][__c] - d_reference[__a][__b][__c]) << std::endl;return 1;}
for(int __a=0;__a<300;__a++) for(int __b=0;__b<16;__b++) if(std::abs((system_matrix[__a][__b] - system_matrix_reference[__a][__b]) / (system_matrix_reference[__a][__b])) > 1 ){std::cerr<<"Mismatch detected in SCALAR implementation"<<std::endl;std::cerr<<"Variable system_matrix:"<<std::endl;std::cerr<<"seed="<<seed<<", __a="<<__a<<", __b="<<__b<<std::endl;std::cerr<<"system_matrix SCALAR=  "<<system_matrix[__a][__b]<<std::endl;std::cerr<<"system_matrix Reference=  "<<system_matrix_reference[__a][__b]<<std::endl;std::cerr<<"system_matrix Rel Difference=  "<< std::abs((system_matrix[__a][__b] - system_matrix_reference[__a][__b]) / (system_matrix_reference[__a][__b])) << std::endl;std::cerr<<"system_matrix Abs Difference=  "<< std::abs(system_matrix[__a][__b] - system_matrix_reference[__a][__b]) << std::endl;return 1;}

}

//=======================================================
//
//               COMPUTE SSE RESULTS
//
//=======================================================

#ifdef ENABLE_SSE_INSTRUCTION_SET
{
typedef T (&refArray1)[3][8][16];typedef T (&refArray2)[16];typedef T (&refArray3)[16];typedef T (&refArray4)[16];typedef T (&refArray5)[16];typedef T (&refArray6)[16];typedef T (&refArray7)[16];typedef T (&refArray8)[16];typedef T (&refArray9)[16];typedef T (&refArray10)[9][16];typedef T (&refArray11)[9][16];typedef T (&refArray12)[3][16];typedef T (&refArray13)[3][16];typedef T (&refArray14)[12][16];typedef T (&refArray15)[3][8][16];typedef T (&refArray16)[300][16];for(int __a=0;__a<9;__a++) for(int __b=0;__b<16;__b++) U[__a][__b]=U_original[__a][__b];for(int __a=0;__a<9;__a++) for(int __b=0;__b<16;__b++) V[__a][__b]=V_original[__a][__b];for(int __a=0;__a<3;__a++) for(int __b=0;__b<16;__b++) Sigma[__a][__b]=Sigma_original[__a][__b];for(int __a=0;__a<3;__a++) for(int __b=0;__b<16;__b++) Q_hat[__a][__b]=Q_hat_original[__a][__b];for(int __a=0;__a<12;__a++) for(int __b=0;__b<16;__b++) dPdF[__a][__b]=dPdF_original[__a][__b];for(int __a=0;__a<3;__a++) for(int __b=0;__b<8;__b++) for(int __c=0;__c<16;__c++) d[__a][__b][__c]=d_original[__a][__b][__c];for(int __a=0;__a<300;__a++) for(int __b=0;__b<16;__b++) system_matrix[__a][__b]=system_matrix_original[__a][__b];for(int i=0; i<16; i+=4) {refArray1 uk=reinterpret_cast<refArray1>(u[0][0][i]);refArray2 pk=reinterpret_cast<refArray2>(p[i]);refArray3 muk=reinterpret_cast<refArray3>(mu[i]);refArray4 mu_stabk=reinterpret_cast<refArray4>(mu_stab[i]);refArray5 kappak=reinterpret_cast<refArray5>(kappa[i]);refArray6 alphak=reinterpret_cast<refArray6>(alpha[i]);refArray7 cutoffk=reinterpret_cast<refArray7>(cutoff[i]);refArray8 one_over_hk=reinterpret_cast<refArray8>(one_over_h[i]);refArray9 cell_volumek=reinterpret_cast<refArray9>(cell_volume[i]);refArray10 Uk=reinterpret_cast<refArray10>(U[0][i]);refArray11 Vk=reinterpret_cast<refArray11>(V[0][i]);refArray12 Sigmak=reinterpret_cast<refArray12>(Sigma[0][i]);refArray13 Q_hatk=reinterpret_cast<refArray13>(Q_hat[0][i]);refArray14 dPdFk=reinterpret_cast<refArray14>(dPdF[0][i]);refArray15 dk=reinterpret_cast<refArray15>(d[0][0][i]);refArray16 system_matrixk=reinterpret_cast<refArray16>(system_matrix[0][i]);Update_Position_Based_State<COROTATED_TAG,__m128,float[16],int[16]>::Run(uk,pk,muk,mu_stabk,kappak,alphak,cutoffk,one_over_hk,cell_volumek,Uk,Vk,Sigmak,Q_hatk,dPdFk,dk,system_matrixk);}
for(int __a=0;__a<9;__a++) for(int __b=0;__b<16;__b++) if(std::abs((U[__a][__b] - U_reference[__a][__b]) / (U_reference[__a][__b])) > 1 ){std::cerr<<"Mismatch detected in SSE implementation"<<std::endl;std::cerr<<"Variable U:"<<std::endl;std::cerr<<"seed="<<seed<<", __a="<<__a<<", __b="<<__b<<std::endl;std::cerr<<"U SSE=  "<<U[__a][__b]<<std::endl;std::cerr<<"U Reference=  "<<U_reference[__a][__b]<<std::endl;std::cerr<<"U Rel Difference=  "<< std::abs((U[__a][__b] - U_reference[__a][__b]) / (U_reference[__a][__b])) << std::endl;std::cerr<<"U Abs Difference=  "<< std::abs(U[__a][__b] - U_reference[__a][__b]) << std::endl;return 1;}
for(int __a=0;__a<9;__a++) for(int __b=0;__b<16;__b++) if(std::abs((V[__a][__b] - V_reference[__a][__b]) / (V_reference[__a][__b])) > 1 ){std::cerr<<"Mismatch detected in SSE implementation"<<std::endl;std::cerr<<"Variable V:"<<std::endl;std::cerr<<"seed="<<seed<<", __a="<<__a<<", __b="<<__b<<std::endl;std::cerr<<"V SSE=  "<<V[__a][__b]<<std::endl;std::cerr<<"V Reference=  "<<V_reference[__a][__b]<<std::endl;std::cerr<<"V Rel Difference=  "<< std::abs((V[__a][__b] - V_reference[__a][__b]) / (V_reference[__a][__b])) << std::endl;std::cerr<<"V Abs Difference=  "<< std::abs(V[__a][__b] - V_reference[__a][__b]) << std::endl;return 1;}
for(int __a=0;__a<3;__a++) for(int __b=0;__b<16;__b++) if(std::abs((Sigma[__a][__b] - Sigma_reference[__a][__b]) / (Sigma_reference[__a][__b])) > 1 ){std::cerr<<"Mismatch detected in SSE implementation"<<std::endl;std::cerr<<"Variable Sigma:"<<std::endl;std::cerr<<"seed="<<seed<<", __a="<<__a<<", __b="<<__b<<std::endl;std::cerr<<"Sigma SSE=  "<<Sigma[__a][__b]<<std::endl;std::cerr<<"Sigma Reference=  "<<Sigma_reference[__a][__b]<<std::endl;std::cerr<<"Sigma Rel Difference=  "<< std::abs((Sigma[__a][__b] - Sigma_reference[__a][__b]) / (Sigma_reference[__a][__b])) << std::endl;std::cerr<<"Sigma Abs Difference=  "<< std::abs(Sigma[__a][__b] - Sigma_reference[__a][__b]) << std::endl;return 1;}
for(int __a=0;__a<3;__a++) for(int __b=0;__b<16;__b++) if(std::abs((Q_hat[__a][__b] - Q_hat_reference[__a][__b]) / (Q_hat_reference[__a][__b])) > 1 ){std::cerr<<"Mismatch detected in SSE implementation"<<std::endl;std::cerr<<"Variable Q_hat:"<<std::endl;std::cerr<<"seed="<<seed<<", __a="<<__a<<", __b="<<__b<<std::endl;std::cerr<<"Q_hat SSE=  "<<Q_hat[__a][__b]<<std::endl;std::cerr<<"Q_hat Reference=  "<<Q_hat_reference[__a][__b]<<std::endl;std::cerr<<"Q_hat Rel Difference=  "<< std::abs((Q_hat[__a][__b] - Q_hat_reference[__a][__b]) / (Q_hat_reference[__a][__b])) << std::endl;std::cerr<<"Q_hat Abs Difference=  "<< std::abs(Q_hat[__a][__b] - Q_hat_reference[__a][__b]) << std::endl;return 1;}
for(int __a=0;__a<12;__a++) for(int __b=0;__b<16;__b++) if(std::abs((dPdF[__a][__b] - dPdF_reference[__a][__b]) / (dPdF_reference[__a][__b])) > 1 ){std::cerr<<"Mismatch detected in SSE implementation"<<std::endl;std::cerr<<"Variable dPdF:"<<std::endl;std::cerr<<"seed="<<seed<<", __a="<<__a<<", __b="<<__b<<std::endl;std::cerr<<"dPdF SSE=  "<<dPdF[__a][__b]<<std::endl;std::cerr<<"dPdF Reference=  "<<dPdF_reference[__a][__b]<<std::endl;std::cerr<<"dPdF Rel Difference=  "<< std::abs((dPdF[__a][__b] - dPdF_reference[__a][__b]) / (dPdF_reference[__a][__b])) << std::endl;std::cerr<<"dPdF Abs Difference=  "<< std::abs(dPdF[__a][__b] - dPdF_reference[__a][__b]) << std::endl;return 1;}
for(int __a=0;__a<3;__a++) for(int __b=0;__b<8;__b++) for(int __c=0;__c<16;__c++) if(std::abs((d[__a][__b][__c] - d_reference[__a][__b][__c]) / (d_reference[__a][__b][__c])) > 1 ){std::cerr<<"Mismatch detected in SSE implementation"<<std::endl;std::cerr<<"Variable d:"<<std::endl;std::cerr<<"seed="<<seed<<", __a="<<__a<<", __b="<<__b<<", __c="<<__c<<std::endl;std::cerr<<"d SSE=  "<<d[__a][__b][__c]<<std::endl;std::cerr<<"d Reference=  "<<d_reference[__a][__b][__c]<<std::endl;std::cerr<<"d Rel Difference=  "<< std::abs((d[__a][__b][__c] - d_reference[__a][__b][__c]) / (d_reference[__a][__b][__c])) << std::endl;std::cerr<<"d Abs Difference=  "<< std::abs(d[__a][__b][__c] - d_reference[__a][__b][__c]) << std::endl;return 1;}
for(int __a=0;__a<300;__a++) for(int __b=0;__b<16;__b++) if(std::abs((system_matrix[__a][__b] - system_matrix_reference[__a][__b]) / (system_matrix_reference[__a][__b])) > 1 ){std::cerr<<"Mismatch detected in SSE implementation"<<std::endl;std::cerr<<"Variable system_matrix:"<<std::endl;std::cerr<<"seed="<<seed<<", __a="<<__a<<", __b="<<__b<<std::endl;std::cerr<<"system_matrix SSE=  "<<system_matrix[__a][__b]<<std::endl;std::cerr<<"system_matrix Reference=  "<<system_matrix_reference[__a][__b]<<std::endl;std::cerr<<"system_matrix Rel Difference=  "<< std::abs((system_matrix[__a][__b] - system_matrix_reference[__a][__b]) / (system_matrix_reference[__a][__b])) << std::endl;std::cerr<<"system_matrix Abs Difference=  "<< std::abs(system_matrix[__a][__b] - system_matrix_reference[__a][__b]) << std::endl;return 1;}

}
#endif

//=======================================================
//
//               COMPUTE AVX RESULTS
//
//=======================================================

#ifdef ENABLE_AVX_INSTRUCTION_SET
{
typedef T (&refArray1)[3][8][16];typedef T (&refArray2)[16];typedef T (&refArray3)[16];typedef T (&refArray4)[16];typedef T (&refArray5)[16];typedef T (&refArray6)[16];typedef T (&refArray7)[16];typedef T (&refArray8)[16];typedef T (&refArray9)[16];typedef T (&refArray10)[9][16];typedef T (&refArray11)[9][16];typedef T (&refArray12)[3][16];typedef T (&refArray13)[3][16];typedef T (&refArray14)[12][16];typedef T (&refArray15)[3][8][16];typedef T (&refArray16)[300][16];for(int __a=0;__a<9;__a++) for(int __b=0;__b<16;__b++) U[__a][__b]=U_original[__a][__b];for(int __a=0;__a<9;__a++) for(int __b=0;__b<16;__b++) V[__a][__b]=V_original[__a][__b];for(int __a=0;__a<3;__a++) for(int __b=0;__b<16;__b++) Sigma[__a][__b]=Sigma_original[__a][__b];for(int __a=0;__a<3;__a++) for(int __b=0;__b<16;__b++) Q_hat[__a][__b]=Q_hat_original[__a][__b];for(int __a=0;__a<12;__a++) for(int __b=0;__b<16;__b++) dPdF[__a][__b]=dPdF_original[__a][__b];for(int __a=0;__a<3;__a++) for(int __b=0;__b<8;__b++) for(int __c=0;__c<16;__c++) d[__a][__b][__c]=d_original[__a][__b][__c];for(int __a=0;__a<300;__a++) for(int __b=0;__b<16;__b++) system_matrix[__a][__b]=system_matrix_original[__a][__b];for(int i=0; i<16; i+=8) {refArray1 uk=reinterpret_cast<refArray1>(u[0][0][i]);refArray2 pk=reinterpret_cast<refArray2>(p[i]);refArray3 muk=reinterpret_cast<refArray3>(mu[i]);refArray4 mu_stabk=reinterpret_cast<refArray4>(mu_stab[i]);refArray5 kappak=reinterpret_cast<refArray5>(kappa[i]);refArray6 alphak=reinterpret_cast<refArray6>(alpha[i]);refArray7 cutoffk=reinterpret_cast<refArray7>(cutoff[i]);refArray8 one_over_hk=reinterpret_cast<refArray8>(one_over_h[i]);refArray9 cell_volumek=reinterpret_cast<refArray9>(cell_volume[i]);refArray10 Uk=reinterpret_cast<refArray10>(U[0][i]);refArray11 Vk=reinterpret_cast<refArray11>(V[0][i]);refArray12 Sigmak=reinterpret_cast<refArray12>(Sigma[0][i]);refArray13 Q_hatk=reinterpret_cast<refArray13>(Q_hat[0][i]);refArray14 dPdFk=reinterpret_cast<refArray14>(dPdF[0][i]);refArray15 dk=reinterpret_cast<refArray15>(d[0][0][i]);refArray16 system_matrixk=reinterpret_cast<refArray16>(system_matrix[0][i]);Update_Position_Based_State<COROTATED_TAG,__m256,float[16],int[16]>::Run(uk,pk,muk,mu_stabk,kappak,alphak,cutoffk,one_over_hk,cell_volumek,Uk,Vk,Sigmak,Q_hatk,dPdFk,dk,system_matrixk);}
for(int __a=0;__a<9;__a++) for(int __b=0;__b<16;__b++) if(std::abs((U[__a][__b] - U_reference[__a][__b]) / (U_reference[__a][__b])) > 1 ){std::cerr<<"Mismatch detected in AVX implementation"<<std::endl;std::cerr<<"Variable U:"<<std::endl;std::cerr<<"seed="<<seed<<", __a="<<__a<<", __b="<<__b<<std::endl;std::cerr<<"U AVX=  "<<U[__a][__b]<<std::endl;std::cerr<<"U Reference=  "<<U_reference[__a][__b]<<std::endl;std::cerr<<"U Rel Difference=  "<< std::abs((U[__a][__b] - U_reference[__a][__b]) / (U_reference[__a][__b])) << std::endl;std::cerr<<"U Abs Difference=  "<< std::abs(U[__a][__b] - U_reference[__a][__b]) << std::endl;return 1;}
for(int __a=0;__a<9;__a++) for(int __b=0;__b<16;__b++) if(std::abs((V[__a][__b] - V_reference[__a][__b]) / (V_reference[__a][__b])) > 1 ){std::cerr<<"Mismatch detected in AVX implementation"<<std::endl;std::cerr<<"Variable V:"<<std::endl;std::cerr<<"seed="<<seed<<", __a="<<__a<<", __b="<<__b<<std::endl;std::cerr<<"V AVX=  "<<V[__a][__b]<<std::endl;std::cerr<<"V Reference=  "<<V_reference[__a][__b]<<std::endl;std::cerr<<"V Rel Difference=  "<< std::abs((V[__a][__b] - V_reference[__a][__b]) / (V_reference[__a][__b])) << std::endl;std::cerr<<"V Abs Difference=  "<< std::abs(V[__a][__b] - V_reference[__a][__b]) << std::endl;return 1;}
for(int __a=0;__a<3;__a++) for(int __b=0;__b<16;__b++) if(std::abs((Sigma[__a][__b] - Sigma_reference[__a][__b]) / (Sigma_reference[__a][__b])) > 1 ){std::cerr<<"Mismatch detected in AVX implementation"<<std::endl;std::cerr<<"Variable Sigma:"<<std::endl;std::cerr<<"seed="<<seed<<", __a="<<__a<<", __b="<<__b<<std::endl;std::cerr<<"Sigma AVX=  "<<Sigma[__a][__b]<<std::endl;std::cerr<<"Sigma Reference=  "<<Sigma_reference[__a][__b]<<std::endl;std::cerr<<"Sigma Rel Difference=  "<< std::abs((Sigma[__a][__b] - Sigma_reference[__a][__b]) / (Sigma_reference[__a][__b])) << std::endl;std::cerr<<"Sigma Abs Difference=  "<< std::abs(Sigma[__a][__b] - Sigma_reference[__a][__b]) << std::endl;return 1;}
for(int __a=0;__a<3;__a++) for(int __b=0;__b<16;__b++) if(std::abs((Q_hat[__a][__b] - Q_hat_reference[__a][__b]) / (Q_hat_reference[__a][__b])) > 1 ){std::cerr<<"Mismatch detected in AVX implementation"<<std::endl;std::cerr<<"Variable Q_hat:"<<std::endl;std::cerr<<"seed="<<seed<<", __a="<<__a<<", __b="<<__b<<std::endl;std::cerr<<"Q_hat AVX=  "<<Q_hat[__a][__b]<<std::endl;std::cerr<<"Q_hat Reference=  "<<Q_hat_reference[__a][__b]<<std::endl;std::cerr<<"Q_hat Rel Difference=  "<< std::abs((Q_hat[__a][__b] - Q_hat_reference[__a][__b]) / (Q_hat_reference[__a][__b])) << std::endl;std::cerr<<"Q_hat Abs Difference=  "<< std::abs(Q_hat[__a][__b] - Q_hat_reference[__a][__b]) << std::endl;return 1;}
for(int __a=0;__a<12;__a++) for(int __b=0;__b<16;__b++) if(std::abs((dPdF[__a][__b] - dPdF_reference[__a][__b]) / (dPdF_reference[__a][__b])) > 1 ){std::cerr<<"Mismatch detected in AVX implementation"<<std::endl;std::cerr<<"Variable dPdF:"<<std::endl;std::cerr<<"seed="<<seed<<", __a="<<__a<<", __b="<<__b<<std::endl;std::cerr<<"dPdF AVX=  "<<dPdF[__a][__b]<<std::endl;std::cerr<<"dPdF Reference=  "<<dPdF_reference[__a][__b]<<std::endl;std::cerr<<"dPdF Rel Difference=  "<< std::abs((dPdF[__a][__b] - dPdF_reference[__a][__b]) / (dPdF_reference[__a][__b])) << std::endl;std::cerr<<"dPdF Abs Difference=  "<< std::abs(dPdF[__a][__b] - dPdF_reference[__a][__b]) << std::endl;return 1;}
for(int __a=0;__a<3;__a++) for(int __b=0;__b<8;__b++) for(int __c=0;__c<16;__c++) if(std::abs((d[__a][__b][__c] - d_reference[__a][__b][__c]) / (d_reference[__a][__b][__c])) > 1 ){std::cerr<<"Mismatch detected in AVX implementation"<<std::endl;std::cerr<<"Variable d:"<<std::endl;std::cerr<<"seed="<<seed<<", __a="<<__a<<", __b="<<__b<<", __c="<<__c<<std::endl;std::cerr<<"d AVX=  "<<d[__a][__b][__c]<<std::endl;std::cerr<<"d Reference=  "<<d_reference[__a][__b][__c]<<std::endl;std::cerr<<"d Rel Difference=  "<< std::abs((d[__a][__b][__c] - d_reference[__a][__b][__c]) / (d_reference[__a][__b][__c])) << std::endl;std::cerr<<"d Abs Difference=  "<< std::abs(d[__a][__b][__c] - d_reference[__a][__b][__c]) << std::endl;return 1;}
for(int __a=0;__a<300;__a++) for(int __b=0;__b<16;__b++) if(std::abs((system_matrix[__a][__b] - system_matrix_reference[__a][__b]) / (system_matrix_reference[__a][__b])) > 1 ){std::cerr<<"Mismatch detected in AVX implementation"<<std::endl;std::cerr<<"Variable system_matrix:"<<std::endl;std::cerr<<"seed="<<seed<<", __a="<<__a<<", __b="<<__b<<std::endl;std::cerr<<"system_matrix AVX=  "<<system_matrix[__a][__b]<<std::endl;std::cerr<<"system_matrix Reference=  "<<system_matrix_reference[__a][__b]<<std::endl;std::cerr<<"system_matrix Rel Difference=  "<< std::abs((system_matrix[__a][__b] - system_matrix_reference[__a][__b]) / (system_matrix_reference[__a][__b])) << std::endl;std::cerr<<"system_matrix Abs Difference=  "<< std::abs(system_matrix[__a][__b] - system_matrix_reference[__a][__b]) << std::endl;return 1;}

}
#endif

//=======================================================
//
//               COMPUTE NEON RESULTS
//
//=======================================================

#ifdef ENABLE_NEON_INSTRUCTION_SET
{
typedef T (&refArray1)[3][8][16];typedef T (&refArray2)[16];typedef T (&refArray3)[16];typedef T (&refArray4)[16];typedef T (&refArray5)[16];typedef T (&refArray6)[16];typedef T (&refArray7)[16];typedef T (&refArray8)[16];typedef T (&refArray9)[16];typedef T (&refArray10)[9][16];typedef T (&refArray11)[9][16];typedef T (&refArray12)[3][16];typedef T (&refArray13)[3][16];typedef T (&refArray14)[12][16];typedef T (&refArray15)[3][8][16];typedef T (&refArray16)[300][16];for(int __a=0;__a<9;__a++) for(int __b=0;__b<16;__b++) U[__a][__b]=U_original[__a][__b];for(int __a=0;__a<9;__a++) for(int __b=0;__b<16;__b++) V[__a][__b]=V_original[__a][__b];for(int __a=0;__a<3;__a++) for(int __b=0;__b<16;__b++) Sigma[__a][__b]=Sigma_original[__a][__b];for(int __a=0;__a<3;__a++) for(int __b=0;__b<16;__b++) Q_hat[__a][__b]=Q_hat_original[__a][__b];for(int __a=0;__a<12;__a++) for(int __b=0;__b<16;__b++) dPdF[__a][__b]=dPdF_original[__a][__b];for(int __a=0;__a<3;__a++) for(int __b=0;__b<8;__b++) for(int __c=0;__c<16;__c++) d[__a][__b][__c]=d_original[__a][__b][__c];for(int __a=0;__a<300;__a++) for(int __b=0;__b<16;__b++) system_matrix[__a][__b]=system_matrix_original[__a][__b];for(int i=0; i<16; i+=4) {refArray1 uk=reinterpret_cast<refArray1>(u[0][0][i]);refArray2 pk=reinterpret_cast<refArray2>(p[i]);refArray3 muk=reinterpret_cast<refArray3>(mu[i]);refArray4 mu_stabk=reinterpret_cast<refArray4>(mu_stab[i]);refArray5 kappak=reinterpret_cast<refArray5>(kappa[i]);refArray6 alphak=reinterpret_cast<refArray6>(alpha[i]);refArray7 cutoffk=reinterpret_cast<refArray7>(cutoff[i]);refArray8 one_over_hk=reinterpret_cast<refArray8>(one_over_h[i]);refArray9 cell_volumek=reinterpret_cast<refArray9>(cell_volume[i]);refArray10 Uk=reinterpret_cast<refArray10>(U[0][i]);refArray11 Vk=reinterpret_cast<refArray11>(V[0][i]);refArray12 Sigmak=reinterpret_cast<refArray12>(Sigma[0][i]);refArray13 Q_hatk=reinterpret_cast<refArray13>(Q_hat[0][i]);refArray14 dPdFk=reinterpret_cast<refArray14>(dPdF[0][i]);refArray15 dk=reinterpret_cast<refArray15>(d[0][0][i]);refArray16 system_matrixk=reinterpret_cast<refArray16>(system_matrix[0][i]);Update_Position_Based_State<COROTATED_TAG,float32x4_t,float[16],int[16]>::Run(uk,pk,muk,mu_stabk,kappak,alphak,cutoffk,one_over_hk,cell_volumek,Uk,Vk,Sigmak,Q_hatk,dPdFk,dk,system_matrixk);}
for(int __a=0;__a<9;__a++) for(int __b=0;__b<16;__b++) if(std::abs((U[__a][__b] - U_reference[__a][__b]) / (U_reference[__a][__b])) > 1 ){std::cerr<<"Mismatch detected in NEON implementation"<<std::endl;std::cerr<<"Variable U:"<<std::endl;std::cerr<<"seed="<<seed<<", __a="<<__a<<", __b="<<__b<<std::endl;std::cerr<<"U NEON=  "<<U[__a][__b]<<std::endl;std::cerr<<"U Reference=  "<<U_reference[__a][__b]<<std::endl;std::cerr<<"U Rel Difference=  "<< std::abs((U[__a][__b] - U_reference[__a][__b]) / (U_reference[__a][__b])) << std::endl;std::cerr<<"U Abs Difference=  "<< std::abs(U[__a][__b] - U_reference[__a][__b]) << std::endl;return 1;}
for(int __a=0;__a<9;__a++) for(int __b=0;__b<16;__b++) if(std::abs((V[__a][__b] - V_reference[__a][__b]) / (V_reference[__a][__b])) > 1 ){std::cerr<<"Mismatch detected in NEON implementation"<<std::endl;std::cerr<<"Variable V:"<<std::endl;std::cerr<<"seed="<<seed<<", __a="<<__a<<", __b="<<__b<<std::endl;std::cerr<<"V NEON=  "<<V[__a][__b]<<std::endl;std::cerr<<"V Reference=  "<<V_reference[__a][__b]<<std::endl;std::cerr<<"V Rel Difference=  "<< std::abs((V[__a][__b] - V_reference[__a][__b]) / (V_reference[__a][__b])) << std::endl;std::cerr<<"V Abs Difference=  "<< std::abs(V[__a][__b] - V_reference[__a][__b]) << std::endl;return 1;}
for(int __a=0;__a<3;__a++) for(int __b=0;__b<16;__b++) if(std::abs((Sigma[__a][__b] - Sigma_reference[__a][__b]) / (Sigma_reference[__a][__b])) > 1 ){std::cerr<<"Mismatch detected in NEON implementation"<<std::endl;std::cerr<<"Variable Sigma:"<<std::endl;std::cerr<<"seed="<<seed<<", __a="<<__a<<", __b="<<__b<<std::endl;std::cerr<<"Sigma NEON=  "<<Sigma[__a][__b]<<std::endl;std::cerr<<"Sigma Reference=  "<<Sigma_reference[__a][__b]<<std::endl;std::cerr<<"Sigma Rel Difference=  "<< std::abs((Sigma[__a][__b] - Sigma_reference[__a][__b]) / (Sigma_reference[__a][__b])) << std::endl;std::cerr<<"Sigma Abs Difference=  "<< std::abs(Sigma[__a][__b] - Sigma_reference[__a][__b]) << std::endl;return 1;}
for(int __a=0;__a<3;__a++) for(int __b=0;__b<16;__b++) if(std::abs((Q_hat[__a][__b] - Q_hat_reference[__a][__b]) / (Q_hat_reference[__a][__b])) > 1 ){std::cerr<<"Mismatch detected in NEON implementation"<<std::endl;std::cerr<<"Variable Q_hat:"<<std::endl;std::cerr<<"seed="<<seed<<", __a="<<__a<<", __b="<<__b<<std::endl;std::cerr<<"Q_hat NEON=  "<<Q_hat[__a][__b]<<std::endl;std::cerr<<"Q_hat Reference=  "<<Q_hat_reference[__a][__b]<<std::endl;std::cerr<<"Q_hat Rel Difference=  "<< std::abs((Q_hat[__a][__b] - Q_hat_reference[__a][__b]) / (Q_hat_reference[__a][__b])) << std::endl;std::cerr<<"Q_hat Abs Difference=  "<< std::abs(Q_hat[__a][__b] - Q_hat_reference[__a][__b]) << std::endl;return 1;}
for(int __a=0;__a<12;__a++) for(int __b=0;__b<16;__b++) if(std::abs((dPdF[__a][__b] - dPdF_reference[__a][__b]) / (dPdF_reference[__a][__b])) > 1 ){std::cerr<<"Mismatch detected in NEON implementation"<<std::endl;std::cerr<<"Variable dPdF:"<<std::endl;std::cerr<<"seed="<<seed<<", __a="<<__a<<", __b="<<__b<<std::endl;std::cerr<<"dPdF NEON=  "<<dPdF[__a][__b]<<std::endl;std::cerr<<"dPdF Reference=  "<<dPdF_reference[__a][__b]<<std::endl;std::cerr<<"dPdF Rel Difference=  "<< std::abs((dPdF[__a][__b] - dPdF_reference[__a][__b]) / (dPdF_reference[__a][__b])) << std::endl;std::cerr<<"dPdF Abs Difference=  "<< std::abs(dPdF[__a][__b] - dPdF_reference[__a][__b]) << std::endl;return 1;}
for(int __a=0;__a<3;__a++) for(int __b=0;__b<8;__b++) for(int __c=0;__c<16;__c++) if(std::abs((d[__a][__b][__c] - d_reference[__a][__b][__c]) / (d_reference[__a][__b][__c])) > 1 ){std::cerr<<"Mismatch detected in NEON implementation"<<std::endl;std::cerr<<"Variable d:"<<std::endl;std::cerr<<"seed="<<seed<<", __a="<<__a<<", __b="<<__b<<", __c="<<__c<<std::endl;std::cerr<<"d NEON=  "<<d[__a][__b][__c]<<std::endl;std::cerr<<"d Reference=  "<<d_reference[__a][__b][__c]<<std::endl;std::cerr<<"d Rel Difference=  "<< std::abs((d[__a][__b][__c] - d_reference[__a][__b][__c]) / (d_reference[__a][__b][__c])) << std::endl;std::cerr<<"d Abs Difference=  "<< std::abs(d[__a][__b][__c] - d_reference[__a][__b][__c]) << std::endl;return 1;}
for(int __a=0;__a<300;__a++) for(int __b=0;__b<16;__b++) if(std::abs((system_matrix[__a][__b] - system_matrix_reference[__a][__b]) / (system_matrix_reference[__a][__b])) > 1 ){std::cerr<<"Mismatch detected in NEON implementation"<<std::endl;std::cerr<<"Variable system_matrix:"<<std::endl;std::cerr<<"seed="<<seed<<", __a="<<__a<<", __b="<<__b<<std::endl;std::cerr<<"system_matrix NEON=  "<<system_matrix[__a][__b]<<std::endl;std::cerr<<"system_matrix Reference=  "<<system_matrix_reference[__a][__b]<<std::endl;std::cerr<<"system_matrix Rel Difference=  "<< std::abs((system_matrix[__a][__b] - system_matrix_reference[__a][__b]) / (system_matrix_reference[__a][__b])) << std::endl;std::cerr<<"system_matrix Abs Difference=  "<< std::abs(system_matrix[__a][__b] - system_matrix_reference[__a][__b]) << std::endl;return 1;}

}
#endif

//=======================================================
//
//               COMPUTE MIC RESULTS
//
//=======================================================

#ifdef ENABLE_MIC_INSTRUCTION_SET
{
typedef T (&refArray1)[3][8][16];typedef T (&refArray2)[16];typedef T (&refArray3)[16];typedef T (&refArray4)[16];typedef T (&refArray5)[16];typedef T (&refArray6)[16];typedef T (&refArray7)[16];typedef T (&refArray8)[16];typedef T (&refArray9)[16];typedef T (&refArray10)[9][16];typedef T (&refArray11)[9][16];typedef T (&refArray12)[3][16];typedef T (&refArray13)[3][16];typedef T (&refArray14)[12][16];typedef T (&refArray15)[3][8][16];typedef T (&refArray16)[300][16];for(int __a=0;__a<9;__a++) for(int __b=0;__b<16;__b++) U[__a][__b]=U_original[__a][__b];for(int __a=0;__a<9;__a++) for(int __b=0;__b<16;__b++) V[__a][__b]=V_original[__a][__b];for(int __a=0;__a<3;__a++) for(int __b=0;__b<16;__b++) Sigma[__a][__b]=Sigma_original[__a][__b];for(int __a=0;__a<3;__a++) for(int __b=0;__b<16;__b++) Q_hat[__a][__b]=Q_hat_original[__a][__b];for(int __a=0;__a<12;__a++) for(int __b=0;__b<16;__b++) dPdF[__a][__b]=dPdF_original[__a][__b];for(int __a=0;__a<3;__a++) for(int __b=0;__b<8;__b++) for(int __c=0;__c<16;__c++) d[__a][__b][__c]=d_original[__a][__b][__c];for(int __a=0;__a<300;__a++) for(int __b=0;__b<16;__b++) system_matrix[__a][__b]=system_matrix_original[__a][__b];for(int i=0; i<16; i+=16) {refArray1 uk=reinterpret_cast<refArray1>(u[0][0][i]);refArray2 pk=reinterpret_cast<refArray2>(p[i]);refArray3 muk=reinterpret_cast<refArray3>(mu[i]);refArray4 mu_stabk=reinterpret_cast<refArray4>(mu_stab[i]);refArray5 kappak=reinterpret_cast<refArray5>(kappa[i]);refArray6 alphak=reinterpret_cast<refArray6>(alpha[i]);refArray7 cutoffk=reinterpret_cast<refArray7>(cutoff[i]);refArray8 one_over_hk=reinterpret_cast<refArray8>(one_over_h[i]);refArray9 cell_volumek=reinterpret_cast<refArray9>(cell_volume[i]);refArray10 Uk=reinterpret_cast<refArray10>(U[0][i]);refArray11 Vk=reinterpret_cast<refArray11>(V[0][i]);refArray12 Sigmak=reinterpret_cast<refArray12>(Sigma[0][i]);refArray13 Q_hatk=reinterpret_cast<refArray13>(Q_hat[0][i]);refArray14 dPdFk=reinterpret_cast<refArray14>(dPdF[0][i]);refArray15 dk=reinterpret_cast<refArray15>(d[0][0][i]);refArray16 system_matrixk=reinterpret_cast<refArray16>(system_matrix[0][i]);Update_Position_Based_State<COROTATED_TAG,__m512,float[16],int[16]>::Run(uk,pk,muk,mu_stabk,kappak,alphak,cutoffk,one_over_hk,cell_volumek,Uk,Vk,Sigmak,Q_hatk,dPdFk,dk,system_matrixk);}
for(int __a=0;__a<9;__a++) for(int __b=0;__b<16;__b++) if(std::abs((U[__a][__b] - U_reference[__a][__b]) / (U_reference[__a][__b])) > 1 ){std::cerr<<"Mismatch detected in MIC implementation"<<std::endl;std::cerr<<"Variable U:"<<std::endl;std::cerr<<"seed="<<seed<<", __a="<<__a<<", __b="<<__b<<std::endl;std::cerr<<"U MIC=  "<<U[__a][__b]<<std::endl;std::cerr<<"U Reference=  "<<U_reference[__a][__b]<<std::endl;std::cerr<<"U Rel Difference=  "<< std::abs((U[__a][__b] - U_reference[__a][__b]) / (U_reference[__a][__b])) << std::endl;std::cerr<<"U Abs Difference=  "<< std::abs(U[__a][__b] - U_reference[__a][__b]) << std::endl;return 1;}
for(int __a=0;__a<9;__a++) for(int __b=0;__b<16;__b++) if(std::abs((V[__a][__b] - V_reference[__a][__b]) / (V_reference[__a][__b])) > 1 ){std::cerr<<"Mismatch detected in MIC implementation"<<std::endl;std::cerr<<"Variable V:"<<std::endl;std::cerr<<"seed="<<seed<<", __a="<<__a<<", __b="<<__b<<std::endl;std::cerr<<"V MIC=  "<<V[__a][__b]<<std::endl;std::cerr<<"V Reference=  "<<V_reference[__a][__b]<<std::endl;std::cerr<<"V Rel Difference=  "<< std::abs((V[__a][__b] - V_reference[__a][__b]) / (V_reference[__a][__b])) << std::endl;std::cerr<<"V Abs Difference=  "<< std::abs(V[__a][__b] - V_reference[__a][__b]) << std::endl;return 1;}
for(int __a=0;__a<3;__a++) for(int __b=0;__b<16;__b++) if(std::abs((Sigma[__a][__b] - Sigma_reference[__a][__b]) / (Sigma_reference[__a][__b])) > 1 ){std::cerr<<"Mismatch detected in MIC implementation"<<std::endl;std::cerr<<"Variable Sigma:"<<std::endl;std::cerr<<"seed="<<seed<<", __a="<<__a<<", __b="<<__b<<std::endl;std::cerr<<"Sigma MIC=  "<<Sigma[__a][__b]<<std::endl;std::cerr<<"Sigma Reference=  "<<Sigma_reference[__a][__b]<<std::endl;std::cerr<<"Sigma Rel Difference=  "<< std::abs((Sigma[__a][__b] - Sigma_reference[__a][__b]) / (Sigma_reference[__a][__b])) << std::endl;std::cerr<<"Sigma Abs Difference=  "<< std::abs(Sigma[__a][__b] - Sigma_reference[__a][__b]) << std::endl;return 1;}
for(int __a=0;__a<3;__a++) for(int __b=0;__b<16;__b++) if(std::abs((Q_hat[__a][__b] - Q_hat_reference[__a][__b]) / (Q_hat_reference[__a][__b])) > 1 ){std::cerr<<"Mismatch detected in MIC implementation"<<std::endl;std::cerr<<"Variable Q_hat:"<<std::endl;std::cerr<<"seed="<<seed<<", __a="<<__a<<", __b="<<__b<<std::endl;std::cerr<<"Q_hat MIC=  "<<Q_hat[__a][__b]<<std::endl;std::cerr<<"Q_hat Reference=  "<<Q_hat_reference[__a][__b]<<std::endl;std::cerr<<"Q_hat Rel Difference=  "<< std::abs((Q_hat[__a][__b] - Q_hat_reference[__a][__b]) / (Q_hat_reference[__a][__b])) << std::endl;std::cerr<<"Q_hat Abs Difference=  "<< std::abs(Q_hat[__a][__b] - Q_hat_reference[__a][__b]) << std::endl;return 1;}
for(int __a=0;__a<12;__a++) for(int __b=0;__b<16;__b++) if(std::abs((dPdF[__a][__b] - dPdF_reference[__a][__b]) / (dPdF_reference[__a][__b])) > 1 ){std::cerr<<"Mismatch detected in MIC implementation"<<std::endl;std::cerr<<"Variable dPdF:"<<std::endl;std::cerr<<"seed="<<seed<<", __a="<<__a<<", __b="<<__b<<std::endl;std::cerr<<"dPdF MIC=  "<<dPdF[__a][__b]<<std::endl;std::cerr<<"dPdF Reference=  "<<dPdF_reference[__a][__b]<<std::endl;std::cerr<<"dPdF Rel Difference=  "<< std::abs((dPdF[__a][__b] - dPdF_reference[__a][__b]) / (dPdF_reference[__a][__b])) << std::endl;std::cerr<<"dPdF Abs Difference=  "<< std::abs(dPdF[__a][__b] - dPdF_reference[__a][__b]) << std::endl;return 1;}
for(int __a=0;__a<3;__a++) for(int __b=0;__b<8;__b++) for(int __c=0;__c<16;__c++) if(std::abs((d[__a][__b][__c] - d_reference[__a][__b][__c]) / (d_reference[__a][__b][__c])) > 1 ){std::cerr<<"Mismatch detected in MIC implementation"<<std::endl;std::cerr<<"Variable d:"<<std::endl;std::cerr<<"seed="<<seed<<", __a="<<__a<<", __b="<<__b<<", __c="<<__c<<std::endl;std::cerr<<"d MIC=  "<<d[__a][__b][__c]<<std::endl;std::cerr<<"d Reference=  "<<d_reference[__a][__b][__c]<<std::endl;std::cerr<<"d Rel Difference=  "<< std::abs((d[__a][__b][__c] - d_reference[__a][__b][__c]) / (d_reference[__a][__b][__c])) << std::endl;std::cerr<<"d Abs Difference=  "<< std::abs(d[__a][__b][__c] - d_reference[__a][__b][__c]) << std::endl;return 1;}
for(int __a=0;__a<300;__a++) for(int __b=0;__b<16;__b++) if(std::abs((system_matrix[__a][__b] - system_matrix_reference[__a][__b]) / (system_matrix_reference[__a][__b])) > 1 ){std::cerr<<"Mismatch detected in MIC implementation"<<std::endl;std::cerr<<"Variable system_matrix:"<<std::endl;std::cerr<<"seed="<<seed<<", __a="<<__a<<", __b="<<__b<<std::endl;std::cerr<<"system_matrix MIC=  "<<system_matrix[__a][__b]<<std::endl;std::cerr<<"system_matrix Reference=  "<<system_matrix_reference[__a][__b]<<std::endl;std::cerr<<"system_matrix Rel Difference=  "<< std::abs((system_matrix[__a][__b] - system_matrix_reference[__a][__b]) / (system_matrix_reference[__a][__b])) << std::endl;std::cerr<<"system_matrix Abs Difference=  "<< std::abs(system_matrix[__a][__b] - system_matrix_reference[__a][__b]) << std::endl;return 1;}

}
#endif

}



    std::cout<<"SIMD check successful!"<<std::endl;

    return 0;

}




#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <sys/time.h>
#include "KernelCommon.h"
struct NEOHOOKEAN_TAG; struct COROTATED_TAG;

#include "Update_Position_Based_State.h"

#define NUM_TRIALS 1000000

template<class T>
T Get_Random(const T a=(T)-1.,const T b=(T)1.)
{
    return ((b-a)*(T)rand())/(T)RAND_MAX+a;
}

struct timeval starttime,stoptime;
void start_timer(){gettimeofday(&starttime,NULL);}
void stop_timer(){gettimeofday(&stoptime,NULL);}
double get_time(){return (double)stoptime.tv_sec-(double)starttime.tv_sec+(double)1e-6*(double)stoptime.tv_usec-(double)1e-6*(double)starttime.tv_usec;}

int main(int argc,char* argv[])
{
    typedef float T;

    std::cout << "Preparing to Run " << NUM_TRIALS << " of all kernels." << std::endl;

    int seed=1;
    if(argc==2) seed=atoi(argv[1]);
    srand(seed);

    

{
std::cout << "Running Stream Test for Update_Position_Based_State " << std::endl;

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
//             COMPUTE SCALAR RESULTS
//
//=======================================================

{
std::cout << "	Running " << NUM_TRIALS << " of SCALAR :  ";
start_timer();
for(int n=0;n<NUM_TRIALS;n++)
{
typedef T (&refArray1)[3][8][16];typedef T (&refArray2)[16];typedef T (&refArray3)[16];typedef T (&refArray4)[16];typedef T (&refArray5)[16];typedef T (&refArray6)[16];typedef T (&refArray7)[16];typedef T (&refArray8)[16];typedef T (&refArray9)[16];typedef T (&refArray10)[9][16];typedef T (&refArray11)[9][16];typedef T (&refArray12)[3][16];typedef T (&refArray13)[3][16];typedef T (&refArray14)[12][16];typedef T (&refArray15)[3][8][16];typedef T (&refArray16)[300][16];for(int i=0; i<16; i+=1) {refArray1 uk=reinterpret_cast<refArray1>(u[0][0][i]);refArray2 pk=reinterpret_cast<refArray2>(p[i]);refArray3 muk=reinterpret_cast<refArray3>(mu[i]);refArray4 mu_stabk=reinterpret_cast<refArray4>(mu_stab[i]);refArray5 kappak=reinterpret_cast<refArray5>(kappa[i]);refArray6 alphak=reinterpret_cast<refArray6>(alpha[i]);refArray7 cutoffk=reinterpret_cast<refArray7>(cutoff[i]);refArray8 one_over_hk=reinterpret_cast<refArray8>(one_over_h[i]);refArray9 cell_volumek=reinterpret_cast<refArray9>(cell_volume[i]);refArray10 Uk=reinterpret_cast<refArray10>(U[0][i]);refArray11 Vk=reinterpret_cast<refArray11>(V[0][i]);refArray12 Sigmak=reinterpret_cast<refArray12>(Sigma[0][i]);refArray13 Q_hatk=reinterpret_cast<refArray13>(Q_hat[0][i]);refArray14 dPdFk=reinterpret_cast<refArray14>(dPdF[0][i]);refArray15 dk=reinterpret_cast<refArray15>(d[0][0][i]);refArray16 system_matrixk=reinterpret_cast<refArray16>(system_matrix[0][i]);Update_Position_Based_State<NEOHOOKEAN_TAG,float,float[16],int[16]>::Run(uk,pk,muk,mu_stabk,kappak,alphak,cutoffk,one_over_hk,cell_volumek,Uk,Vk,Sigmak,Q_hatk,dPdFk,dk,system_matrixk);}

}
stop_timer();
std::cout << get_time() << "s"<< std::endl;
}

//=======================================================
//
//             COMPUTE SSE RESULTS
//
//=======================================================

#ifdef ENABLE_SSE_INSTRUCTION_SET
{
std::cout << "	Running " << NUM_TRIALS << " of SSE :  ";
start_timer();
for(int n=0;n<NUM_TRIALS;n++)
{
typedef T (&refArray1)[3][8][16];typedef T (&refArray2)[16];typedef T (&refArray3)[16];typedef T (&refArray4)[16];typedef T (&refArray5)[16];typedef T (&refArray6)[16];typedef T (&refArray7)[16];typedef T (&refArray8)[16];typedef T (&refArray9)[16];typedef T (&refArray10)[9][16];typedef T (&refArray11)[9][16];typedef T (&refArray12)[3][16];typedef T (&refArray13)[3][16];typedef T (&refArray14)[12][16];typedef T (&refArray15)[3][8][16];typedef T (&refArray16)[300][16];for(int i=0; i<16; i+=4) {refArray1 uk=reinterpret_cast<refArray1>(u[0][0][i]);refArray2 pk=reinterpret_cast<refArray2>(p[i]);refArray3 muk=reinterpret_cast<refArray3>(mu[i]);refArray4 mu_stabk=reinterpret_cast<refArray4>(mu_stab[i]);refArray5 kappak=reinterpret_cast<refArray5>(kappa[i]);refArray6 alphak=reinterpret_cast<refArray6>(alpha[i]);refArray7 cutoffk=reinterpret_cast<refArray7>(cutoff[i]);refArray8 one_over_hk=reinterpret_cast<refArray8>(one_over_h[i]);refArray9 cell_volumek=reinterpret_cast<refArray9>(cell_volume[i]);refArray10 Uk=reinterpret_cast<refArray10>(U[0][i]);refArray11 Vk=reinterpret_cast<refArray11>(V[0][i]);refArray12 Sigmak=reinterpret_cast<refArray12>(Sigma[0][i]);refArray13 Q_hatk=reinterpret_cast<refArray13>(Q_hat[0][i]);refArray14 dPdFk=reinterpret_cast<refArray14>(dPdF[0][i]);refArray15 dk=reinterpret_cast<refArray15>(d[0][0][i]);refArray16 system_matrixk=reinterpret_cast<refArray16>(system_matrix[0][i]);Update_Position_Based_State<NEOHOOKEAN_TAG,__m128,float[16],int[16]>::Run(uk,pk,muk,mu_stabk,kappak,alphak,cutoffk,one_over_hk,cell_volumek,Uk,Vk,Sigmak,Q_hatk,dPdFk,dk,system_matrixk);}

}
stop_timer();
std::cout << get_time() << "s"<< std::endl;
}
#endif

//=======================================================
//
//             COMPUTE AVX RESULTS
//
//=======================================================

#ifdef ENABLE_AVX_INSTRUCTION_SET
{
std::cout << "	Running " << NUM_TRIALS << " of AVX :  ";
start_timer();
for(int n=0;n<NUM_TRIALS;n++)
{
typedef T (&refArray1)[3][8][16];typedef T (&refArray2)[16];typedef T (&refArray3)[16];typedef T (&refArray4)[16];typedef T (&refArray5)[16];typedef T (&refArray6)[16];typedef T (&refArray7)[16];typedef T (&refArray8)[16];typedef T (&refArray9)[16];typedef T (&refArray10)[9][16];typedef T (&refArray11)[9][16];typedef T (&refArray12)[3][16];typedef T (&refArray13)[3][16];typedef T (&refArray14)[12][16];typedef T (&refArray15)[3][8][16];typedef T (&refArray16)[300][16];for(int i=0; i<16; i+=8) {refArray1 uk=reinterpret_cast<refArray1>(u[0][0][i]);refArray2 pk=reinterpret_cast<refArray2>(p[i]);refArray3 muk=reinterpret_cast<refArray3>(mu[i]);refArray4 mu_stabk=reinterpret_cast<refArray4>(mu_stab[i]);refArray5 kappak=reinterpret_cast<refArray5>(kappa[i]);refArray6 alphak=reinterpret_cast<refArray6>(alpha[i]);refArray7 cutoffk=reinterpret_cast<refArray7>(cutoff[i]);refArray8 one_over_hk=reinterpret_cast<refArray8>(one_over_h[i]);refArray9 cell_volumek=reinterpret_cast<refArray9>(cell_volume[i]);refArray10 Uk=reinterpret_cast<refArray10>(U[0][i]);refArray11 Vk=reinterpret_cast<refArray11>(V[0][i]);refArray12 Sigmak=reinterpret_cast<refArray12>(Sigma[0][i]);refArray13 Q_hatk=reinterpret_cast<refArray13>(Q_hat[0][i]);refArray14 dPdFk=reinterpret_cast<refArray14>(dPdF[0][i]);refArray15 dk=reinterpret_cast<refArray15>(d[0][0][i]);refArray16 system_matrixk=reinterpret_cast<refArray16>(system_matrix[0][i]);Update_Position_Based_State<NEOHOOKEAN_TAG,__m256,float[16],int[16]>::Run(uk,pk,muk,mu_stabk,kappak,alphak,cutoffk,one_over_hk,cell_volumek,Uk,Vk,Sigmak,Q_hatk,dPdFk,dk,system_matrixk);}

}
stop_timer();
std::cout << get_time() << "s"<< std::endl;
}
#endif

//=======================================================
//
//             COMPUTE NEON RESULTS
//
//=======================================================

#ifdef ENABLE_NEON_INSTRUCTION_SET
{
std::cout << "	Running " << NUM_TRIALS << " of NEON :  ";
start_timer();
for(int n=0;n<NUM_TRIALS;n++)
{
typedef T (&refArray1)[3][8][16];typedef T (&refArray2)[16];typedef T (&refArray3)[16];typedef T (&refArray4)[16];typedef T (&refArray5)[16];typedef T (&refArray6)[16];typedef T (&refArray7)[16];typedef T (&refArray8)[16];typedef T (&refArray9)[16];typedef T (&refArray10)[9][16];typedef T (&refArray11)[9][16];typedef T (&refArray12)[3][16];typedef T (&refArray13)[3][16];typedef T (&refArray14)[12][16];typedef T (&refArray15)[3][8][16];typedef T (&refArray16)[300][16];for(int i=0; i<16; i+=4) {refArray1 uk=reinterpret_cast<refArray1>(u[0][0][i]);refArray2 pk=reinterpret_cast<refArray2>(p[i]);refArray3 muk=reinterpret_cast<refArray3>(mu[i]);refArray4 mu_stabk=reinterpret_cast<refArray4>(mu_stab[i]);refArray5 kappak=reinterpret_cast<refArray5>(kappa[i]);refArray6 alphak=reinterpret_cast<refArray6>(alpha[i]);refArray7 cutoffk=reinterpret_cast<refArray7>(cutoff[i]);refArray8 one_over_hk=reinterpret_cast<refArray8>(one_over_h[i]);refArray9 cell_volumek=reinterpret_cast<refArray9>(cell_volume[i]);refArray10 Uk=reinterpret_cast<refArray10>(U[0][i]);refArray11 Vk=reinterpret_cast<refArray11>(V[0][i]);refArray12 Sigmak=reinterpret_cast<refArray12>(Sigma[0][i]);refArray13 Q_hatk=reinterpret_cast<refArray13>(Q_hat[0][i]);refArray14 dPdFk=reinterpret_cast<refArray14>(dPdF[0][i]);refArray15 dk=reinterpret_cast<refArray15>(d[0][0][i]);refArray16 system_matrixk=reinterpret_cast<refArray16>(system_matrix[0][i]);Update_Position_Based_State<NEOHOOKEAN_TAG,float32x4_t,float[16],int[16]>::Run(uk,pk,muk,mu_stabk,kappak,alphak,cutoffk,one_over_hk,cell_volumek,Uk,Vk,Sigmak,Q_hatk,dPdFk,dk,system_matrixk);}

}
stop_timer();
std::cout << get_time() << "s"<< std::endl;
}
#endif

//=======================================================
//
//             COMPUTE MIC RESULTS
//
//=======================================================

#ifdef ENABLE_MIC_INSTRUCTION_SET
{
std::cout << "	Running " << NUM_TRIALS << " of MIC :  ";
start_timer();
for(int n=0;n<NUM_TRIALS;n++)
{
typedef T (&refArray1)[3][8][16];typedef T (&refArray2)[16];typedef T (&refArray3)[16];typedef T (&refArray4)[16];typedef T (&refArray5)[16];typedef T (&refArray6)[16];typedef T (&refArray7)[16];typedef T (&refArray8)[16];typedef T (&refArray9)[16];typedef T (&refArray10)[9][16];typedef T (&refArray11)[9][16];typedef T (&refArray12)[3][16];typedef T (&refArray13)[3][16];typedef T (&refArray14)[12][16];typedef T (&refArray15)[3][8][16];typedef T (&refArray16)[300][16];for(int i=0; i<16; i+=16) {refArray1 uk=reinterpret_cast<refArray1>(u[0][0][i]);refArray2 pk=reinterpret_cast<refArray2>(p[i]);refArray3 muk=reinterpret_cast<refArray3>(mu[i]);refArray4 mu_stabk=reinterpret_cast<refArray4>(mu_stab[i]);refArray5 kappak=reinterpret_cast<refArray5>(kappa[i]);refArray6 alphak=reinterpret_cast<refArray6>(alpha[i]);refArray7 cutoffk=reinterpret_cast<refArray7>(cutoff[i]);refArray8 one_over_hk=reinterpret_cast<refArray8>(one_over_h[i]);refArray9 cell_volumek=reinterpret_cast<refArray9>(cell_volume[i]);refArray10 Uk=reinterpret_cast<refArray10>(U[0][i]);refArray11 Vk=reinterpret_cast<refArray11>(V[0][i]);refArray12 Sigmak=reinterpret_cast<refArray12>(Sigma[0][i]);refArray13 Q_hatk=reinterpret_cast<refArray13>(Q_hat[0][i]);refArray14 dPdFk=reinterpret_cast<refArray14>(dPdF[0][i]);refArray15 dk=reinterpret_cast<refArray15>(d[0][0][i]);refArray16 system_matrixk=reinterpret_cast<refArray16>(system_matrix[0][i]);Update_Position_Based_State<NEOHOOKEAN_TAG,__m512,float[16],int[16]>::Run(uk,pk,muk,mu_stabk,kappak,alphak,cutoffk,one_over_hk,cell_volumek,Uk,Vk,Sigmak,Q_hatk,dPdFk,dk,system_matrixk);}

}
stop_timer();
std::cout << get_time() << "s"<< std::endl;
}
#endif

}



{
std::cout << "Running Stream Test for Update_Position_Based_State " << std::endl;

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
//             COMPUTE SCALAR RESULTS
//
//=======================================================

{
std::cout << "	Running " << NUM_TRIALS << " of SCALAR :  ";
start_timer();
for(int n=0;n<NUM_TRIALS;n++)
{
typedef T (&refArray1)[3][8][16];typedef T (&refArray2)[16];typedef T (&refArray3)[16];typedef T (&refArray4)[16];typedef T (&refArray5)[16];typedef T (&refArray6)[16];typedef T (&refArray7)[16];typedef T (&refArray8)[16];typedef T (&refArray9)[16];typedef T (&refArray10)[9][16];typedef T (&refArray11)[9][16];typedef T (&refArray12)[3][16];typedef T (&refArray13)[3][16];typedef T (&refArray14)[12][16];typedef T (&refArray15)[3][8][16];typedef T (&refArray16)[300][16];for(int i=0; i<16; i+=1) {refArray1 uk=reinterpret_cast<refArray1>(u[0][0][i]);refArray2 pk=reinterpret_cast<refArray2>(p[i]);refArray3 muk=reinterpret_cast<refArray3>(mu[i]);refArray4 mu_stabk=reinterpret_cast<refArray4>(mu_stab[i]);refArray5 kappak=reinterpret_cast<refArray5>(kappa[i]);refArray6 alphak=reinterpret_cast<refArray6>(alpha[i]);refArray7 cutoffk=reinterpret_cast<refArray7>(cutoff[i]);refArray8 one_over_hk=reinterpret_cast<refArray8>(one_over_h[i]);refArray9 cell_volumek=reinterpret_cast<refArray9>(cell_volume[i]);refArray10 Uk=reinterpret_cast<refArray10>(U[0][i]);refArray11 Vk=reinterpret_cast<refArray11>(V[0][i]);refArray12 Sigmak=reinterpret_cast<refArray12>(Sigma[0][i]);refArray13 Q_hatk=reinterpret_cast<refArray13>(Q_hat[0][i]);refArray14 dPdFk=reinterpret_cast<refArray14>(dPdF[0][i]);refArray15 dk=reinterpret_cast<refArray15>(d[0][0][i]);refArray16 system_matrixk=reinterpret_cast<refArray16>(system_matrix[0][i]);Update_Position_Based_State<COROTATED_TAG,float,float[16],int[16]>::Run(uk,pk,muk,mu_stabk,kappak,alphak,cutoffk,one_over_hk,cell_volumek,Uk,Vk,Sigmak,Q_hatk,dPdFk,dk,system_matrixk);}

}
stop_timer();
std::cout << get_time() << "s"<< std::endl;
}

//=======================================================
//
//             COMPUTE SSE RESULTS
//
//=======================================================

#ifdef ENABLE_SSE_INSTRUCTION_SET
{
std::cout << "	Running " << NUM_TRIALS << " of SSE :  ";
start_timer();
for(int n=0;n<NUM_TRIALS;n++)
{
typedef T (&refArray1)[3][8][16];typedef T (&refArray2)[16];typedef T (&refArray3)[16];typedef T (&refArray4)[16];typedef T (&refArray5)[16];typedef T (&refArray6)[16];typedef T (&refArray7)[16];typedef T (&refArray8)[16];typedef T (&refArray9)[16];typedef T (&refArray10)[9][16];typedef T (&refArray11)[9][16];typedef T (&refArray12)[3][16];typedef T (&refArray13)[3][16];typedef T (&refArray14)[12][16];typedef T (&refArray15)[3][8][16];typedef T (&refArray16)[300][16];for(int i=0; i<16; i+=4) {refArray1 uk=reinterpret_cast<refArray1>(u[0][0][i]);refArray2 pk=reinterpret_cast<refArray2>(p[i]);refArray3 muk=reinterpret_cast<refArray3>(mu[i]);refArray4 mu_stabk=reinterpret_cast<refArray4>(mu_stab[i]);refArray5 kappak=reinterpret_cast<refArray5>(kappa[i]);refArray6 alphak=reinterpret_cast<refArray6>(alpha[i]);refArray7 cutoffk=reinterpret_cast<refArray7>(cutoff[i]);refArray8 one_over_hk=reinterpret_cast<refArray8>(one_over_h[i]);refArray9 cell_volumek=reinterpret_cast<refArray9>(cell_volume[i]);refArray10 Uk=reinterpret_cast<refArray10>(U[0][i]);refArray11 Vk=reinterpret_cast<refArray11>(V[0][i]);refArray12 Sigmak=reinterpret_cast<refArray12>(Sigma[0][i]);refArray13 Q_hatk=reinterpret_cast<refArray13>(Q_hat[0][i]);refArray14 dPdFk=reinterpret_cast<refArray14>(dPdF[0][i]);refArray15 dk=reinterpret_cast<refArray15>(d[0][0][i]);refArray16 system_matrixk=reinterpret_cast<refArray16>(system_matrix[0][i]);Update_Position_Based_State<COROTATED_TAG,__m128,float[16],int[16]>::Run(uk,pk,muk,mu_stabk,kappak,alphak,cutoffk,one_over_hk,cell_volumek,Uk,Vk,Sigmak,Q_hatk,dPdFk,dk,system_matrixk);}

}
stop_timer();
std::cout << get_time() << "s"<< std::endl;
}
#endif

//=======================================================
//
//             COMPUTE AVX RESULTS
//
//=======================================================

#ifdef ENABLE_AVX_INSTRUCTION_SET
{
std::cout << "	Running " << NUM_TRIALS << " of AVX :  ";
start_timer();
for(int n=0;n<NUM_TRIALS;n++)
{
typedef T (&refArray1)[3][8][16];typedef T (&refArray2)[16];typedef T (&refArray3)[16];typedef T (&refArray4)[16];typedef T (&refArray5)[16];typedef T (&refArray6)[16];typedef T (&refArray7)[16];typedef T (&refArray8)[16];typedef T (&refArray9)[16];typedef T (&refArray10)[9][16];typedef T (&refArray11)[9][16];typedef T (&refArray12)[3][16];typedef T (&refArray13)[3][16];typedef T (&refArray14)[12][16];typedef T (&refArray15)[3][8][16];typedef T (&refArray16)[300][16];for(int i=0; i<16; i+=8) {refArray1 uk=reinterpret_cast<refArray1>(u[0][0][i]);refArray2 pk=reinterpret_cast<refArray2>(p[i]);refArray3 muk=reinterpret_cast<refArray3>(mu[i]);refArray4 mu_stabk=reinterpret_cast<refArray4>(mu_stab[i]);refArray5 kappak=reinterpret_cast<refArray5>(kappa[i]);refArray6 alphak=reinterpret_cast<refArray6>(alpha[i]);refArray7 cutoffk=reinterpret_cast<refArray7>(cutoff[i]);refArray8 one_over_hk=reinterpret_cast<refArray8>(one_over_h[i]);refArray9 cell_volumek=reinterpret_cast<refArray9>(cell_volume[i]);refArray10 Uk=reinterpret_cast<refArray10>(U[0][i]);refArray11 Vk=reinterpret_cast<refArray11>(V[0][i]);refArray12 Sigmak=reinterpret_cast<refArray12>(Sigma[0][i]);refArray13 Q_hatk=reinterpret_cast<refArray13>(Q_hat[0][i]);refArray14 dPdFk=reinterpret_cast<refArray14>(dPdF[0][i]);refArray15 dk=reinterpret_cast<refArray15>(d[0][0][i]);refArray16 system_matrixk=reinterpret_cast<refArray16>(system_matrix[0][i]);Update_Position_Based_State<COROTATED_TAG,__m256,float[16],int[16]>::Run(uk,pk,muk,mu_stabk,kappak,alphak,cutoffk,one_over_hk,cell_volumek,Uk,Vk,Sigmak,Q_hatk,dPdFk,dk,system_matrixk);}

}
stop_timer();
std::cout << get_time() << "s"<< std::endl;
}
#endif

//=======================================================
//
//             COMPUTE NEON RESULTS
//
//=======================================================

#ifdef ENABLE_NEON_INSTRUCTION_SET
{
std::cout << "	Running " << NUM_TRIALS << " of NEON :  ";
start_timer();
for(int n=0;n<NUM_TRIALS;n++)
{
typedef T (&refArray1)[3][8][16];typedef T (&refArray2)[16];typedef T (&refArray3)[16];typedef T (&refArray4)[16];typedef T (&refArray5)[16];typedef T (&refArray6)[16];typedef T (&refArray7)[16];typedef T (&refArray8)[16];typedef T (&refArray9)[16];typedef T (&refArray10)[9][16];typedef T (&refArray11)[9][16];typedef T (&refArray12)[3][16];typedef T (&refArray13)[3][16];typedef T (&refArray14)[12][16];typedef T (&refArray15)[3][8][16];typedef T (&refArray16)[300][16];for(int i=0; i<16; i+=4) {refArray1 uk=reinterpret_cast<refArray1>(u[0][0][i]);refArray2 pk=reinterpret_cast<refArray2>(p[i]);refArray3 muk=reinterpret_cast<refArray3>(mu[i]);refArray4 mu_stabk=reinterpret_cast<refArray4>(mu_stab[i]);refArray5 kappak=reinterpret_cast<refArray5>(kappa[i]);refArray6 alphak=reinterpret_cast<refArray6>(alpha[i]);refArray7 cutoffk=reinterpret_cast<refArray7>(cutoff[i]);refArray8 one_over_hk=reinterpret_cast<refArray8>(one_over_h[i]);refArray9 cell_volumek=reinterpret_cast<refArray9>(cell_volume[i]);refArray10 Uk=reinterpret_cast<refArray10>(U[0][i]);refArray11 Vk=reinterpret_cast<refArray11>(V[0][i]);refArray12 Sigmak=reinterpret_cast<refArray12>(Sigma[0][i]);refArray13 Q_hatk=reinterpret_cast<refArray13>(Q_hat[0][i]);refArray14 dPdFk=reinterpret_cast<refArray14>(dPdF[0][i]);refArray15 dk=reinterpret_cast<refArray15>(d[0][0][i]);refArray16 system_matrixk=reinterpret_cast<refArray16>(system_matrix[0][i]);Update_Position_Based_State<COROTATED_TAG,float32x4_t,float[16],int[16]>::Run(uk,pk,muk,mu_stabk,kappak,alphak,cutoffk,one_over_hk,cell_volumek,Uk,Vk,Sigmak,Q_hatk,dPdFk,dk,system_matrixk);}

}
stop_timer();
std::cout << get_time() << "s"<< std::endl;
}
#endif

//=======================================================
//
//             COMPUTE MIC RESULTS
//
//=======================================================

#ifdef ENABLE_MIC_INSTRUCTION_SET
{
std::cout << "	Running " << NUM_TRIALS << " of MIC :  ";
start_timer();
for(int n=0;n<NUM_TRIALS;n++)
{
typedef T (&refArray1)[3][8][16];typedef T (&refArray2)[16];typedef T (&refArray3)[16];typedef T (&refArray4)[16];typedef T (&refArray5)[16];typedef T (&refArray6)[16];typedef T (&refArray7)[16];typedef T (&refArray8)[16];typedef T (&refArray9)[16];typedef T (&refArray10)[9][16];typedef T (&refArray11)[9][16];typedef T (&refArray12)[3][16];typedef T (&refArray13)[3][16];typedef T (&refArray14)[12][16];typedef T (&refArray15)[3][8][16];typedef T (&refArray16)[300][16];for(int i=0; i<16; i+=16) {refArray1 uk=reinterpret_cast<refArray1>(u[0][0][i]);refArray2 pk=reinterpret_cast<refArray2>(p[i]);refArray3 muk=reinterpret_cast<refArray3>(mu[i]);refArray4 mu_stabk=reinterpret_cast<refArray4>(mu_stab[i]);refArray5 kappak=reinterpret_cast<refArray5>(kappa[i]);refArray6 alphak=reinterpret_cast<refArray6>(alpha[i]);refArray7 cutoffk=reinterpret_cast<refArray7>(cutoff[i]);refArray8 one_over_hk=reinterpret_cast<refArray8>(one_over_h[i]);refArray9 cell_volumek=reinterpret_cast<refArray9>(cell_volume[i]);refArray10 Uk=reinterpret_cast<refArray10>(U[0][i]);refArray11 Vk=reinterpret_cast<refArray11>(V[0][i]);refArray12 Sigmak=reinterpret_cast<refArray12>(Sigma[0][i]);refArray13 Q_hatk=reinterpret_cast<refArray13>(Q_hat[0][i]);refArray14 dPdFk=reinterpret_cast<refArray14>(dPdF[0][i]);refArray15 dk=reinterpret_cast<refArray15>(d[0][0][i]);refArray16 system_matrixk=reinterpret_cast<refArray16>(system_matrix[0][i]);Update_Position_Based_State<COROTATED_TAG,__m512,float[16],int[16]>::Run(uk,pk,muk,mu_stabk,kappak,alphak,cutoffk,one_over_hk,cell_volumek,Uk,Vk,Sigmak,Q_hatk,dPdFk,dk,system_matrixk);}

}
stop_timer();
std::cout << get_time() << "s"<< std::endl;
}
#endif

}



    return 0;

}




#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <sys/time.h>
#include "KernelCommon.h"
struct NEOHOOKEAN_TAG; struct COROTATED_TAG;

#include "Update_Position_Based_State.h"

#include <Thread_Queueing/PTHREAD_QUEUE.h>
#include <Kernel_Serial_Base_Helper.h>


template<class T>
T Get_Random(const T a=(T)-1.,const T b=(T)1.)
{
    return ((b-a)*(T)rand())/(T)RAND_MAX+a;
}

struct timeval starttime,stoptime;
void start_timer(){gettimeofday(&starttime,NULL);}
void stop_timer(){gettimeofday(&stoptime,NULL);}
double get_time(){return (double)stoptime.tv_sec-(double)starttime.tv_sec+(double)1e-6*(double)stoptime.tv_usec-(double)1e-6*(double)starttime.tv_usec;}


template <int SIZE> 
class Update_Position_Based_State_SCALAR_NEOHOOKEAN
    {
    private:
        // Generate Variables Here
        float * _local_u; float * _local_p; float * _local_mu; float * _local_mu_stab; float * _local_kappa; float * _local_alpha; float * _local_cutoff; float * _local_one_over_h; float * _local_cell_volume; float * _local_U; float * _local_V; float * _local_Sigma; float * _local_Q_hat; float * _local_dPdF; float * _local_d; float * _local_system_matrix;

    public:
        explicit Update_Position_Based_State_SCALAR_NEOHOOKEAN(float * u_in,float * p_in,float * mu_in,float * mu_stab_in,float * kappa_in,float * alpha_in,float * cutoff_in,float * one_over_h_in,float * cell_volume_in,float * U_in,float * V_in,float * Sigma_in,float * Q_hat_in,float * dPdF_in,float * d_in,float * system_matrix_in) : _local_u(u_in),_local_p(p_in),_local_mu(mu_in),_local_mu_stab(mu_stab_in),_local_kappa(kappa_in),_local_alpha(alpha_in),_local_cutoff(cutoff_in),_local_one_over_h(one_over_h_in),_local_cell_volume(cell_volume_in),_local_U(U_in),_local_V(V_in),_local_Sigma(Sigma_in),_local_Q_hat(Q_hat_in),_local_dPdF(dPdF_in),_local_d(d_in),_local_system_matrix(system_matrix_in) {}
        void Execute(int index)
        {
        // full array typedefs
        //typedef int (&refArray)[SIZE][3][8];
        typedef float (&fullArray1)[SIZE][3][8][16];typedef float (&fullArray2)[SIZE][16];typedef float (&fullArray3)[SIZE][16];typedef float (&fullArray4)[SIZE][16];typedef float (&fullArray5)[SIZE][16];typedef float (&fullArray6)[SIZE][16];typedef float (&fullArray7)[SIZE][16];typedef float (&fullArray8)[SIZE][16];typedef float (&fullArray9)[SIZE][16];typedef float (&fullArray10)[SIZE][9][16];typedef float (&fullArray11)[SIZE][9][16];typedef float (&fullArray12)[SIZE][3][16];typedef float (&fullArray13)[SIZE][3][16];typedef float (&fullArray14)[SIZE][12][16];typedef float (&fullArray15)[SIZE][3][8][16];typedef float (&fullArray16)[SIZE][300][16];

        // full array extractions
        //refArray _rA = reinterpret_cast < refArray >(*_A);
        fullArray1 _ru =reinterpret_cast<fullArray1>(*_local_u);fullArray2 _rp =reinterpret_cast<fullArray2>(*_local_p);fullArray3 _rmu =reinterpret_cast<fullArray3>(*_local_mu);fullArray4 _rmu_stab =reinterpret_cast<fullArray4>(*_local_mu_stab);fullArray5 _rkappa =reinterpret_cast<fullArray5>(*_local_kappa);fullArray6 _ralpha =reinterpret_cast<fullArray6>(*_local_alpha);fullArray7 _rcutoff =reinterpret_cast<fullArray7>(*_local_cutoff);fullArray8 _rone_over_h =reinterpret_cast<fullArray8>(*_local_one_over_h);fullArray9 _rcell_volume =reinterpret_cast<fullArray9>(*_local_cell_volume);fullArray10 _rU =reinterpret_cast<fullArray10>(*_local_U);fullArray11 _rV =reinterpret_cast<fullArray11>(*_local_V);fullArray12 _rSigma =reinterpret_cast<fullArray12>(*_local_Sigma);fullArray13 _rQ_hat =reinterpret_cast<fullArray13>(*_local_Q_hat);fullArray14 _rdPdF =reinterpret_cast<fullArray14>(*_local_dPdF);fullArray15 _rd =reinterpret_cast<fullArray15>(*_local_d);fullArray16 _rsystem_matrix =reinterpret_cast<fullArray16>(*_local_system_matrix);

        const int ChunkSize = 1;
        // chunk typedef
        //typedef int (&refArray1)[3][8];
        typedef float (&refArray1)[3][8][16];typedef float (&refArray2)[16];typedef float (&refArray3)[16];typedef float (&refArray4)[16];typedef float (&refArray5)[16];typedef float (&refArray6)[16];typedef float (&refArray7)[16];typedef float (&refArray8)[16];typedef float (&refArray9)[16];typedef float (&refArray10)[9][16];typedef float (&refArray11)[9][16];typedef float (&refArray12)[3][16];typedef float (&refArray13)[3][16];typedef float (&refArray14)[12][16];typedef float (&refArray15)[3][8][16];typedef float (&refArray16)[300][16];

        for( int chunk_offset=0; chunk_offset<16; chunk_offset+=ChunkSize)
            {
                //refArray1 Ak = reinterpret_cast < refArray1 > (_rA[entry][0][chunk_offset]);
                refArray1 uk=reinterpret_cast<refArray1>(_ru[index][0][0][chunk_offset]);refArray2 pk=reinterpret_cast<refArray2>(_rp[index][chunk_offset]);refArray3 muk=reinterpret_cast<refArray3>(_rmu[index][chunk_offset]);refArray4 mu_stabk=reinterpret_cast<refArray4>(_rmu_stab[index][chunk_offset]);refArray5 kappak=reinterpret_cast<refArray5>(_rkappa[index][chunk_offset]);refArray6 alphak=reinterpret_cast<refArray6>(_ralpha[index][chunk_offset]);refArray7 cutoffk=reinterpret_cast<refArray7>(_rcutoff[index][chunk_offset]);refArray8 one_over_hk=reinterpret_cast<refArray8>(_rone_over_h[index][chunk_offset]);refArray9 cell_volumek=reinterpret_cast<refArray9>(_rcell_volume[index][chunk_offset]);refArray10 Uk=reinterpret_cast<refArray10>(_rU[index][0][chunk_offset]);refArray11 Vk=reinterpret_cast<refArray11>(_rV[index][0][chunk_offset]);refArray12 Sigmak=reinterpret_cast<refArray12>(_rSigma[index][0][chunk_offset]);refArray13 Q_hatk=reinterpret_cast<refArray13>(_rQ_hat[index][0][chunk_offset]);refArray14 dPdFk=reinterpret_cast<refArray14>(_rdPdF[index][0][chunk_offset]);refArray15 dk=reinterpret_cast<refArray15>(_rd[index][0][0][chunk_offset]);refArray16 system_matrixk=reinterpret_cast<refArray16>(_rsystem_matrix[index][0][chunk_offset]);
               
                Update_Position_Based_State<NEOHOOKEAN_TAG,float,float[16],int[16]>::Run(uk,pk,muk,mu_stabk,kappak,alphak,cutoffk,one_over_hk,cell_volumek,Uk,Vk,Sigmak,Q_hatk,dPdFk,dk,system_matrixk);
             }

        }
    };

template <int SIZE> 
class Update_Position_Based_State_SCALAR_COROTATED
    {
    private:
        // Generate Variables Here
        float * _local_u; float * _local_p; float * _local_mu; float * _local_mu_stab; float * _local_kappa; float * _local_alpha; float * _local_cutoff; float * _local_one_over_h; float * _local_cell_volume; float * _local_U; float * _local_V; float * _local_Sigma; float * _local_Q_hat; float * _local_dPdF; float * _local_d; float * _local_system_matrix;

    public:
        explicit Update_Position_Based_State_SCALAR_COROTATED(float * u_in,float * p_in,float * mu_in,float * mu_stab_in,float * kappa_in,float * alpha_in,float * cutoff_in,float * one_over_h_in,float * cell_volume_in,float * U_in,float * V_in,float * Sigma_in,float * Q_hat_in,float * dPdF_in,float * d_in,float * system_matrix_in) : _local_u(u_in),_local_p(p_in),_local_mu(mu_in),_local_mu_stab(mu_stab_in),_local_kappa(kappa_in),_local_alpha(alpha_in),_local_cutoff(cutoff_in),_local_one_over_h(one_over_h_in),_local_cell_volume(cell_volume_in),_local_U(U_in),_local_V(V_in),_local_Sigma(Sigma_in),_local_Q_hat(Q_hat_in),_local_dPdF(dPdF_in),_local_d(d_in),_local_system_matrix(system_matrix_in) {}
        void Execute(int index)
        {
        // full array typedefs
        //typedef int (&refArray)[SIZE][3][8];
        typedef float (&fullArray1)[SIZE][3][8][16];typedef float (&fullArray2)[SIZE][16];typedef float (&fullArray3)[SIZE][16];typedef float (&fullArray4)[SIZE][16];typedef float (&fullArray5)[SIZE][16];typedef float (&fullArray6)[SIZE][16];typedef float (&fullArray7)[SIZE][16];typedef float (&fullArray8)[SIZE][16];typedef float (&fullArray9)[SIZE][16];typedef float (&fullArray10)[SIZE][9][16];typedef float (&fullArray11)[SIZE][9][16];typedef float (&fullArray12)[SIZE][3][16];typedef float (&fullArray13)[SIZE][3][16];typedef float (&fullArray14)[SIZE][12][16];typedef float (&fullArray15)[SIZE][3][8][16];typedef float (&fullArray16)[SIZE][300][16];

        // full array extractions
        //refArray _rA = reinterpret_cast < refArray >(*_A);
        fullArray1 _ru =reinterpret_cast<fullArray1>(*_local_u);fullArray2 _rp =reinterpret_cast<fullArray2>(*_local_p);fullArray3 _rmu =reinterpret_cast<fullArray3>(*_local_mu);fullArray4 _rmu_stab =reinterpret_cast<fullArray4>(*_local_mu_stab);fullArray5 _rkappa =reinterpret_cast<fullArray5>(*_local_kappa);fullArray6 _ralpha =reinterpret_cast<fullArray6>(*_local_alpha);fullArray7 _rcutoff =reinterpret_cast<fullArray7>(*_local_cutoff);fullArray8 _rone_over_h =reinterpret_cast<fullArray8>(*_local_one_over_h);fullArray9 _rcell_volume =reinterpret_cast<fullArray9>(*_local_cell_volume);fullArray10 _rU =reinterpret_cast<fullArray10>(*_local_U);fullArray11 _rV =reinterpret_cast<fullArray11>(*_local_V);fullArray12 _rSigma =reinterpret_cast<fullArray12>(*_local_Sigma);fullArray13 _rQ_hat =reinterpret_cast<fullArray13>(*_local_Q_hat);fullArray14 _rdPdF =reinterpret_cast<fullArray14>(*_local_dPdF);fullArray15 _rd =reinterpret_cast<fullArray15>(*_local_d);fullArray16 _rsystem_matrix =reinterpret_cast<fullArray16>(*_local_system_matrix);

        const int ChunkSize = 1;
        // chunk typedef
        //typedef int (&refArray1)[3][8];
        typedef float (&refArray1)[3][8][16];typedef float (&refArray2)[16];typedef float (&refArray3)[16];typedef float (&refArray4)[16];typedef float (&refArray5)[16];typedef float (&refArray6)[16];typedef float (&refArray7)[16];typedef float (&refArray8)[16];typedef float (&refArray9)[16];typedef float (&refArray10)[9][16];typedef float (&refArray11)[9][16];typedef float (&refArray12)[3][16];typedef float (&refArray13)[3][16];typedef float (&refArray14)[12][16];typedef float (&refArray15)[3][8][16];typedef float (&refArray16)[300][16];

        for( int chunk_offset=0; chunk_offset<16; chunk_offset+=ChunkSize)
            {
                //refArray1 Ak = reinterpret_cast < refArray1 > (_rA[entry][0][chunk_offset]);
                refArray1 uk=reinterpret_cast<refArray1>(_ru[index][0][0][chunk_offset]);refArray2 pk=reinterpret_cast<refArray2>(_rp[index][chunk_offset]);refArray3 muk=reinterpret_cast<refArray3>(_rmu[index][chunk_offset]);refArray4 mu_stabk=reinterpret_cast<refArray4>(_rmu_stab[index][chunk_offset]);refArray5 kappak=reinterpret_cast<refArray5>(_rkappa[index][chunk_offset]);refArray6 alphak=reinterpret_cast<refArray6>(_ralpha[index][chunk_offset]);refArray7 cutoffk=reinterpret_cast<refArray7>(_rcutoff[index][chunk_offset]);refArray8 one_over_hk=reinterpret_cast<refArray8>(_rone_over_h[index][chunk_offset]);refArray9 cell_volumek=reinterpret_cast<refArray9>(_rcell_volume[index][chunk_offset]);refArray10 Uk=reinterpret_cast<refArray10>(_rU[index][0][chunk_offset]);refArray11 Vk=reinterpret_cast<refArray11>(_rV[index][0][chunk_offset]);refArray12 Sigmak=reinterpret_cast<refArray12>(_rSigma[index][0][chunk_offset]);refArray13 Q_hatk=reinterpret_cast<refArray13>(_rQ_hat[index][0][chunk_offset]);refArray14 dPdFk=reinterpret_cast<refArray14>(_rdPdF[index][0][chunk_offset]);refArray15 dk=reinterpret_cast<refArray15>(_rd[index][0][0][chunk_offset]);refArray16 system_matrixk=reinterpret_cast<refArray16>(_rsystem_matrix[index][0][chunk_offset]);
               
                Update_Position_Based_State<COROTATED_TAG,float,float[16],int[16]>::Run(uk,pk,muk,mu_stabk,kappak,alphak,cutoffk,one_over_hk,cell_volumek,Uk,Vk,Sigmak,Q_hatk,dPdFk,dk,system_matrixk);
             }

        }
    };


#ifdef ENABLE_SSE_INSTRUCTION_SET

template <int SIZE> 
class Update_Position_Based_State_SSE_NEOHOOKEAN
    {
    private:
        // Generate Variables Here
        float * _local_u; float * _local_p; float * _local_mu; float * _local_mu_stab; float * _local_kappa; float * _local_alpha; float * _local_cutoff; float * _local_one_over_h; float * _local_cell_volume; float * _local_U; float * _local_V; float * _local_Sigma; float * _local_Q_hat; float * _local_dPdF; float * _local_d; float * _local_system_matrix;

    public:
        explicit Update_Position_Based_State_SSE_NEOHOOKEAN(float * u_in,float * p_in,float * mu_in,float * mu_stab_in,float * kappa_in,float * alpha_in,float * cutoff_in,float * one_over_h_in,float * cell_volume_in,float * U_in,float * V_in,float * Sigma_in,float * Q_hat_in,float * dPdF_in,float * d_in,float * system_matrix_in) : _local_u(u_in),_local_p(p_in),_local_mu(mu_in),_local_mu_stab(mu_stab_in),_local_kappa(kappa_in),_local_alpha(alpha_in),_local_cutoff(cutoff_in),_local_one_over_h(one_over_h_in),_local_cell_volume(cell_volume_in),_local_U(U_in),_local_V(V_in),_local_Sigma(Sigma_in),_local_Q_hat(Q_hat_in),_local_dPdF(dPdF_in),_local_d(d_in),_local_system_matrix(system_matrix_in) {}
        void Execute(int index)
        {
        // full array typedefs
        //typedef int (&refArray)[SIZE][3][8];
        typedef float (&fullArray1)[SIZE][3][8][16];typedef float (&fullArray2)[SIZE][16];typedef float (&fullArray3)[SIZE][16];typedef float (&fullArray4)[SIZE][16];typedef float (&fullArray5)[SIZE][16];typedef float (&fullArray6)[SIZE][16];typedef float (&fullArray7)[SIZE][16];typedef float (&fullArray8)[SIZE][16];typedef float (&fullArray9)[SIZE][16];typedef float (&fullArray10)[SIZE][9][16];typedef float (&fullArray11)[SIZE][9][16];typedef float (&fullArray12)[SIZE][3][16];typedef float (&fullArray13)[SIZE][3][16];typedef float (&fullArray14)[SIZE][12][16];typedef float (&fullArray15)[SIZE][3][8][16];typedef float (&fullArray16)[SIZE][300][16];

        // full array extractions
        //refArray _rA = reinterpret_cast < refArray >(*_A);
        fullArray1 _ru =reinterpret_cast<fullArray1>(*_local_u);fullArray2 _rp =reinterpret_cast<fullArray2>(*_local_p);fullArray3 _rmu =reinterpret_cast<fullArray3>(*_local_mu);fullArray4 _rmu_stab =reinterpret_cast<fullArray4>(*_local_mu_stab);fullArray5 _rkappa =reinterpret_cast<fullArray5>(*_local_kappa);fullArray6 _ralpha =reinterpret_cast<fullArray6>(*_local_alpha);fullArray7 _rcutoff =reinterpret_cast<fullArray7>(*_local_cutoff);fullArray8 _rone_over_h =reinterpret_cast<fullArray8>(*_local_one_over_h);fullArray9 _rcell_volume =reinterpret_cast<fullArray9>(*_local_cell_volume);fullArray10 _rU =reinterpret_cast<fullArray10>(*_local_U);fullArray11 _rV =reinterpret_cast<fullArray11>(*_local_V);fullArray12 _rSigma =reinterpret_cast<fullArray12>(*_local_Sigma);fullArray13 _rQ_hat =reinterpret_cast<fullArray13>(*_local_Q_hat);fullArray14 _rdPdF =reinterpret_cast<fullArray14>(*_local_dPdF);fullArray15 _rd =reinterpret_cast<fullArray15>(*_local_d);fullArray16 _rsystem_matrix =reinterpret_cast<fullArray16>(*_local_system_matrix);

        const int ChunkSize = 4;
        // chunk typedef
        //typedef int (&refArray1)[3][8];
        typedef float (&refArray1)[3][8][16];typedef float (&refArray2)[16];typedef float (&refArray3)[16];typedef float (&refArray4)[16];typedef float (&refArray5)[16];typedef float (&refArray6)[16];typedef float (&refArray7)[16];typedef float (&refArray8)[16];typedef float (&refArray9)[16];typedef float (&refArray10)[9][16];typedef float (&refArray11)[9][16];typedef float (&refArray12)[3][16];typedef float (&refArray13)[3][16];typedef float (&refArray14)[12][16];typedef float (&refArray15)[3][8][16];typedef float (&refArray16)[300][16];

        for( int chunk_offset=0; chunk_offset<16; chunk_offset+=ChunkSize)
            {
                //refArray1 Ak = reinterpret_cast < refArray1 > (_rA[entry][0][chunk_offset]);
                refArray1 uk=reinterpret_cast<refArray1>(_ru[index][0][0][chunk_offset]);refArray2 pk=reinterpret_cast<refArray2>(_rp[index][chunk_offset]);refArray3 muk=reinterpret_cast<refArray3>(_rmu[index][chunk_offset]);refArray4 mu_stabk=reinterpret_cast<refArray4>(_rmu_stab[index][chunk_offset]);refArray5 kappak=reinterpret_cast<refArray5>(_rkappa[index][chunk_offset]);refArray6 alphak=reinterpret_cast<refArray6>(_ralpha[index][chunk_offset]);refArray7 cutoffk=reinterpret_cast<refArray7>(_rcutoff[index][chunk_offset]);refArray8 one_over_hk=reinterpret_cast<refArray8>(_rone_over_h[index][chunk_offset]);refArray9 cell_volumek=reinterpret_cast<refArray9>(_rcell_volume[index][chunk_offset]);refArray10 Uk=reinterpret_cast<refArray10>(_rU[index][0][chunk_offset]);refArray11 Vk=reinterpret_cast<refArray11>(_rV[index][0][chunk_offset]);refArray12 Sigmak=reinterpret_cast<refArray12>(_rSigma[index][0][chunk_offset]);refArray13 Q_hatk=reinterpret_cast<refArray13>(_rQ_hat[index][0][chunk_offset]);refArray14 dPdFk=reinterpret_cast<refArray14>(_rdPdF[index][0][chunk_offset]);refArray15 dk=reinterpret_cast<refArray15>(_rd[index][0][0][chunk_offset]);refArray16 system_matrixk=reinterpret_cast<refArray16>(_rsystem_matrix[index][0][chunk_offset]);
               
                Update_Position_Based_State<NEOHOOKEAN_TAG,__m128,float[16],int[16]>::Run(uk,pk,muk,mu_stabk,kappak,alphak,cutoffk,one_over_hk,cell_volumek,Uk,Vk,Sigmak,Q_hatk,dPdFk,dk,system_matrixk);
             }

        }
    };

template <int SIZE> 
class Update_Position_Based_State_SSE_COROTATED
    {
    private:
        // Generate Variables Here
        float * _local_u; float * _local_p; float * _local_mu; float * _local_mu_stab; float * _local_kappa; float * _local_alpha; float * _local_cutoff; float * _local_one_over_h; float * _local_cell_volume; float * _local_U; float * _local_V; float * _local_Sigma; float * _local_Q_hat; float * _local_dPdF; float * _local_d; float * _local_system_matrix;

    public:
        explicit Update_Position_Based_State_SSE_COROTATED(float * u_in,float * p_in,float * mu_in,float * mu_stab_in,float * kappa_in,float * alpha_in,float * cutoff_in,float * one_over_h_in,float * cell_volume_in,float * U_in,float * V_in,float * Sigma_in,float * Q_hat_in,float * dPdF_in,float * d_in,float * system_matrix_in) : _local_u(u_in),_local_p(p_in),_local_mu(mu_in),_local_mu_stab(mu_stab_in),_local_kappa(kappa_in),_local_alpha(alpha_in),_local_cutoff(cutoff_in),_local_one_over_h(one_over_h_in),_local_cell_volume(cell_volume_in),_local_U(U_in),_local_V(V_in),_local_Sigma(Sigma_in),_local_Q_hat(Q_hat_in),_local_dPdF(dPdF_in),_local_d(d_in),_local_system_matrix(system_matrix_in) {}
        void Execute(int index)
        {
        // full array typedefs
        //typedef int (&refArray)[SIZE][3][8];
        typedef float (&fullArray1)[SIZE][3][8][16];typedef float (&fullArray2)[SIZE][16];typedef float (&fullArray3)[SIZE][16];typedef float (&fullArray4)[SIZE][16];typedef float (&fullArray5)[SIZE][16];typedef float (&fullArray6)[SIZE][16];typedef float (&fullArray7)[SIZE][16];typedef float (&fullArray8)[SIZE][16];typedef float (&fullArray9)[SIZE][16];typedef float (&fullArray10)[SIZE][9][16];typedef float (&fullArray11)[SIZE][9][16];typedef float (&fullArray12)[SIZE][3][16];typedef float (&fullArray13)[SIZE][3][16];typedef float (&fullArray14)[SIZE][12][16];typedef float (&fullArray15)[SIZE][3][8][16];typedef float (&fullArray16)[SIZE][300][16];

        // full array extractions
        //refArray _rA = reinterpret_cast < refArray >(*_A);
        fullArray1 _ru =reinterpret_cast<fullArray1>(*_local_u);fullArray2 _rp =reinterpret_cast<fullArray2>(*_local_p);fullArray3 _rmu =reinterpret_cast<fullArray3>(*_local_mu);fullArray4 _rmu_stab =reinterpret_cast<fullArray4>(*_local_mu_stab);fullArray5 _rkappa =reinterpret_cast<fullArray5>(*_local_kappa);fullArray6 _ralpha =reinterpret_cast<fullArray6>(*_local_alpha);fullArray7 _rcutoff =reinterpret_cast<fullArray7>(*_local_cutoff);fullArray8 _rone_over_h =reinterpret_cast<fullArray8>(*_local_one_over_h);fullArray9 _rcell_volume =reinterpret_cast<fullArray9>(*_local_cell_volume);fullArray10 _rU =reinterpret_cast<fullArray10>(*_local_U);fullArray11 _rV =reinterpret_cast<fullArray11>(*_local_V);fullArray12 _rSigma =reinterpret_cast<fullArray12>(*_local_Sigma);fullArray13 _rQ_hat =reinterpret_cast<fullArray13>(*_local_Q_hat);fullArray14 _rdPdF =reinterpret_cast<fullArray14>(*_local_dPdF);fullArray15 _rd =reinterpret_cast<fullArray15>(*_local_d);fullArray16 _rsystem_matrix =reinterpret_cast<fullArray16>(*_local_system_matrix);

        const int ChunkSize = 4;
        // chunk typedef
        //typedef int (&refArray1)[3][8];
        typedef float (&refArray1)[3][8][16];typedef float (&refArray2)[16];typedef float (&refArray3)[16];typedef float (&refArray4)[16];typedef float (&refArray5)[16];typedef float (&refArray6)[16];typedef float (&refArray7)[16];typedef float (&refArray8)[16];typedef float (&refArray9)[16];typedef float (&refArray10)[9][16];typedef float (&refArray11)[9][16];typedef float (&refArray12)[3][16];typedef float (&refArray13)[3][16];typedef float (&refArray14)[12][16];typedef float (&refArray15)[3][8][16];typedef float (&refArray16)[300][16];

        for( int chunk_offset=0; chunk_offset<16; chunk_offset+=ChunkSize)
            {
                //refArray1 Ak = reinterpret_cast < refArray1 > (_rA[entry][0][chunk_offset]);
                refArray1 uk=reinterpret_cast<refArray1>(_ru[index][0][0][chunk_offset]);refArray2 pk=reinterpret_cast<refArray2>(_rp[index][chunk_offset]);refArray3 muk=reinterpret_cast<refArray3>(_rmu[index][chunk_offset]);refArray4 mu_stabk=reinterpret_cast<refArray4>(_rmu_stab[index][chunk_offset]);refArray5 kappak=reinterpret_cast<refArray5>(_rkappa[index][chunk_offset]);refArray6 alphak=reinterpret_cast<refArray6>(_ralpha[index][chunk_offset]);refArray7 cutoffk=reinterpret_cast<refArray7>(_rcutoff[index][chunk_offset]);refArray8 one_over_hk=reinterpret_cast<refArray8>(_rone_over_h[index][chunk_offset]);refArray9 cell_volumek=reinterpret_cast<refArray9>(_rcell_volume[index][chunk_offset]);refArray10 Uk=reinterpret_cast<refArray10>(_rU[index][0][chunk_offset]);refArray11 Vk=reinterpret_cast<refArray11>(_rV[index][0][chunk_offset]);refArray12 Sigmak=reinterpret_cast<refArray12>(_rSigma[index][0][chunk_offset]);refArray13 Q_hatk=reinterpret_cast<refArray13>(_rQ_hat[index][0][chunk_offset]);refArray14 dPdFk=reinterpret_cast<refArray14>(_rdPdF[index][0][chunk_offset]);refArray15 dk=reinterpret_cast<refArray15>(_rd[index][0][0][chunk_offset]);refArray16 system_matrixk=reinterpret_cast<refArray16>(_rsystem_matrix[index][0][chunk_offset]);
               
                Update_Position_Based_State<COROTATED_TAG,__m128,float[16],int[16]>::Run(uk,pk,muk,mu_stabk,kappak,alphak,cutoffk,one_over_hk,cell_volumek,Uk,Vk,Sigmak,Q_hatk,dPdFk,dk,system_matrixk);
             }

        }
    };

#endif

#ifdef ENABLE_AVX_INSTRUCTION_SET

template <int SIZE> 
class Update_Position_Based_State_AVX_NEOHOOKEAN
    {
    private:
        // Generate Variables Here
        float * _local_u; float * _local_p; float * _local_mu; float * _local_mu_stab; float * _local_kappa; float * _local_alpha; float * _local_cutoff; float * _local_one_over_h; float * _local_cell_volume; float * _local_U; float * _local_V; float * _local_Sigma; float * _local_Q_hat; float * _local_dPdF; float * _local_d; float * _local_system_matrix;

    public:
        explicit Update_Position_Based_State_AVX_NEOHOOKEAN(float * u_in,float * p_in,float * mu_in,float * mu_stab_in,float * kappa_in,float * alpha_in,float * cutoff_in,float * one_over_h_in,float * cell_volume_in,float * U_in,float * V_in,float * Sigma_in,float * Q_hat_in,float * dPdF_in,float * d_in,float * system_matrix_in) : _local_u(u_in),_local_p(p_in),_local_mu(mu_in),_local_mu_stab(mu_stab_in),_local_kappa(kappa_in),_local_alpha(alpha_in),_local_cutoff(cutoff_in),_local_one_over_h(one_over_h_in),_local_cell_volume(cell_volume_in),_local_U(U_in),_local_V(V_in),_local_Sigma(Sigma_in),_local_Q_hat(Q_hat_in),_local_dPdF(dPdF_in),_local_d(d_in),_local_system_matrix(system_matrix_in) {}
        void Execute(int index)
        {
        // full array typedefs
        //typedef int (&refArray)[SIZE][3][8];
        typedef float (&fullArray1)[SIZE][3][8][16];typedef float (&fullArray2)[SIZE][16];typedef float (&fullArray3)[SIZE][16];typedef float (&fullArray4)[SIZE][16];typedef float (&fullArray5)[SIZE][16];typedef float (&fullArray6)[SIZE][16];typedef float (&fullArray7)[SIZE][16];typedef float (&fullArray8)[SIZE][16];typedef float (&fullArray9)[SIZE][16];typedef float (&fullArray10)[SIZE][9][16];typedef float (&fullArray11)[SIZE][9][16];typedef float (&fullArray12)[SIZE][3][16];typedef float (&fullArray13)[SIZE][3][16];typedef float (&fullArray14)[SIZE][12][16];typedef float (&fullArray15)[SIZE][3][8][16];typedef float (&fullArray16)[SIZE][300][16];

        // full array extractions
        //refArray _rA = reinterpret_cast < refArray >(*_A);
        fullArray1 _ru =reinterpret_cast<fullArray1>(*_local_u);fullArray2 _rp =reinterpret_cast<fullArray2>(*_local_p);fullArray3 _rmu =reinterpret_cast<fullArray3>(*_local_mu);fullArray4 _rmu_stab =reinterpret_cast<fullArray4>(*_local_mu_stab);fullArray5 _rkappa =reinterpret_cast<fullArray5>(*_local_kappa);fullArray6 _ralpha =reinterpret_cast<fullArray6>(*_local_alpha);fullArray7 _rcutoff =reinterpret_cast<fullArray7>(*_local_cutoff);fullArray8 _rone_over_h =reinterpret_cast<fullArray8>(*_local_one_over_h);fullArray9 _rcell_volume =reinterpret_cast<fullArray9>(*_local_cell_volume);fullArray10 _rU =reinterpret_cast<fullArray10>(*_local_U);fullArray11 _rV =reinterpret_cast<fullArray11>(*_local_V);fullArray12 _rSigma =reinterpret_cast<fullArray12>(*_local_Sigma);fullArray13 _rQ_hat =reinterpret_cast<fullArray13>(*_local_Q_hat);fullArray14 _rdPdF =reinterpret_cast<fullArray14>(*_local_dPdF);fullArray15 _rd =reinterpret_cast<fullArray15>(*_local_d);fullArray16 _rsystem_matrix =reinterpret_cast<fullArray16>(*_local_system_matrix);

        const int ChunkSize = 8;
        // chunk typedef
        //typedef int (&refArray1)[3][8];
        typedef float (&refArray1)[3][8][16];typedef float (&refArray2)[16];typedef float (&refArray3)[16];typedef float (&refArray4)[16];typedef float (&refArray5)[16];typedef float (&refArray6)[16];typedef float (&refArray7)[16];typedef float (&refArray8)[16];typedef float (&refArray9)[16];typedef float (&refArray10)[9][16];typedef float (&refArray11)[9][16];typedef float (&refArray12)[3][16];typedef float (&refArray13)[3][16];typedef float (&refArray14)[12][16];typedef float (&refArray15)[3][8][16];typedef float (&refArray16)[300][16];

        for( int chunk_offset=0; chunk_offset<16; chunk_offset+=ChunkSize)
            {
                //refArray1 Ak = reinterpret_cast < refArray1 > (_rA[entry][0][chunk_offset]);
                refArray1 uk=reinterpret_cast<refArray1>(_ru[index][0][0][chunk_offset]);refArray2 pk=reinterpret_cast<refArray2>(_rp[index][chunk_offset]);refArray3 muk=reinterpret_cast<refArray3>(_rmu[index][chunk_offset]);refArray4 mu_stabk=reinterpret_cast<refArray4>(_rmu_stab[index][chunk_offset]);refArray5 kappak=reinterpret_cast<refArray5>(_rkappa[index][chunk_offset]);refArray6 alphak=reinterpret_cast<refArray6>(_ralpha[index][chunk_offset]);refArray7 cutoffk=reinterpret_cast<refArray7>(_rcutoff[index][chunk_offset]);refArray8 one_over_hk=reinterpret_cast<refArray8>(_rone_over_h[index][chunk_offset]);refArray9 cell_volumek=reinterpret_cast<refArray9>(_rcell_volume[index][chunk_offset]);refArray10 Uk=reinterpret_cast<refArray10>(_rU[index][0][chunk_offset]);refArray11 Vk=reinterpret_cast<refArray11>(_rV[index][0][chunk_offset]);refArray12 Sigmak=reinterpret_cast<refArray12>(_rSigma[index][0][chunk_offset]);refArray13 Q_hatk=reinterpret_cast<refArray13>(_rQ_hat[index][0][chunk_offset]);refArray14 dPdFk=reinterpret_cast<refArray14>(_rdPdF[index][0][chunk_offset]);refArray15 dk=reinterpret_cast<refArray15>(_rd[index][0][0][chunk_offset]);refArray16 system_matrixk=reinterpret_cast<refArray16>(_rsystem_matrix[index][0][chunk_offset]);
               
                Update_Position_Based_State<NEOHOOKEAN_TAG,__m256,float[16],int[16]>::Run(uk,pk,muk,mu_stabk,kappak,alphak,cutoffk,one_over_hk,cell_volumek,Uk,Vk,Sigmak,Q_hatk,dPdFk,dk,system_matrixk);
             }

        }
    };

template <int SIZE> 
class Update_Position_Based_State_AVX_COROTATED
    {
    private:
        // Generate Variables Here
        float * _local_u; float * _local_p; float * _local_mu; float * _local_mu_stab; float * _local_kappa; float * _local_alpha; float * _local_cutoff; float * _local_one_over_h; float * _local_cell_volume; float * _local_U; float * _local_V; float * _local_Sigma; float * _local_Q_hat; float * _local_dPdF; float * _local_d; float * _local_system_matrix;

    public:
        explicit Update_Position_Based_State_AVX_COROTATED(float * u_in,float * p_in,float * mu_in,float * mu_stab_in,float * kappa_in,float * alpha_in,float * cutoff_in,float * one_over_h_in,float * cell_volume_in,float * U_in,float * V_in,float * Sigma_in,float * Q_hat_in,float * dPdF_in,float * d_in,float * system_matrix_in) : _local_u(u_in),_local_p(p_in),_local_mu(mu_in),_local_mu_stab(mu_stab_in),_local_kappa(kappa_in),_local_alpha(alpha_in),_local_cutoff(cutoff_in),_local_one_over_h(one_over_h_in),_local_cell_volume(cell_volume_in),_local_U(U_in),_local_V(V_in),_local_Sigma(Sigma_in),_local_Q_hat(Q_hat_in),_local_dPdF(dPdF_in),_local_d(d_in),_local_system_matrix(system_matrix_in) {}
        void Execute(int index)
        {
        // full array typedefs
        //typedef int (&refArray)[SIZE][3][8];
        typedef float (&fullArray1)[SIZE][3][8][16];typedef float (&fullArray2)[SIZE][16];typedef float (&fullArray3)[SIZE][16];typedef float (&fullArray4)[SIZE][16];typedef float (&fullArray5)[SIZE][16];typedef float (&fullArray6)[SIZE][16];typedef float (&fullArray7)[SIZE][16];typedef float (&fullArray8)[SIZE][16];typedef float (&fullArray9)[SIZE][16];typedef float (&fullArray10)[SIZE][9][16];typedef float (&fullArray11)[SIZE][9][16];typedef float (&fullArray12)[SIZE][3][16];typedef float (&fullArray13)[SIZE][3][16];typedef float (&fullArray14)[SIZE][12][16];typedef float (&fullArray15)[SIZE][3][8][16];typedef float (&fullArray16)[SIZE][300][16];

        // full array extractions
        //refArray _rA = reinterpret_cast < refArray >(*_A);
        fullArray1 _ru =reinterpret_cast<fullArray1>(*_local_u);fullArray2 _rp =reinterpret_cast<fullArray2>(*_local_p);fullArray3 _rmu =reinterpret_cast<fullArray3>(*_local_mu);fullArray4 _rmu_stab =reinterpret_cast<fullArray4>(*_local_mu_stab);fullArray5 _rkappa =reinterpret_cast<fullArray5>(*_local_kappa);fullArray6 _ralpha =reinterpret_cast<fullArray6>(*_local_alpha);fullArray7 _rcutoff =reinterpret_cast<fullArray7>(*_local_cutoff);fullArray8 _rone_over_h =reinterpret_cast<fullArray8>(*_local_one_over_h);fullArray9 _rcell_volume =reinterpret_cast<fullArray9>(*_local_cell_volume);fullArray10 _rU =reinterpret_cast<fullArray10>(*_local_U);fullArray11 _rV =reinterpret_cast<fullArray11>(*_local_V);fullArray12 _rSigma =reinterpret_cast<fullArray12>(*_local_Sigma);fullArray13 _rQ_hat =reinterpret_cast<fullArray13>(*_local_Q_hat);fullArray14 _rdPdF =reinterpret_cast<fullArray14>(*_local_dPdF);fullArray15 _rd =reinterpret_cast<fullArray15>(*_local_d);fullArray16 _rsystem_matrix =reinterpret_cast<fullArray16>(*_local_system_matrix);

        const int ChunkSize = 8;
        // chunk typedef
        //typedef int (&refArray1)[3][8];
        typedef float (&refArray1)[3][8][16];typedef float (&refArray2)[16];typedef float (&refArray3)[16];typedef float (&refArray4)[16];typedef float (&refArray5)[16];typedef float (&refArray6)[16];typedef float (&refArray7)[16];typedef float (&refArray8)[16];typedef float (&refArray9)[16];typedef float (&refArray10)[9][16];typedef float (&refArray11)[9][16];typedef float (&refArray12)[3][16];typedef float (&refArray13)[3][16];typedef float (&refArray14)[12][16];typedef float (&refArray15)[3][8][16];typedef float (&refArray16)[300][16];

        for( int chunk_offset=0; chunk_offset<16; chunk_offset+=ChunkSize)
            {
                //refArray1 Ak = reinterpret_cast < refArray1 > (_rA[entry][0][chunk_offset]);
                refArray1 uk=reinterpret_cast<refArray1>(_ru[index][0][0][chunk_offset]);refArray2 pk=reinterpret_cast<refArray2>(_rp[index][chunk_offset]);refArray3 muk=reinterpret_cast<refArray3>(_rmu[index][chunk_offset]);refArray4 mu_stabk=reinterpret_cast<refArray4>(_rmu_stab[index][chunk_offset]);refArray5 kappak=reinterpret_cast<refArray5>(_rkappa[index][chunk_offset]);refArray6 alphak=reinterpret_cast<refArray6>(_ralpha[index][chunk_offset]);refArray7 cutoffk=reinterpret_cast<refArray7>(_rcutoff[index][chunk_offset]);refArray8 one_over_hk=reinterpret_cast<refArray8>(_rone_over_h[index][chunk_offset]);refArray9 cell_volumek=reinterpret_cast<refArray9>(_rcell_volume[index][chunk_offset]);refArray10 Uk=reinterpret_cast<refArray10>(_rU[index][0][chunk_offset]);refArray11 Vk=reinterpret_cast<refArray11>(_rV[index][0][chunk_offset]);refArray12 Sigmak=reinterpret_cast<refArray12>(_rSigma[index][0][chunk_offset]);refArray13 Q_hatk=reinterpret_cast<refArray13>(_rQ_hat[index][0][chunk_offset]);refArray14 dPdFk=reinterpret_cast<refArray14>(_rdPdF[index][0][chunk_offset]);refArray15 dk=reinterpret_cast<refArray15>(_rd[index][0][0][chunk_offset]);refArray16 system_matrixk=reinterpret_cast<refArray16>(_rsystem_matrix[index][0][chunk_offset]);
               
                Update_Position_Based_State<COROTATED_TAG,__m256,float[16],int[16]>::Run(uk,pk,muk,mu_stabk,kappak,alphak,cutoffk,one_over_hk,cell_volumek,Uk,Vk,Sigmak,Q_hatk,dPdFk,dk,system_matrixk);
             }

        }
    };

#endif

#ifdef ENABLE_NEON_INSTRUCTION_SET

template <int SIZE> 
class Update_Position_Based_State_NEON_NEOHOOKEAN
    {
    private:
        // Generate Variables Here
        float * _local_u; float * _local_p; float * _local_mu; float * _local_mu_stab; float * _local_kappa; float * _local_alpha; float * _local_cutoff; float * _local_one_over_h; float * _local_cell_volume; float * _local_U; float * _local_V; float * _local_Sigma; float * _local_Q_hat; float * _local_dPdF; float * _local_d; float * _local_system_matrix;

    public:
        explicit Update_Position_Based_State_NEON_NEOHOOKEAN(float * u_in,float * p_in,float * mu_in,float * mu_stab_in,float * kappa_in,float * alpha_in,float * cutoff_in,float * one_over_h_in,float * cell_volume_in,float * U_in,float * V_in,float * Sigma_in,float * Q_hat_in,float * dPdF_in,float * d_in,float * system_matrix_in) : _local_u(u_in),_local_p(p_in),_local_mu(mu_in),_local_mu_stab(mu_stab_in),_local_kappa(kappa_in),_local_alpha(alpha_in),_local_cutoff(cutoff_in),_local_one_over_h(one_over_h_in),_local_cell_volume(cell_volume_in),_local_U(U_in),_local_V(V_in),_local_Sigma(Sigma_in),_local_Q_hat(Q_hat_in),_local_dPdF(dPdF_in),_local_d(d_in),_local_system_matrix(system_matrix_in) {}
        void Execute(int index)
        {
        // full array typedefs
        //typedef int (&refArray)[SIZE][3][8];
        typedef float (&fullArray1)[SIZE][3][8][16];typedef float (&fullArray2)[SIZE][16];typedef float (&fullArray3)[SIZE][16];typedef float (&fullArray4)[SIZE][16];typedef float (&fullArray5)[SIZE][16];typedef float (&fullArray6)[SIZE][16];typedef float (&fullArray7)[SIZE][16];typedef float (&fullArray8)[SIZE][16];typedef float (&fullArray9)[SIZE][16];typedef float (&fullArray10)[SIZE][9][16];typedef float (&fullArray11)[SIZE][9][16];typedef float (&fullArray12)[SIZE][3][16];typedef float (&fullArray13)[SIZE][3][16];typedef float (&fullArray14)[SIZE][12][16];typedef float (&fullArray15)[SIZE][3][8][16];typedef float (&fullArray16)[SIZE][300][16];

        // full array extractions
        //refArray _rA = reinterpret_cast < refArray >(*_A);
        fullArray1 _ru =reinterpret_cast<fullArray1>(*_local_u);fullArray2 _rp =reinterpret_cast<fullArray2>(*_local_p);fullArray3 _rmu =reinterpret_cast<fullArray3>(*_local_mu);fullArray4 _rmu_stab =reinterpret_cast<fullArray4>(*_local_mu_stab);fullArray5 _rkappa =reinterpret_cast<fullArray5>(*_local_kappa);fullArray6 _ralpha =reinterpret_cast<fullArray6>(*_local_alpha);fullArray7 _rcutoff =reinterpret_cast<fullArray7>(*_local_cutoff);fullArray8 _rone_over_h =reinterpret_cast<fullArray8>(*_local_one_over_h);fullArray9 _rcell_volume =reinterpret_cast<fullArray9>(*_local_cell_volume);fullArray10 _rU =reinterpret_cast<fullArray10>(*_local_U);fullArray11 _rV =reinterpret_cast<fullArray11>(*_local_V);fullArray12 _rSigma =reinterpret_cast<fullArray12>(*_local_Sigma);fullArray13 _rQ_hat =reinterpret_cast<fullArray13>(*_local_Q_hat);fullArray14 _rdPdF =reinterpret_cast<fullArray14>(*_local_dPdF);fullArray15 _rd =reinterpret_cast<fullArray15>(*_local_d);fullArray16 _rsystem_matrix =reinterpret_cast<fullArray16>(*_local_system_matrix);

        const int ChunkSize = 4;
        // chunk typedef
        //typedef int (&refArray1)[3][8];
        typedef float (&refArray1)[3][8][16];typedef float (&refArray2)[16];typedef float (&refArray3)[16];typedef float (&refArray4)[16];typedef float (&refArray5)[16];typedef float (&refArray6)[16];typedef float (&refArray7)[16];typedef float (&refArray8)[16];typedef float (&refArray9)[16];typedef float (&refArray10)[9][16];typedef float (&refArray11)[9][16];typedef float (&refArray12)[3][16];typedef float (&refArray13)[3][16];typedef float (&refArray14)[12][16];typedef float (&refArray15)[3][8][16];typedef float (&refArray16)[300][16];

        for( int chunk_offset=0; chunk_offset<16; chunk_offset+=ChunkSize)
            {
                //refArray1 Ak = reinterpret_cast < refArray1 > (_rA[entry][0][chunk_offset]);
                refArray1 uk=reinterpret_cast<refArray1>(_ru[index][0][0][chunk_offset]);refArray2 pk=reinterpret_cast<refArray2>(_rp[index][chunk_offset]);refArray3 muk=reinterpret_cast<refArray3>(_rmu[index][chunk_offset]);refArray4 mu_stabk=reinterpret_cast<refArray4>(_rmu_stab[index][chunk_offset]);refArray5 kappak=reinterpret_cast<refArray5>(_rkappa[index][chunk_offset]);refArray6 alphak=reinterpret_cast<refArray6>(_ralpha[index][chunk_offset]);refArray7 cutoffk=reinterpret_cast<refArray7>(_rcutoff[index][chunk_offset]);refArray8 one_over_hk=reinterpret_cast<refArray8>(_rone_over_h[index][chunk_offset]);refArray9 cell_volumek=reinterpret_cast<refArray9>(_rcell_volume[index][chunk_offset]);refArray10 Uk=reinterpret_cast<refArray10>(_rU[index][0][chunk_offset]);refArray11 Vk=reinterpret_cast<refArray11>(_rV[index][0][chunk_offset]);refArray12 Sigmak=reinterpret_cast<refArray12>(_rSigma[index][0][chunk_offset]);refArray13 Q_hatk=reinterpret_cast<refArray13>(_rQ_hat[index][0][chunk_offset]);refArray14 dPdFk=reinterpret_cast<refArray14>(_rdPdF[index][0][chunk_offset]);refArray15 dk=reinterpret_cast<refArray15>(_rd[index][0][0][chunk_offset]);refArray16 system_matrixk=reinterpret_cast<refArray16>(_rsystem_matrix[index][0][chunk_offset]);
               
                Update_Position_Based_State<NEOHOOKEAN_TAG,float32x4_t,float[16],int[16]>::Run(uk,pk,muk,mu_stabk,kappak,alphak,cutoffk,one_over_hk,cell_volumek,Uk,Vk,Sigmak,Q_hatk,dPdFk,dk,system_matrixk);
             }

        }
    };

template <int SIZE> 
class Update_Position_Based_State_NEON_COROTATED
    {
    private:
        // Generate Variables Here
        float * _local_u; float * _local_p; float * _local_mu; float * _local_mu_stab; float * _local_kappa; float * _local_alpha; float * _local_cutoff; float * _local_one_over_h; float * _local_cell_volume; float * _local_U; float * _local_V; float * _local_Sigma; float * _local_Q_hat; float * _local_dPdF; float * _local_d; float * _local_system_matrix;

    public:
        explicit Update_Position_Based_State_NEON_COROTATED(float * u_in,float * p_in,float * mu_in,float * mu_stab_in,float * kappa_in,float * alpha_in,float * cutoff_in,float * one_over_h_in,float * cell_volume_in,float * U_in,float * V_in,float * Sigma_in,float * Q_hat_in,float * dPdF_in,float * d_in,float * system_matrix_in) : _local_u(u_in),_local_p(p_in),_local_mu(mu_in),_local_mu_stab(mu_stab_in),_local_kappa(kappa_in),_local_alpha(alpha_in),_local_cutoff(cutoff_in),_local_one_over_h(one_over_h_in),_local_cell_volume(cell_volume_in),_local_U(U_in),_local_V(V_in),_local_Sigma(Sigma_in),_local_Q_hat(Q_hat_in),_local_dPdF(dPdF_in),_local_d(d_in),_local_system_matrix(system_matrix_in) {}
        void Execute(int index)
        {
        // full array typedefs
        //typedef int (&refArray)[SIZE][3][8];
        typedef float (&fullArray1)[SIZE][3][8][16];typedef float (&fullArray2)[SIZE][16];typedef float (&fullArray3)[SIZE][16];typedef float (&fullArray4)[SIZE][16];typedef float (&fullArray5)[SIZE][16];typedef float (&fullArray6)[SIZE][16];typedef float (&fullArray7)[SIZE][16];typedef float (&fullArray8)[SIZE][16];typedef float (&fullArray9)[SIZE][16];typedef float (&fullArray10)[SIZE][9][16];typedef float (&fullArray11)[SIZE][9][16];typedef float (&fullArray12)[SIZE][3][16];typedef float (&fullArray13)[SIZE][3][16];typedef float (&fullArray14)[SIZE][12][16];typedef float (&fullArray15)[SIZE][3][8][16];typedef float (&fullArray16)[SIZE][300][16];

        // full array extractions
        //refArray _rA = reinterpret_cast < refArray >(*_A);
        fullArray1 _ru =reinterpret_cast<fullArray1>(*_local_u);fullArray2 _rp =reinterpret_cast<fullArray2>(*_local_p);fullArray3 _rmu =reinterpret_cast<fullArray3>(*_local_mu);fullArray4 _rmu_stab =reinterpret_cast<fullArray4>(*_local_mu_stab);fullArray5 _rkappa =reinterpret_cast<fullArray5>(*_local_kappa);fullArray6 _ralpha =reinterpret_cast<fullArray6>(*_local_alpha);fullArray7 _rcutoff =reinterpret_cast<fullArray7>(*_local_cutoff);fullArray8 _rone_over_h =reinterpret_cast<fullArray8>(*_local_one_over_h);fullArray9 _rcell_volume =reinterpret_cast<fullArray9>(*_local_cell_volume);fullArray10 _rU =reinterpret_cast<fullArray10>(*_local_U);fullArray11 _rV =reinterpret_cast<fullArray11>(*_local_V);fullArray12 _rSigma =reinterpret_cast<fullArray12>(*_local_Sigma);fullArray13 _rQ_hat =reinterpret_cast<fullArray13>(*_local_Q_hat);fullArray14 _rdPdF =reinterpret_cast<fullArray14>(*_local_dPdF);fullArray15 _rd =reinterpret_cast<fullArray15>(*_local_d);fullArray16 _rsystem_matrix =reinterpret_cast<fullArray16>(*_local_system_matrix);

        const int ChunkSize = 4;
        // chunk typedef
        //typedef int (&refArray1)[3][8];
        typedef float (&refArray1)[3][8][16];typedef float (&refArray2)[16];typedef float (&refArray3)[16];typedef float (&refArray4)[16];typedef float (&refArray5)[16];typedef float (&refArray6)[16];typedef float (&refArray7)[16];typedef float (&refArray8)[16];typedef float (&refArray9)[16];typedef float (&refArray10)[9][16];typedef float (&refArray11)[9][16];typedef float (&refArray12)[3][16];typedef float (&refArray13)[3][16];typedef float (&refArray14)[12][16];typedef float (&refArray15)[3][8][16];typedef float (&refArray16)[300][16];

        for( int chunk_offset=0; chunk_offset<16; chunk_offset+=ChunkSize)
            {
                //refArray1 Ak = reinterpret_cast < refArray1 > (_rA[entry][0][chunk_offset]);
                refArray1 uk=reinterpret_cast<refArray1>(_ru[index][0][0][chunk_offset]);refArray2 pk=reinterpret_cast<refArray2>(_rp[index][chunk_offset]);refArray3 muk=reinterpret_cast<refArray3>(_rmu[index][chunk_offset]);refArray4 mu_stabk=reinterpret_cast<refArray4>(_rmu_stab[index][chunk_offset]);refArray5 kappak=reinterpret_cast<refArray5>(_rkappa[index][chunk_offset]);refArray6 alphak=reinterpret_cast<refArray6>(_ralpha[index][chunk_offset]);refArray7 cutoffk=reinterpret_cast<refArray7>(_rcutoff[index][chunk_offset]);refArray8 one_over_hk=reinterpret_cast<refArray8>(_rone_over_h[index][chunk_offset]);refArray9 cell_volumek=reinterpret_cast<refArray9>(_rcell_volume[index][chunk_offset]);refArray10 Uk=reinterpret_cast<refArray10>(_rU[index][0][chunk_offset]);refArray11 Vk=reinterpret_cast<refArray11>(_rV[index][0][chunk_offset]);refArray12 Sigmak=reinterpret_cast<refArray12>(_rSigma[index][0][chunk_offset]);refArray13 Q_hatk=reinterpret_cast<refArray13>(_rQ_hat[index][0][chunk_offset]);refArray14 dPdFk=reinterpret_cast<refArray14>(_rdPdF[index][0][chunk_offset]);refArray15 dk=reinterpret_cast<refArray15>(_rd[index][0][0][chunk_offset]);refArray16 system_matrixk=reinterpret_cast<refArray16>(_rsystem_matrix[index][0][chunk_offset]);
               
                Update_Position_Based_State<COROTATED_TAG,float32x4_t,float[16],int[16]>::Run(uk,pk,muk,mu_stabk,kappak,alphak,cutoffk,one_over_hk,cell_volumek,Uk,Vk,Sigmak,Q_hatk,dPdFk,dk,system_matrixk);
             }

        }
    };

#endif

#ifdef ENABLE_MIC_INSTRUCTION_SET

template <int SIZE> 
class Update_Position_Based_State_MIC_NEOHOOKEAN
    {
    private:
        // Generate Variables Here
        float * _local_u; float * _local_p; float * _local_mu; float * _local_mu_stab; float * _local_kappa; float * _local_alpha; float * _local_cutoff; float * _local_one_over_h; float * _local_cell_volume; float * _local_U; float * _local_V; float * _local_Sigma; float * _local_Q_hat; float * _local_dPdF; float * _local_d; float * _local_system_matrix;

    public:
        explicit Update_Position_Based_State_MIC_NEOHOOKEAN(float * u_in,float * p_in,float * mu_in,float * mu_stab_in,float * kappa_in,float * alpha_in,float * cutoff_in,float * one_over_h_in,float * cell_volume_in,float * U_in,float * V_in,float * Sigma_in,float * Q_hat_in,float * dPdF_in,float * d_in,float * system_matrix_in) : _local_u(u_in),_local_p(p_in),_local_mu(mu_in),_local_mu_stab(mu_stab_in),_local_kappa(kappa_in),_local_alpha(alpha_in),_local_cutoff(cutoff_in),_local_one_over_h(one_over_h_in),_local_cell_volume(cell_volume_in),_local_U(U_in),_local_V(V_in),_local_Sigma(Sigma_in),_local_Q_hat(Q_hat_in),_local_dPdF(dPdF_in),_local_d(d_in),_local_system_matrix(system_matrix_in) {}
        void Execute(int index)
        {
        // full array typedefs
        //typedef int (&refArray)[SIZE][3][8];
        typedef float (&fullArray1)[SIZE][3][8][16];typedef float (&fullArray2)[SIZE][16];typedef float (&fullArray3)[SIZE][16];typedef float (&fullArray4)[SIZE][16];typedef float (&fullArray5)[SIZE][16];typedef float (&fullArray6)[SIZE][16];typedef float (&fullArray7)[SIZE][16];typedef float (&fullArray8)[SIZE][16];typedef float (&fullArray9)[SIZE][16];typedef float (&fullArray10)[SIZE][9][16];typedef float (&fullArray11)[SIZE][9][16];typedef float (&fullArray12)[SIZE][3][16];typedef float (&fullArray13)[SIZE][3][16];typedef float (&fullArray14)[SIZE][12][16];typedef float (&fullArray15)[SIZE][3][8][16];typedef float (&fullArray16)[SIZE][300][16];

        // full array extractions
        //refArray _rA = reinterpret_cast < refArray >(*_A);
        fullArray1 _ru =reinterpret_cast<fullArray1>(*_local_u);fullArray2 _rp =reinterpret_cast<fullArray2>(*_local_p);fullArray3 _rmu =reinterpret_cast<fullArray3>(*_local_mu);fullArray4 _rmu_stab =reinterpret_cast<fullArray4>(*_local_mu_stab);fullArray5 _rkappa =reinterpret_cast<fullArray5>(*_local_kappa);fullArray6 _ralpha =reinterpret_cast<fullArray6>(*_local_alpha);fullArray7 _rcutoff =reinterpret_cast<fullArray7>(*_local_cutoff);fullArray8 _rone_over_h =reinterpret_cast<fullArray8>(*_local_one_over_h);fullArray9 _rcell_volume =reinterpret_cast<fullArray9>(*_local_cell_volume);fullArray10 _rU =reinterpret_cast<fullArray10>(*_local_U);fullArray11 _rV =reinterpret_cast<fullArray11>(*_local_V);fullArray12 _rSigma =reinterpret_cast<fullArray12>(*_local_Sigma);fullArray13 _rQ_hat =reinterpret_cast<fullArray13>(*_local_Q_hat);fullArray14 _rdPdF =reinterpret_cast<fullArray14>(*_local_dPdF);fullArray15 _rd =reinterpret_cast<fullArray15>(*_local_d);fullArray16 _rsystem_matrix =reinterpret_cast<fullArray16>(*_local_system_matrix);

        const int ChunkSize = 16;
        // chunk typedef
        //typedef int (&refArray1)[3][8];
        typedef float (&refArray1)[3][8][16];typedef float (&refArray2)[16];typedef float (&refArray3)[16];typedef float (&refArray4)[16];typedef float (&refArray5)[16];typedef float (&refArray6)[16];typedef float (&refArray7)[16];typedef float (&refArray8)[16];typedef float (&refArray9)[16];typedef float (&refArray10)[9][16];typedef float (&refArray11)[9][16];typedef float (&refArray12)[3][16];typedef float (&refArray13)[3][16];typedef float (&refArray14)[12][16];typedef float (&refArray15)[3][8][16];typedef float (&refArray16)[300][16];

        for( int chunk_offset=0; chunk_offset<16; chunk_offset+=ChunkSize)
            {
                //refArray1 Ak = reinterpret_cast < refArray1 > (_rA[entry][0][chunk_offset]);
                refArray1 uk=reinterpret_cast<refArray1>(_ru[index][0][0][chunk_offset]);refArray2 pk=reinterpret_cast<refArray2>(_rp[index][chunk_offset]);refArray3 muk=reinterpret_cast<refArray3>(_rmu[index][chunk_offset]);refArray4 mu_stabk=reinterpret_cast<refArray4>(_rmu_stab[index][chunk_offset]);refArray5 kappak=reinterpret_cast<refArray5>(_rkappa[index][chunk_offset]);refArray6 alphak=reinterpret_cast<refArray6>(_ralpha[index][chunk_offset]);refArray7 cutoffk=reinterpret_cast<refArray7>(_rcutoff[index][chunk_offset]);refArray8 one_over_hk=reinterpret_cast<refArray8>(_rone_over_h[index][chunk_offset]);refArray9 cell_volumek=reinterpret_cast<refArray9>(_rcell_volume[index][chunk_offset]);refArray10 Uk=reinterpret_cast<refArray10>(_rU[index][0][chunk_offset]);refArray11 Vk=reinterpret_cast<refArray11>(_rV[index][0][chunk_offset]);refArray12 Sigmak=reinterpret_cast<refArray12>(_rSigma[index][0][chunk_offset]);refArray13 Q_hatk=reinterpret_cast<refArray13>(_rQ_hat[index][0][chunk_offset]);refArray14 dPdFk=reinterpret_cast<refArray14>(_rdPdF[index][0][chunk_offset]);refArray15 dk=reinterpret_cast<refArray15>(_rd[index][0][0][chunk_offset]);refArray16 system_matrixk=reinterpret_cast<refArray16>(_rsystem_matrix[index][0][chunk_offset]);
               
                Update_Position_Based_State<NEOHOOKEAN_TAG,__m512,float[16],int[16]>::Run(uk,pk,muk,mu_stabk,kappak,alphak,cutoffk,one_over_hk,cell_volumek,Uk,Vk,Sigmak,Q_hatk,dPdFk,dk,system_matrixk);
             }

        }
    };

template <int SIZE> 
class Update_Position_Based_State_MIC_COROTATED
    {
    private:
        // Generate Variables Here
        float * _local_u; float * _local_p; float * _local_mu; float * _local_mu_stab; float * _local_kappa; float * _local_alpha; float * _local_cutoff; float * _local_one_over_h; float * _local_cell_volume; float * _local_U; float * _local_V; float * _local_Sigma; float * _local_Q_hat; float * _local_dPdF; float * _local_d; float * _local_system_matrix;

    public:
        explicit Update_Position_Based_State_MIC_COROTATED(float * u_in,float * p_in,float * mu_in,float * mu_stab_in,float * kappa_in,float * alpha_in,float * cutoff_in,float * one_over_h_in,float * cell_volume_in,float * U_in,float * V_in,float * Sigma_in,float * Q_hat_in,float * dPdF_in,float * d_in,float * system_matrix_in) : _local_u(u_in),_local_p(p_in),_local_mu(mu_in),_local_mu_stab(mu_stab_in),_local_kappa(kappa_in),_local_alpha(alpha_in),_local_cutoff(cutoff_in),_local_one_over_h(one_over_h_in),_local_cell_volume(cell_volume_in),_local_U(U_in),_local_V(V_in),_local_Sigma(Sigma_in),_local_Q_hat(Q_hat_in),_local_dPdF(dPdF_in),_local_d(d_in),_local_system_matrix(system_matrix_in) {}
        void Execute(int index)
        {
        // full array typedefs
        //typedef int (&refArray)[SIZE][3][8];
        typedef float (&fullArray1)[SIZE][3][8][16];typedef float (&fullArray2)[SIZE][16];typedef float (&fullArray3)[SIZE][16];typedef float (&fullArray4)[SIZE][16];typedef float (&fullArray5)[SIZE][16];typedef float (&fullArray6)[SIZE][16];typedef float (&fullArray7)[SIZE][16];typedef float (&fullArray8)[SIZE][16];typedef float (&fullArray9)[SIZE][16];typedef float (&fullArray10)[SIZE][9][16];typedef float (&fullArray11)[SIZE][9][16];typedef float (&fullArray12)[SIZE][3][16];typedef float (&fullArray13)[SIZE][3][16];typedef float (&fullArray14)[SIZE][12][16];typedef float (&fullArray15)[SIZE][3][8][16];typedef float (&fullArray16)[SIZE][300][16];

        // full array extractions
        //refArray _rA = reinterpret_cast < refArray >(*_A);
        fullArray1 _ru =reinterpret_cast<fullArray1>(*_local_u);fullArray2 _rp =reinterpret_cast<fullArray2>(*_local_p);fullArray3 _rmu =reinterpret_cast<fullArray3>(*_local_mu);fullArray4 _rmu_stab =reinterpret_cast<fullArray4>(*_local_mu_stab);fullArray5 _rkappa =reinterpret_cast<fullArray5>(*_local_kappa);fullArray6 _ralpha =reinterpret_cast<fullArray6>(*_local_alpha);fullArray7 _rcutoff =reinterpret_cast<fullArray7>(*_local_cutoff);fullArray8 _rone_over_h =reinterpret_cast<fullArray8>(*_local_one_over_h);fullArray9 _rcell_volume =reinterpret_cast<fullArray9>(*_local_cell_volume);fullArray10 _rU =reinterpret_cast<fullArray10>(*_local_U);fullArray11 _rV =reinterpret_cast<fullArray11>(*_local_V);fullArray12 _rSigma =reinterpret_cast<fullArray12>(*_local_Sigma);fullArray13 _rQ_hat =reinterpret_cast<fullArray13>(*_local_Q_hat);fullArray14 _rdPdF =reinterpret_cast<fullArray14>(*_local_dPdF);fullArray15 _rd =reinterpret_cast<fullArray15>(*_local_d);fullArray16 _rsystem_matrix =reinterpret_cast<fullArray16>(*_local_system_matrix);

        const int ChunkSize = 16;
        // chunk typedef
        //typedef int (&refArray1)[3][8];
        typedef float (&refArray1)[3][8][16];typedef float (&refArray2)[16];typedef float (&refArray3)[16];typedef float (&refArray4)[16];typedef float (&refArray5)[16];typedef float (&refArray6)[16];typedef float (&refArray7)[16];typedef float (&refArray8)[16];typedef float (&refArray9)[16];typedef float (&refArray10)[9][16];typedef float (&refArray11)[9][16];typedef float (&refArray12)[3][16];typedef float (&refArray13)[3][16];typedef float (&refArray14)[12][16];typedef float (&refArray15)[3][8][16];typedef float (&refArray16)[300][16];

        for( int chunk_offset=0; chunk_offset<16; chunk_offset+=ChunkSize)
            {
                //refArray1 Ak = reinterpret_cast < refArray1 > (_rA[entry][0][chunk_offset]);
                refArray1 uk=reinterpret_cast<refArray1>(_ru[index][0][0][chunk_offset]);refArray2 pk=reinterpret_cast<refArray2>(_rp[index][chunk_offset]);refArray3 muk=reinterpret_cast<refArray3>(_rmu[index][chunk_offset]);refArray4 mu_stabk=reinterpret_cast<refArray4>(_rmu_stab[index][chunk_offset]);refArray5 kappak=reinterpret_cast<refArray5>(_rkappa[index][chunk_offset]);refArray6 alphak=reinterpret_cast<refArray6>(_ralpha[index][chunk_offset]);refArray7 cutoffk=reinterpret_cast<refArray7>(_rcutoff[index][chunk_offset]);refArray8 one_over_hk=reinterpret_cast<refArray8>(_rone_over_h[index][chunk_offset]);refArray9 cell_volumek=reinterpret_cast<refArray9>(_rcell_volume[index][chunk_offset]);refArray10 Uk=reinterpret_cast<refArray10>(_rU[index][0][chunk_offset]);refArray11 Vk=reinterpret_cast<refArray11>(_rV[index][0][chunk_offset]);refArray12 Sigmak=reinterpret_cast<refArray12>(_rSigma[index][0][chunk_offset]);refArray13 Q_hatk=reinterpret_cast<refArray13>(_rQ_hat[index][0][chunk_offset]);refArray14 dPdFk=reinterpret_cast<refArray14>(_rdPdF[index][0][chunk_offset]);refArray15 dk=reinterpret_cast<refArray15>(_rd[index][0][0][chunk_offset]);refArray16 system_matrixk=reinterpret_cast<refArray16>(_rsystem_matrix[index][0][chunk_offset]);
               
                Update_Position_Based_State<COROTATED_TAG,__m512,float[16],int[16]>::Run(uk,pk,muk,mu_stabk,kappak,alphak,cutoffk,one_over_hk,cell_volumek,Uk,Vk,Sigmak,Q_hatk,dPdFk,dk,system_matrixk);
             }

        }
    };

#endif

int main(int argc,char* argv[])
{
    typedef float T;

    int seed=1;
    int threads=1;
    int threads_max=1;
    int passes=1;
    const int data_size=1000000;
    if(argc>=2) threads=atoi(argv[1]);
    if(argc>=3) threads_max=atoi(argv[2]);
    if(argc>=4) passes=atoi(argv[3]);
    srand(seed);

    pthread_queue=new PTHREAD_QUEUE(threads_max);

    std::cout << "Preparing to Run " << data_size << " of all kernels with " << threads << " threads." << std::endl;

    

{
std::cout << "Running Thread Test for Update_Position_Based_State " << std::endl;

//=======================================================
//
//        DEFINE ALL VARIABLES USED BY KERNEL
//
//=======================================================
std::cout << "\nAllocating all data: ";
std::cout.flush();

start_timer();
typedef T (&u_type)[data_size][3][8][16]; u_type u = reinterpret_cast<u_type>(*((T*)(_mm_malloc(data_size*3*8*16*sizeof(T),64))));
typedef T (&p_type)[data_size][16]; p_type p = reinterpret_cast<p_type>(*((T*)(_mm_malloc(data_size*16*sizeof(T),64))));
typedef T (&mu_type)[data_size][16]; mu_type mu = reinterpret_cast<mu_type>(*((T*)(_mm_malloc(data_size*16*sizeof(T),64))));
typedef T (&mu_stab_type)[data_size][16]; mu_stab_type mu_stab = reinterpret_cast<mu_stab_type>(*((T*)(_mm_malloc(data_size*16*sizeof(T),64))));
typedef T (&kappa_type)[data_size][16]; kappa_type kappa = reinterpret_cast<kappa_type>(*((T*)(_mm_malloc(data_size*16*sizeof(T),64))));
typedef T (&alpha_type)[data_size][16]; alpha_type alpha = reinterpret_cast<alpha_type>(*((T*)(_mm_malloc(data_size*16*sizeof(T),64))));
typedef T (&cutoff_type)[data_size][16]; cutoff_type cutoff = reinterpret_cast<cutoff_type>(*((T*)(_mm_malloc(data_size*16*sizeof(T),64))));
typedef T (&one_over_h_type)[data_size][16]; one_over_h_type one_over_h = reinterpret_cast<one_over_h_type>(*((T*)(_mm_malloc(data_size*16*sizeof(T),64))));
typedef T (&cell_volume_type)[data_size][16]; cell_volume_type cell_volume = reinterpret_cast<cell_volume_type>(*((T*)(_mm_malloc(data_size*16*sizeof(T),64))));
typedef T (&U_type)[data_size][9][16]; U_type U = reinterpret_cast<U_type>(*((T*)(_mm_malloc(data_size*9*16*sizeof(T),64))));
typedef T (&V_type)[data_size][9][16]; V_type V = reinterpret_cast<V_type>(*((T*)(_mm_malloc(data_size*9*16*sizeof(T),64))));
typedef T (&Sigma_type)[data_size][3][16]; Sigma_type Sigma = reinterpret_cast<Sigma_type>(*((T*)(_mm_malloc(data_size*3*16*sizeof(T),64))));
typedef T (&Q_hat_type)[data_size][3][16]; Q_hat_type Q_hat = reinterpret_cast<Q_hat_type>(*((T*)(_mm_malloc(data_size*3*16*sizeof(T),64))));
typedef T (&dPdF_type)[data_size][12][16]; dPdF_type dPdF = reinterpret_cast<dPdF_type>(*((T*)(_mm_malloc(data_size*12*16*sizeof(T),64))));
typedef T (&d_type)[data_size][3][8][16]; d_type d = reinterpret_cast<d_type>(*((T*)(_mm_malloc(data_size*3*8*16*sizeof(T),64))));
typedef T (&system_matrix_type)[data_size][300][16]; system_matrix_type system_matrix = reinterpret_cast<system_matrix_type>(*((T*)(_mm_malloc(data_size*300*16*sizeof(T),64))));
  

for(int __a=0;__a<data_size;__a++) for(int __b=0;__b<3;__b++) for(int __c=0;__c<8;__c++) for(int __d=0;__d<16;__d++){ 
u[__a][__b][__c][__d]=Get_Random<float>();}for(int __a=0;__a<data_size;__a++) for(int __b=0;__b<16;__b++){ 
p[__a][__b]=Get_Random<float>();}for(int __a=0;__a<data_size;__a++) for(int __b=0;__b<16;__b++){ 
mu[__a][__b]=Get_Random<float>();}for(int __a=0;__a<data_size;__a++) for(int __b=0;__b<16;__b++){ 
mu_stab[__a][__b]=Get_Random<float>();}for(int __a=0;__a<data_size;__a++) for(int __b=0;__b<16;__b++){ 
kappa[__a][__b]=Get_Random<float>();}for(int __a=0;__a<data_size;__a++) for(int __b=0;__b<16;__b++){ 
alpha[__a][__b]=Get_Random<float>();}for(int __a=0;__a<data_size;__a++) for(int __b=0;__b<16;__b++){ 
cutoff[__a][__b]=Get_Random<float>();}for(int __a=0;__a<data_size;__a++) for(int __b=0;__b<16;__b++){ 
one_over_h[__a][__b]=Get_Random<float>();}for(int __a=0;__a<data_size;__a++) for(int __b=0;__b<16;__b++){ 
cell_volume[__a][__b]=Get_Random<float>();}for(int __a=0;__a<data_size;__a++) for(int __b=0;__b<9;__b++) for(int __c=0;__c<16;__c++){ 
U[__a][__b][__c]=Get_Random<float>();}for(int __a=0;__a<data_size;__a++) for(int __b=0;__b<9;__b++) for(int __c=0;__c<16;__c++){ 
V[__a][__b][__c]=Get_Random<float>();}for(int __a=0;__a<data_size;__a++) for(int __b=0;__b<3;__b++) for(int __c=0;__c<16;__c++){ 
Sigma[__a][__b][__c]=Get_Random<float>();}for(int __a=0;__a<data_size;__a++) for(int __b=0;__b<3;__b++) for(int __c=0;__c<16;__c++){ 
Q_hat[__a][__b][__c]=Get_Random<float>();}for(int __a=0;__a<data_size;__a++) for(int __b=0;__b<12;__b++) for(int __c=0;__c<16;__c++){ 
dPdF[__a][__b][__c]=Get_Random<float>();}for(int __a=0;__a<data_size;__a++) for(int __b=0;__b<3;__b++) for(int __c=0;__c<8;__c++) for(int __d=0;__d<16;__d++){ 
d[__a][__b][__c][__d]=Get_Random<float>();}for(int __a=0;__a<data_size;__a++) for(int __b=0;__b<300;__b++) for(int __c=0;__c<16;__c++){ 
system_matrix[__a][__b][__c]=Get_Random<float>();}
stop_timer();

std::cout << get_time() << "s\n\n"<< std::endl;


//=======================================================
//
//             COMPUTE SCALAR RESULTS
//
//=======================================================

{
std::cout << "	Running " << data_size << " of SCALAR :  "<< std::endl;


    Update_Position_Based_State_SCALAR_NEOHOOKEAN<data_size> op((float *)&u,(float *)&p,(float *)&mu,(float *)&mu_stab,(float *)&kappa,(float *)&alpha,(float *)&cutoff,(float *)&one_over_h,(float *)&cell_volume,(float *)&U,(float *)&V,(float *)&Sigma,(float *)&Q_hat,(float *)&dPdF,(float *)&d,(float *)&system_matrix);

    for( int t = threads; t <= threads_max; t+=std::max<int>(((threads_max - threads)/30), 1)){
         std::cout << "Running Test with " << t << " threads." << std::endl;
         MT_Streaming_Kernels::Kernel_Serial_Base_Helper<Update_Position_Based_State_SCALAR_NEOHOOKEAN<data_size> > helper(op,data_size,t);

         double min_time = 10000000;
         double max_time = -1;
         double avg_time = 0;
                  
         for(int i=0; i<passes; i++){   
           start_timer();
           helper.Run_Parallel();
           stop_timer();
           std::cout << get_time() << "s"<< std::endl;
           min_time = std::min<double>( min_time, get_time() );
           max_time = std::max<double>( max_time, get_time() );
           avg_time += get_time();
         }
         avg_time = avg_time / passes;
         std::cout << "Min pass time: " << min_time << std::endl;
         std::cout << "Max pass time: " << max_time << std::endl;
         std::cout << "Avg pass time: " << avg_time << std::endl;
     }




}

//=======================================================
//
//             COMPUTE SSE RESULTS
//
//=======================================================
#ifdef ENABLE_SSE_INSTRUCTION_SET
{
std::cout << "	Running " << data_size << " of SSE :  "<< std::endl;


    Update_Position_Based_State_SSE_NEOHOOKEAN<data_size> op((float *)&u,(float *)&p,(float *)&mu,(float *)&mu_stab,(float *)&kappa,(float *)&alpha,(float *)&cutoff,(float *)&one_over_h,(float *)&cell_volume,(float *)&U,(float *)&V,(float *)&Sigma,(float *)&Q_hat,(float *)&dPdF,(float *)&d,(float *)&system_matrix);

    for( int t = threads; t <= threads_max; t+=std::max<int>(((threads_max - threads)/30), 1)){
         std::cout << "Running Test with " << t << " threads." << std::endl;
         MT_Streaming_Kernels::Kernel_Serial_Base_Helper<Update_Position_Based_State_SSE_NEOHOOKEAN<data_size> > helper(op,data_size,t);

         double min_time = 10000000;
         double max_time = -1;
         double avg_time = 0;
                  
         for(int i=0; i<passes; i++){   
           start_timer();
           helper.Run_Parallel();
           stop_timer();
           std::cout << get_time() << "s"<< std::endl;
           min_time = std::min<double>( min_time, get_time() );
           max_time = std::max<double>( max_time, get_time() );
           avg_time += get_time();
         }
         avg_time = avg_time / passes;
         std::cout << "Min pass time: " << min_time << std::endl;
         std::cout << "Max pass time: " << max_time << std::endl;
         std::cout << "Avg pass time: " << avg_time << std::endl;
     }




}
#endif

//=======================================================
//
//             COMPUTE AVX RESULTS
//
//=======================================================
#ifdef ENABLE_AVX_INSTRUCTION_SET
{
std::cout << "	Running " << data_size << " of AVX :  "<< std::endl;


    Update_Position_Based_State_AVX_NEOHOOKEAN<data_size> op((float *)&u,(float *)&p,(float *)&mu,(float *)&mu_stab,(float *)&kappa,(float *)&alpha,(float *)&cutoff,(float *)&one_over_h,(float *)&cell_volume,(float *)&U,(float *)&V,(float *)&Sigma,(float *)&Q_hat,(float *)&dPdF,(float *)&d,(float *)&system_matrix);

    for( int t = threads; t <= threads_max; t+=std::max<int>(((threads_max - threads)/30), 1)){
         std::cout << "Running Test with " << t << " threads." << std::endl;
         MT_Streaming_Kernels::Kernel_Serial_Base_Helper<Update_Position_Based_State_AVX_NEOHOOKEAN<data_size> > helper(op,data_size,t);

         double min_time = 10000000;
         double max_time = -1;
         double avg_time = 0;
                  
         for(int i=0; i<passes; i++){   
           start_timer();
           helper.Run_Parallel();
           stop_timer();
           std::cout << get_time() << "s"<< std::endl;
           min_time = std::min<double>( min_time, get_time() );
           max_time = std::max<double>( max_time, get_time() );
           avg_time += get_time();
         }
         avg_time = avg_time / passes;
         std::cout << "Min pass time: " << min_time << std::endl;
         std::cout << "Max pass time: " << max_time << std::endl;
         std::cout << "Avg pass time: " << avg_time << std::endl;
     }




}
#endif

//=======================================================
//
//             COMPUTE NEON RESULTS
//
//=======================================================
#ifdef ENABLE_NEON_INSTRUCTION_SET
{
std::cout << "	Running " << data_size << " of NEON :  "<< std::endl;


    Update_Position_Based_State_NEON_NEOHOOKEAN<data_size> op((float *)&u,(float *)&p,(float *)&mu,(float *)&mu_stab,(float *)&kappa,(float *)&alpha,(float *)&cutoff,(float *)&one_over_h,(float *)&cell_volume,(float *)&U,(float *)&V,(float *)&Sigma,(float *)&Q_hat,(float *)&dPdF,(float *)&d,(float *)&system_matrix);

    for( int t = threads; t <= threads_max; t+=std::max<int>(((threads_max - threads)/30), 1)){
         std::cout << "Running Test with " << t << " threads." << std::endl;
         MT_Streaming_Kernels::Kernel_Serial_Base_Helper<Update_Position_Based_State_NEON_NEOHOOKEAN<data_size> > helper(op,data_size,t);

         double min_time = 10000000;
         double max_time = -1;
         double avg_time = 0;
                  
         for(int i=0; i<passes; i++){   
           start_timer();
           helper.Run_Parallel();
           stop_timer();
           std::cout << get_time() << "s"<< std::endl;
           min_time = std::min<double>( min_time, get_time() );
           max_time = std::max<double>( max_time, get_time() );
           avg_time += get_time();
         }
         avg_time = avg_time / passes;
         std::cout << "Min pass time: " << min_time << std::endl;
         std::cout << "Max pass time: " << max_time << std::endl;
         std::cout << "Avg pass time: " << avg_time << std::endl;
     }




}
#endif

//=======================================================
//
//             COMPUTE MIC RESULTS
//
//=======================================================
#ifdef ENABLE_MIC_INSTRUCTION_SET
{
std::cout << "	Running " << data_size << " of MIC :  "<< std::endl;


    Update_Position_Based_State_MIC_NEOHOOKEAN<data_size> op((float *)&u,(float *)&p,(float *)&mu,(float *)&mu_stab,(float *)&kappa,(float *)&alpha,(float *)&cutoff,(float *)&one_over_h,(float *)&cell_volume,(float *)&U,(float *)&V,(float *)&Sigma,(float *)&Q_hat,(float *)&dPdF,(float *)&d,(float *)&system_matrix);

    for( int t = threads; t <= threads_max; t+=std::max<int>(((threads_max - threads)/30), 1)){
         std::cout << "Running Test with " << t << " threads." << std::endl;
         MT_Streaming_Kernels::Kernel_Serial_Base_Helper<Update_Position_Based_State_MIC_NEOHOOKEAN<data_size> > helper(op,data_size,t);

         double min_time = 10000000;
         double max_time = -1;
         double avg_time = 0;
                  
         for(int i=0; i<passes; i++){   
           start_timer();
           helper.Run_Parallel();
           stop_timer();
           std::cout << get_time() << "s"<< std::endl;
           min_time = std::min<double>( min_time, get_time() );
           max_time = std::max<double>( max_time, get_time() );
           avg_time += get_time();
         }
         avg_time = avg_time / passes;
         std::cout << "Min pass time: " << min_time << std::endl;
         std::cout << "Max pass time: " << max_time << std::endl;
         std::cout << "Avg pass time: " << avg_time << std::endl;
     }




}
#endif

//=======================================================
//
//        FREE MEMORY USED BY ALL VARIABLES
//
//=======================================================
std::cout << "\nFreeing all data: " << std::endl;
std::cout.flush();

_mm_free( reinterpret_cast< void* >( u ));
_mm_free( reinterpret_cast< void* >( p ));
_mm_free( reinterpret_cast< void* >( mu ));
_mm_free( reinterpret_cast< void* >( mu_stab ));
_mm_free( reinterpret_cast< void* >( kappa ));
_mm_free( reinterpret_cast< void* >( alpha ));
_mm_free( reinterpret_cast< void* >( cutoff ));
_mm_free( reinterpret_cast< void* >( one_over_h ));
_mm_free( reinterpret_cast< void* >( cell_volume ));
_mm_free( reinterpret_cast< void* >( U ));
_mm_free( reinterpret_cast< void* >( V ));
_mm_free( reinterpret_cast< void* >( Sigma ));
_mm_free( reinterpret_cast< void* >( Q_hat ));
_mm_free( reinterpret_cast< void* >( dPdF ));
_mm_free( reinterpret_cast< void* >( d ));
_mm_free( reinterpret_cast< void* >( system_matrix ));


}


{
std::cout << "Running Thread Test for Update_Position_Based_State " << std::endl;

//=======================================================
//
//        DEFINE ALL VARIABLES USED BY KERNEL
//
//=======================================================
std::cout << "\nAllocating all data: ";
std::cout.flush();

start_timer();
typedef T (&u_type)[data_size][3][8][16]; u_type u = reinterpret_cast<u_type>(*((T*)(_mm_malloc(data_size*3*8*16*sizeof(T),64))));
typedef T (&p_type)[data_size][16]; p_type p = reinterpret_cast<p_type>(*((T*)(_mm_malloc(data_size*16*sizeof(T),64))));
typedef T (&mu_type)[data_size][16]; mu_type mu = reinterpret_cast<mu_type>(*((T*)(_mm_malloc(data_size*16*sizeof(T),64))));
typedef T (&mu_stab_type)[data_size][16]; mu_stab_type mu_stab = reinterpret_cast<mu_stab_type>(*((T*)(_mm_malloc(data_size*16*sizeof(T),64))));
typedef T (&kappa_type)[data_size][16]; kappa_type kappa = reinterpret_cast<kappa_type>(*((T*)(_mm_malloc(data_size*16*sizeof(T),64))));
typedef T (&alpha_type)[data_size][16]; alpha_type alpha = reinterpret_cast<alpha_type>(*((T*)(_mm_malloc(data_size*16*sizeof(T),64))));
typedef T (&cutoff_type)[data_size][16]; cutoff_type cutoff = reinterpret_cast<cutoff_type>(*((T*)(_mm_malloc(data_size*16*sizeof(T),64))));
typedef T (&one_over_h_type)[data_size][16]; one_over_h_type one_over_h = reinterpret_cast<one_over_h_type>(*((T*)(_mm_malloc(data_size*16*sizeof(T),64))));
typedef T (&cell_volume_type)[data_size][16]; cell_volume_type cell_volume = reinterpret_cast<cell_volume_type>(*((T*)(_mm_malloc(data_size*16*sizeof(T),64))));
typedef T (&U_type)[data_size][9][16]; U_type U = reinterpret_cast<U_type>(*((T*)(_mm_malloc(data_size*9*16*sizeof(T),64))));
typedef T (&V_type)[data_size][9][16]; V_type V = reinterpret_cast<V_type>(*((T*)(_mm_malloc(data_size*9*16*sizeof(T),64))));
typedef T (&Sigma_type)[data_size][3][16]; Sigma_type Sigma = reinterpret_cast<Sigma_type>(*((T*)(_mm_malloc(data_size*3*16*sizeof(T),64))));
typedef T (&Q_hat_type)[data_size][3][16]; Q_hat_type Q_hat = reinterpret_cast<Q_hat_type>(*((T*)(_mm_malloc(data_size*3*16*sizeof(T),64))));
typedef T (&dPdF_type)[data_size][12][16]; dPdF_type dPdF = reinterpret_cast<dPdF_type>(*((T*)(_mm_malloc(data_size*12*16*sizeof(T),64))));
typedef T (&d_type)[data_size][3][8][16]; d_type d = reinterpret_cast<d_type>(*((T*)(_mm_malloc(data_size*3*8*16*sizeof(T),64))));
typedef T (&system_matrix_type)[data_size][300][16]; system_matrix_type system_matrix = reinterpret_cast<system_matrix_type>(*((T*)(_mm_malloc(data_size*300*16*sizeof(T),64))));
  

for(int __a=0;__a<data_size;__a++) for(int __b=0;__b<3;__b++) for(int __c=0;__c<8;__c++) for(int __d=0;__d<16;__d++){ 
u[__a][__b][__c][__d]=Get_Random<float>();}for(int __a=0;__a<data_size;__a++) for(int __b=0;__b<16;__b++){ 
p[__a][__b]=Get_Random<float>();}for(int __a=0;__a<data_size;__a++) for(int __b=0;__b<16;__b++){ 
mu[__a][__b]=Get_Random<float>();}for(int __a=0;__a<data_size;__a++) for(int __b=0;__b<16;__b++){ 
mu_stab[__a][__b]=Get_Random<float>();}for(int __a=0;__a<data_size;__a++) for(int __b=0;__b<16;__b++){ 
kappa[__a][__b]=Get_Random<float>();}for(int __a=0;__a<data_size;__a++) for(int __b=0;__b<16;__b++){ 
alpha[__a][__b]=Get_Random<float>();}for(int __a=0;__a<data_size;__a++) for(int __b=0;__b<16;__b++){ 
cutoff[__a][__b]=Get_Random<float>();}for(int __a=0;__a<data_size;__a++) for(int __b=0;__b<16;__b++){ 
one_over_h[__a][__b]=Get_Random<float>();}for(int __a=0;__a<data_size;__a++) for(int __b=0;__b<16;__b++){ 
cell_volume[__a][__b]=Get_Random<float>();}for(int __a=0;__a<data_size;__a++) for(int __b=0;__b<9;__b++) for(int __c=0;__c<16;__c++){ 
U[__a][__b][__c]=Get_Random<float>();}for(int __a=0;__a<data_size;__a++) for(int __b=0;__b<9;__b++) for(int __c=0;__c<16;__c++){ 
V[__a][__b][__c]=Get_Random<float>();}for(int __a=0;__a<data_size;__a++) for(int __b=0;__b<3;__b++) for(int __c=0;__c<16;__c++){ 
Sigma[__a][__b][__c]=Get_Random<float>();}for(int __a=0;__a<data_size;__a++) for(int __b=0;__b<3;__b++) for(int __c=0;__c<16;__c++){ 
Q_hat[__a][__b][__c]=Get_Random<float>();}for(int __a=0;__a<data_size;__a++) for(int __b=0;__b<12;__b++) for(int __c=0;__c<16;__c++){ 
dPdF[__a][__b][__c]=Get_Random<float>();}for(int __a=0;__a<data_size;__a++) for(int __b=0;__b<3;__b++) for(int __c=0;__c<8;__c++) for(int __d=0;__d<16;__d++){ 
d[__a][__b][__c][__d]=Get_Random<float>();}for(int __a=0;__a<data_size;__a++) for(int __b=0;__b<300;__b++) for(int __c=0;__c<16;__c++){ 
system_matrix[__a][__b][__c]=Get_Random<float>();}
stop_timer();

std::cout << get_time() << "s\n\n"<< std::endl;


//=======================================================
//
//             COMPUTE SCALAR RESULTS
//
//=======================================================

{
std::cout << "	Running " << data_size << " of SCALAR :  "<< std::endl;


    Update_Position_Based_State_SCALAR_COROTATED<data_size> op((float *)&u,(float *)&p,(float *)&mu,(float *)&mu_stab,(float *)&kappa,(float *)&alpha,(float *)&cutoff,(float *)&one_over_h,(float *)&cell_volume,(float *)&U,(float *)&V,(float *)&Sigma,(float *)&Q_hat,(float *)&dPdF,(float *)&d,(float *)&system_matrix);

    for( int t = threads; t <= threads_max; t+=std::max<int>(((threads_max - threads)/30), 1)){
         std::cout << "Running Test with " << t << " threads." << std::endl;
         MT_Streaming_Kernels::Kernel_Serial_Base_Helper<Update_Position_Based_State_SCALAR_COROTATED<data_size> > helper(op,data_size,t);

         double min_time = 10000000;
         double max_time = -1;
         double avg_time = 0;
                  
         for(int i=0; i<passes; i++){   
           start_timer();
           helper.Run_Parallel();
           stop_timer();
           std::cout << get_time() << "s"<< std::endl;
           min_time = std::min<double>( min_time, get_time() );
           max_time = std::max<double>( max_time, get_time() );
           avg_time += get_time();
         }
         avg_time = avg_time / passes;
         std::cout << "Min pass time: " << min_time << std::endl;
         std::cout << "Max pass time: " << max_time << std::endl;
         std::cout << "Avg pass time: " << avg_time << std::endl;
     }




}

//=======================================================
//
//             COMPUTE SSE RESULTS
//
//=======================================================
#ifdef ENABLE_SSE_INSTRUCTION_SET
{
std::cout << "	Running " << data_size << " of SSE :  "<< std::endl;


    Update_Position_Based_State_SSE_COROTATED<data_size> op((float *)&u,(float *)&p,(float *)&mu,(float *)&mu_stab,(float *)&kappa,(float *)&alpha,(float *)&cutoff,(float *)&one_over_h,(float *)&cell_volume,(float *)&U,(float *)&V,(float *)&Sigma,(float *)&Q_hat,(float *)&dPdF,(float *)&d,(float *)&system_matrix);

    for( int t = threads; t <= threads_max; t+=std::max<int>(((threads_max - threads)/30), 1)){
         std::cout << "Running Test with " << t << " threads." << std::endl;
         MT_Streaming_Kernels::Kernel_Serial_Base_Helper<Update_Position_Based_State_SSE_COROTATED<data_size> > helper(op,data_size,t);

         double min_time = 10000000;
         double max_time = -1;
         double avg_time = 0;
                  
         for(int i=0; i<passes; i++){   
           start_timer();
           helper.Run_Parallel();
           stop_timer();
           std::cout << get_time() << "s"<< std::endl;
           min_time = std::min<double>( min_time, get_time() );
           max_time = std::max<double>( max_time, get_time() );
           avg_time += get_time();
         }
         avg_time = avg_time / passes;
         std::cout << "Min pass time: " << min_time << std::endl;
         std::cout << "Max pass time: " << max_time << std::endl;
         std::cout << "Avg pass time: " << avg_time << std::endl;
     }




}
#endif

//=======================================================
//
//             COMPUTE AVX RESULTS
//
//=======================================================
#ifdef ENABLE_AVX_INSTRUCTION_SET
{
std::cout << "	Running " << data_size << " of AVX :  "<< std::endl;


    Update_Position_Based_State_AVX_COROTATED<data_size> op((float *)&u,(float *)&p,(float *)&mu,(float *)&mu_stab,(float *)&kappa,(float *)&alpha,(float *)&cutoff,(float *)&one_over_h,(float *)&cell_volume,(float *)&U,(float *)&V,(float *)&Sigma,(float *)&Q_hat,(float *)&dPdF,(float *)&d,(float *)&system_matrix);

    for( int t = threads; t <= threads_max; t+=std::max<int>(((threads_max - threads)/30), 1)){
         std::cout << "Running Test with " << t << " threads." << std::endl;
         MT_Streaming_Kernels::Kernel_Serial_Base_Helper<Update_Position_Based_State_AVX_COROTATED<data_size> > helper(op,data_size,t);

         double min_time = 10000000;
         double max_time = -1;
         double avg_time = 0;
                  
         for(int i=0; i<passes; i++){   
           start_timer();
           helper.Run_Parallel();
           stop_timer();
           std::cout << get_time() << "s"<< std::endl;
           min_time = std::min<double>( min_time, get_time() );
           max_time = std::max<double>( max_time, get_time() );
           avg_time += get_time();
         }
         avg_time = avg_time / passes;
         std::cout << "Min pass time: " << min_time << std::endl;
         std::cout << "Max pass time: " << max_time << std::endl;
         std::cout << "Avg pass time: " << avg_time << std::endl;
     }




}
#endif

//=======================================================
//
//             COMPUTE NEON RESULTS
//
//=======================================================
#ifdef ENABLE_NEON_INSTRUCTION_SET
{
std::cout << "	Running " << data_size << " of NEON :  "<< std::endl;


    Update_Position_Based_State_NEON_COROTATED<data_size> op((float *)&u,(float *)&p,(float *)&mu,(float *)&mu_stab,(float *)&kappa,(float *)&alpha,(float *)&cutoff,(float *)&one_over_h,(float *)&cell_volume,(float *)&U,(float *)&V,(float *)&Sigma,(float *)&Q_hat,(float *)&dPdF,(float *)&d,(float *)&system_matrix);

    for( int t = threads; t <= threads_max; t+=std::max<int>(((threads_max - threads)/30), 1)){
         std::cout << "Running Test with " << t << " threads." << std::endl;
         MT_Streaming_Kernels::Kernel_Serial_Base_Helper<Update_Position_Based_State_NEON_COROTATED<data_size> > helper(op,data_size,t);

         double min_time = 10000000;
         double max_time = -1;
         double avg_time = 0;
                  
         for(int i=0; i<passes; i++){   
           start_timer();
           helper.Run_Parallel();
           stop_timer();
           std::cout << get_time() << "s"<< std::endl;
           min_time = std::min<double>( min_time, get_time() );
           max_time = std::max<double>( max_time, get_time() );
           avg_time += get_time();
         }
         avg_time = avg_time / passes;
         std::cout << "Min pass time: " << min_time << std::endl;
         std::cout << "Max pass time: " << max_time << std::endl;
         std::cout << "Avg pass time: " << avg_time << std::endl;
     }




}
#endif

//=======================================================
//
//             COMPUTE MIC RESULTS
//
//=======================================================
#ifdef ENABLE_MIC_INSTRUCTION_SET
{
std::cout << "	Running " << data_size << " of MIC :  "<< std::endl;


    Update_Position_Based_State_MIC_COROTATED<data_size> op((float *)&u,(float *)&p,(float *)&mu,(float *)&mu_stab,(float *)&kappa,(float *)&alpha,(float *)&cutoff,(float *)&one_over_h,(float *)&cell_volume,(float *)&U,(float *)&V,(float *)&Sigma,(float *)&Q_hat,(float *)&dPdF,(float *)&d,(float *)&system_matrix);

    for( int t = threads; t <= threads_max; t+=std::max<int>(((threads_max - threads)/30), 1)){
         std::cout << "Running Test with " << t << " threads." << std::endl;
         MT_Streaming_Kernels::Kernel_Serial_Base_Helper<Update_Position_Based_State_MIC_COROTATED<data_size> > helper(op,data_size,t);

         double min_time = 10000000;
         double max_time = -1;
         double avg_time = 0;
                  
         for(int i=0; i<passes; i++){   
           start_timer();
           helper.Run_Parallel();
           stop_timer();
           std::cout << get_time() << "s"<< std::endl;
           min_time = std::min<double>( min_time, get_time() );
           max_time = std::max<double>( max_time, get_time() );
           avg_time += get_time();
         }
         avg_time = avg_time / passes;
         std::cout << "Min pass time: " << min_time << std::endl;
         std::cout << "Max pass time: " << max_time << std::endl;
         std::cout << "Avg pass time: " << avg_time << std::endl;
     }




}
#endif

//=======================================================
//
//        FREE MEMORY USED BY ALL VARIABLES
//
//=======================================================
std::cout << "\nFreeing all data: " << std::endl;
std::cout.flush();

_mm_free( reinterpret_cast< void* >( u ));
_mm_free( reinterpret_cast< void* >( p ));
_mm_free( reinterpret_cast< void* >( mu ));
_mm_free( reinterpret_cast< void* >( mu_stab ));
_mm_free( reinterpret_cast< void* >( kappa ));
_mm_free( reinterpret_cast< void* >( alpha ));
_mm_free( reinterpret_cast< void* >( cutoff ));
_mm_free( reinterpret_cast< void* >( one_over_h ));
_mm_free( reinterpret_cast< void* >( cell_volume ));
_mm_free( reinterpret_cast< void* >( U ));
_mm_free( reinterpret_cast< void* >( V ));
_mm_free( reinterpret_cast< void* >( Sigma ));
_mm_free( reinterpret_cast< void* >( Q_hat ));
_mm_free( reinterpret_cast< void* >( dPdF ));
_mm_free( reinterpret_cast< void* >( d ));
_mm_free( reinterpret_cast< void* >( system_matrix ));


}


    return 0;

}




#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <sys/time.h>
#include "KernelCommon.h"
struct NEOHOOKEAN_TAG;
struct COROTATED_TAG;

#include "Add_Force_Differential.h"

#include <Thread_Queueing/PTHREAD_QUEUE.h>
#include <Kernel_Serial_Base_Helper.h>


template < class T > T Get_Random (const T a = (T) - 1., const T b = (T) 1.)
{
    return 3.14;
  return ((b - a) * (T) rand ()) / (T) RAND_MAX + a;
}

struct timeval starttime, stoptime;
void
start_timer ()
{
  gettimeofday (&starttime, NULL);
}

void
stop_timer ()
{
  gettimeofday (&stoptime, NULL);
}

double
get_time ()
{
  return (double) stoptime.tv_sec - (double) starttime.tv_sec +
    (double) 1e-6 *(double) stoptime.tv_usec -
    (double) 1e-6 *(double) starttime.tv_usec;
}


template < int SIZE > class Add_Force_Differential_SCALAR_None
{
private:
  // Generate Variables Here
  float *_local_du;
  float *_local_dp;
  float *_local_alpha_squared_over_kappa;
  float *_local_alpha;
  float *_local_one_over_h;
  float *_local_cell_volume;
  float *_local_Q_hat;
  float *_local_U;
  float *_local_V;
  float *_local_dPdF;
  float *_local_df;
  float *_local_dq;

public:
    explicit Add_Force_Differential_SCALAR_None (float *du_in, float *dp_in,
                                                 float
                                                 *alpha_squared_over_kappa_in,
                                                 float *alpha_in,
                                                 float *one_over_h_in,
                                                 float *cell_volume_in,
                                                 float *Q_hat_in, float *U_in,
                                                 float *V_in, float *dPdF_in,
                                                 float *df_in,
                                                 float
                                                 *dq_in):_local_du (du_in),
    _local_dp (dp_in),
    _local_alpha_squared_over_kappa (alpha_squared_over_kappa_in),
    _local_alpha (alpha_in), _local_one_over_h (one_over_h_in),
    _local_cell_volume (cell_volume_in), _local_Q_hat (Q_hat_in),
    _local_U (U_in), _local_V (V_in), _local_dPdF (dPdF_in), _local_df (df_in),
    _local_dq (dq_in)
  {
  }
  void Execute (int index)
  {
    // full array typedefs
    //typedef int (&refArray)[SIZE][3][8];
    typedef float (&fullArray1)[SIZE][3][8][16];
    typedef float (&fullArray2)[SIZE][16];
    typedef float (&fullArray3)[SIZE][16];
    typedef float (&fullArray4)[SIZE][16];
    typedef float (&fullArray5)[SIZE][16];
    typedef float (&fullArray6)[SIZE][16];
    typedef float (&fullArray7)[SIZE][3][16];
    typedef float (&fullArray8)[SIZE][9][16];
    typedef float (&fullArray9)[SIZE][9][16];
    typedef float (&fullArray10)[SIZE][12][16];
    typedef float (&fullArray11)[SIZE][3][8][16];
    typedef float (&fullArray12)[SIZE][16];

    // full array extractions
    //refArray _rA = reinterpret_cast < refArray >(*_A);
    fullArray1 _rdu = reinterpret_cast < fullArray1 > (*_local_du);
    fullArray2 _rdp = reinterpret_cast < fullArray2 > (*_local_dp);
    fullArray3 _ralpha_squared_over_kappa =
      reinterpret_cast < fullArray3 > (*_local_alpha_squared_over_kappa);
    fullArray4 _ralpha = reinterpret_cast < fullArray4 > (*_local_alpha);
    fullArray5 _rone_over_h =
      reinterpret_cast < fullArray5 > (*_local_one_over_h);
    fullArray6 _rcell_volume =
      reinterpret_cast < fullArray6 > (*_local_cell_volume);
    fullArray7 _rQ_hat = reinterpret_cast < fullArray7 > (*_local_Q_hat);
    fullArray8 _rU = reinterpret_cast < fullArray8 > (*_local_U);
    fullArray9 _rV = reinterpret_cast < fullArray9 > (*_local_V);
    fullArray10 _rdPdF = reinterpret_cast < fullArray10 > (*_local_dPdF);
    fullArray11 _rdf = reinterpret_cast < fullArray11 > (*_local_df);
    fullArray12 _rdq = reinterpret_cast < fullArray12 > (*_local_dq);

    const int ChunkSize = 1;
    // chunk typedef
    //typedef int (&refArray1)[3][8];
    typedef float (&refArray1)[3][8][16];
    typedef float (&refArray2)[16];
    typedef float (&refArray3)[16];
    typedef float (&refArray4)[16];
    typedef float (&refArray5)[16];
    typedef float (&refArray6)[16];
    typedef float (&refArray7)[3][16];
    typedef float (&refArray8)[9][16];
    typedef float (&refArray9)[9][16];
    typedef float (&refArray10)[12][16];
    typedef float (&refArray11)[3][8][16];
    typedef float (&refArray12)[16];

    for (int chunk_offset = 0; chunk_offset < 16; chunk_offset += ChunkSize)
      {
        //refArray1 Ak = reinterpret_cast < refArray1 > (_rA[entry][0][chunk_offset]);
        refArray1 duk =
          reinterpret_cast < refArray1 > (_rdu[index][0][0][chunk_offset]);
        refArray2 dpk =
          reinterpret_cast < refArray2 > (_rdp[index][chunk_offset]);
        refArray3 alpha_squared_over_kappak =
          reinterpret_cast < refArray3 >
          (_ralpha_squared_over_kappa[index][chunk_offset]);
        refArray4 alphak =
          reinterpret_cast < refArray4 > (_ralpha[index][chunk_offset]);
        refArray5 one_over_hk =
          reinterpret_cast < refArray5 > (_rone_over_h[index][chunk_offset]);
        refArray6 cell_volumek =
          reinterpret_cast < refArray6 > (_rcell_volume[index][chunk_offset]);
        refArray7 Q_hatk =
          reinterpret_cast < refArray7 > (_rQ_hat[index][0][chunk_offset]);
        refArray8 Uk =
          reinterpret_cast < refArray8 > (_rU[index][0][chunk_offset]);
        refArray9 Vk =
          reinterpret_cast < refArray9 > (_rV[index][0][chunk_offset]);
        refArray10 dPdFk =
          reinterpret_cast < refArray10 > (_rdPdF[index][0][chunk_offset]);
        refArray11 dfk =
          reinterpret_cast < refArray11 > (_rdf[index][0][0][chunk_offset]);
        refArray12 dqk =
          reinterpret_cast < refArray12 > (_rdq[index][chunk_offset]);

        Add_Force_Differential < float, float[16], int[16] > (duk, dpk,
                                                              alpha_squared_over_kappak,
                                                              alphak,
                                                              one_over_hk,
                                                              cell_volumek,
                                                              Q_hatk, Uk, Vk,
                                                              dPdFk, dfk, dqk);
      }

  }
};


#ifdef ENABLE_SSE_INSTRUCTION_SET

template < int SIZE > class Add_Force_Differential_SSE_None
{
private:
  // Generate Variables Here
  float *_local_du;
  float *_local_dp;
  float *_local_alpha_squared_over_kappa;
  float *_local_alpha;
  float *_local_one_over_h;
  float *_local_cell_volume;
  float *_local_Q_hat;
  float *_local_U;
  float *_local_V;
  float *_local_dPdF;
  float *_local_df;
  float *_local_dq;

public:
    explicit Add_Force_Differential_SSE_None (float *du_in, float *dp_in,
                                              float
                                              *alpha_squared_over_kappa_in,
                                              float *alpha_in,
                                              float *one_over_h_in,
                                              float *cell_volume_in,
                                              float *Q_hat_in, float *U_in,
                                              float *V_in, float *dPdF_in,
                                              float *df_in,
                                              float *dq_in):_local_du (du_in),
    _local_dp (dp_in),
    _local_alpha_squared_over_kappa (alpha_squared_over_kappa_in),
    _local_alpha (alpha_in), _local_one_over_h (one_over_h_in),
    _local_cell_volume (cell_volume_in), _local_Q_hat (Q_hat_in),
    _local_U (U_in), _local_V (V_in), _local_dPdF (dPdF_in), _local_df (df_in),
    _local_dq (dq_in)
  {
  }
  void Execute (int index)
  {
    // full array typedefs
    //typedef int (&refArray)[SIZE][3][8];
    typedef float (&fullArray1)[SIZE][3][8][16];
    typedef float (&fullArray2)[SIZE][16];
    typedef float (&fullArray3)[SIZE][16];
    typedef float (&fullArray4)[SIZE][16];
    typedef float (&fullArray5)[SIZE][16];
    typedef float (&fullArray6)[SIZE][16];
    typedef float (&fullArray7)[SIZE][3][16];
    typedef float (&fullArray8)[SIZE][9][16];
    typedef float (&fullArray9)[SIZE][9][16];
    typedef float (&fullArray10)[SIZE][12][16];
    typedef float (&fullArray11)[SIZE][3][8][16];
    typedef float (&fullArray12)[SIZE][16];

    // full array extractions
    //refArray _rA = reinterpret_cast < refArray >(*_A);
    fullArray1 _rdu = reinterpret_cast < fullArray1 > (*_local_du);
    fullArray2 _rdp = reinterpret_cast < fullArray2 > (*_local_dp);
    fullArray3 _ralpha_squared_over_kappa =
      reinterpret_cast < fullArray3 > (*_local_alpha_squared_over_kappa);
    fullArray4 _ralpha = reinterpret_cast < fullArray4 > (*_local_alpha);
    fullArray5 _rone_over_h =
      reinterpret_cast < fullArray5 > (*_local_one_over_h);
    fullArray6 _rcell_volume =
      reinterpret_cast < fullArray6 > (*_local_cell_volume);
    fullArray7 _rQ_hat = reinterpret_cast < fullArray7 > (*_local_Q_hat);
    fullArray8 _rU = reinterpret_cast < fullArray8 > (*_local_U);
    fullArray9 _rV = reinterpret_cast < fullArray9 > (*_local_V);
    fullArray10 _rdPdF = reinterpret_cast < fullArray10 > (*_local_dPdF);
    fullArray11 _rdf = reinterpret_cast < fullArray11 > (*_local_df);
    fullArray12 _rdq = reinterpret_cast < fullArray12 > (*_local_dq);

    const int ChunkSize = 4;
    // chunk typedef
    //typedef int (&refArray1)[3][8];
    typedef float (&refArray1)[3][8][16];
    typedef float (&refArray2)[16];
    typedef float (&refArray3)[16];
    typedef float (&refArray4)[16];
    typedef float (&refArray5)[16];
    typedef float (&refArray6)[16];
    typedef float (&refArray7)[3][16];
    typedef float (&refArray8)[9][16];
    typedef float (&refArray9)[9][16];
    typedef float (&refArray10)[12][16];
    typedef float (&refArray11)[3][8][16];
    typedef float (&refArray12)[16];

    for (int chunk_offset = 0; chunk_offset < 16; chunk_offset += ChunkSize)
      {
        //refArray1 Ak = reinterpret_cast < refArray1 > (_rA[entry][0][chunk_offset]);
        refArray1 duk =
          reinterpret_cast < refArray1 > (_rdu[index][0][0][chunk_offset]);
        refArray2 dpk =
          reinterpret_cast < refArray2 > (_rdp[index][chunk_offset]);
        refArray3 alpha_squared_over_kappak =
          reinterpret_cast < refArray3 >
          (_ralpha_squared_over_kappa[index][chunk_offset]);
        refArray4 alphak =
          reinterpret_cast < refArray4 > (_ralpha[index][chunk_offset]);
        refArray5 one_over_hk =
          reinterpret_cast < refArray5 > (_rone_over_h[index][chunk_offset]);
        refArray6 cell_volumek =
          reinterpret_cast < refArray6 > (_rcell_volume[index][chunk_offset]);
        refArray7 Q_hatk =
          reinterpret_cast < refArray7 > (_rQ_hat[index][0][chunk_offset]);
        refArray8 Uk =
          reinterpret_cast < refArray8 > (_rU[index][0][chunk_offset]);
        refArray9 Vk =
          reinterpret_cast < refArray9 > (_rV[index][0][chunk_offset]);
        refArray10 dPdFk =
          reinterpret_cast < refArray10 > (_rdPdF[index][0][chunk_offset]);
        refArray11 dfk =
          reinterpret_cast < refArray11 > (_rdf[index][0][0][chunk_offset]);
        refArray12 dqk =
          reinterpret_cast < refArray12 > (_rdq[index][chunk_offset]);

        Add_Force_Differential < __m128, float[16], int[16] > (duk, dpk,
                                                               alpha_squared_over_kappak,
                                                               alphak,
                                                               one_over_hk,
                                                               cell_volumek,
                                                               Q_hatk, Uk, Vk,
                                                               dPdFk, dfk, dqk);
      }

  }
};

#endif

#ifdef ENABLE_AVX_INSTRUCTION_SET

template < int SIZE > class Add_Force_Differential_AVX_None
{
private:
  // Generate Variables Here
  float *_local_du;
  float *_local_dp;
  float *_local_alpha_squared_over_kappa;
  float *_local_alpha;
  float *_local_one_over_h;
  float *_local_cell_volume;
  float *_local_Q_hat;
  float *_local_U;
  float *_local_V;
  float *_local_dPdF;
  float *_local_df;
  float *_local_dq;

public:
    explicit Add_Force_Differential_AVX_None (float *du_in, float *dp_in,
                                              float
                                              *alpha_squared_over_kappa_in,
                                              float *alpha_in,
                                              float *one_over_h_in,
                                              float *cell_volume_in,
                                              float *Q_hat_in, float *U_in,
                                              float *V_in, float *dPdF_in,
                                              float *df_in,
                                              float *dq_in):_local_du (du_in),
    _local_dp (dp_in),
    _local_alpha_squared_over_kappa (alpha_squared_over_kappa_in),
    _local_alpha (alpha_in), _local_one_over_h (one_over_h_in),
    _local_cell_volume (cell_volume_in), _local_Q_hat (Q_hat_in),
    _local_U (U_in), _local_V (V_in), _local_dPdF (dPdF_in), _local_df (df_in),
    _local_dq (dq_in)
  {
  }
  void Execute (int index)
  {
    // full array typedefs
    //typedef int (&refArray)[SIZE][3][8];
    typedef float (&fullArray1)[SIZE][3][8][16];
    typedef float (&fullArray2)[SIZE][16];
    typedef float (&fullArray3)[SIZE][16];
    typedef float (&fullArray4)[SIZE][16];
    typedef float (&fullArray5)[SIZE][16];
    typedef float (&fullArray6)[SIZE][16];
    typedef float (&fullArray7)[SIZE][3][16];
    typedef float (&fullArray8)[SIZE][9][16];
    typedef float (&fullArray9)[SIZE][9][16];
    typedef float (&fullArray10)[SIZE][12][16];
    typedef float (&fullArray11)[SIZE][3][8][16];
    typedef float (&fullArray12)[SIZE][16];

    // full array extractions
    //refArray _rA = reinterpret_cast < refArray >(*_A);
    fullArray1 _rdu = reinterpret_cast < fullArray1 > (*_local_du);
    fullArray2 _rdp = reinterpret_cast < fullArray2 > (*_local_dp);
    fullArray3 _ralpha_squared_over_kappa =
      reinterpret_cast < fullArray3 > (*_local_alpha_squared_over_kappa);
    fullArray4 _ralpha = reinterpret_cast < fullArray4 > (*_local_alpha);
    fullArray5 _rone_over_h =
      reinterpret_cast < fullArray5 > (*_local_one_over_h);
    fullArray6 _rcell_volume =
      reinterpret_cast < fullArray6 > (*_local_cell_volume);
    fullArray7 _rQ_hat = reinterpret_cast < fullArray7 > (*_local_Q_hat);
    fullArray8 _rU = reinterpret_cast < fullArray8 > (*_local_U);
    fullArray9 _rV = reinterpret_cast < fullArray9 > (*_local_V);
    fullArray10 _rdPdF = reinterpret_cast < fullArray10 > (*_local_dPdF);
    fullArray11 _rdf = reinterpret_cast < fullArray11 > (*_local_df);
    fullArray12 _rdq = reinterpret_cast < fullArray12 > (*_local_dq);

    const int ChunkSize = 8;
    // chunk typedef
    //typedef int (&refArray1)[3][8];
    typedef float (&refArray1)[3][8][16];
    typedef float (&refArray2)[16];
    typedef float (&refArray3)[16];
    typedef float (&refArray4)[16];
    typedef float (&refArray5)[16];
    typedef float (&refArray6)[16];
    typedef float (&refArray7)[3][16];
    typedef float (&refArray8)[9][16];
    typedef float (&refArray9)[9][16];
    typedef float (&refArray10)[12][16];
    typedef float (&refArray11)[3][8][16];
    typedef float (&refArray12)[16];

    for (int chunk_offset = 0; chunk_offset < 16; chunk_offset += ChunkSize)
      {
        //refArray1 Ak = reinterpret_cast < refArray1 > (_rA[entry][0][chunk_offset]);
        refArray1 duk =
          reinterpret_cast < refArray1 > (_rdu[index][0][0][chunk_offset]);
        refArray2 dpk =
          reinterpret_cast < refArray2 > (_rdp[index][chunk_offset]);
        refArray3 alpha_squared_over_kappak =
          reinterpret_cast < refArray3 >
          (_ralpha_squared_over_kappa[index][chunk_offset]);
        refArray4 alphak =
          reinterpret_cast < refArray4 > (_ralpha[index][chunk_offset]);
        refArray5 one_over_hk =
          reinterpret_cast < refArray5 > (_rone_over_h[index][chunk_offset]);
        refArray6 cell_volumek =
          reinterpret_cast < refArray6 > (_rcell_volume[index][chunk_offset]);
        refArray7 Q_hatk =
          reinterpret_cast < refArray7 > (_rQ_hat[index][0][chunk_offset]);
        refArray8 Uk =
          reinterpret_cast < refArray8 > (_rU[index][0][chunk_offset]);
        refArray9 Vk =
          reinterpret_cast < refArray9 > (_rV[index][0][chunk_offset]);
        refArray10 dPdFk =
          reinterpret_cast < refArray10 > (_rdPdF[index][0][chunk_offset]);
        refArray11 dfk =
          reinterpret_cast < refArray11 > (_rdf[index][0][0][chunk_offset]);
        refArray12 dqk =
          reinterpret_cast < refArray12 > (_rdq[index][chunk_offset]);

        Add_Force_Differential < __m256, float[16], int[16] > (duk, dpk,
                                                               alpha_squared_over_kappak,
                                                               alphak,
                                                               one_over_hk,
                                                               cell_volumek,
                                                               Q_hatk, Uk, Vk,
                                                               dPdFk, dfk, dqk);
      }

  }
};

#endif

#ifdef ENABLE_NEON_INSTRUCTION_SET

template < int SIZE > class Add_Force_Differential_NEON_None
{
private:
  // Generate Variables Here
  float *_local_du;
  float *_local_dp;
  float *_local_alpha_squared_over_kappa;
  float *_local_alpha;
  float *_local_one_over_h;
  float *_local_cell_volume;
  float *_local_Q_hat;
  float *_local_U;
  float *_local_V;
  float *_local_dPdF;
  float *_local_df;
  float *_local_dq;

public:
    explicit Add_Force_Differential_NEON_None (float *du_in, float *dp_in,
                                               float
                                               *alpha_squared_over_kappa_in,
                                               float *alpha_in,
                                               float *one_over_h_in,
                                               float *cell_volume_in,
                                               float *Q_hat_in, float *U_in,
                                               float *V_in, float *dPdF_in,
                                               float *df_in,
                                               float *dq_in):_local_du (du_in),
    _local_dp (dp_in),
    _local_alpha_squared_over_kappa (alpha_squared_over_kappa_in),
    _local_alpha (alpha_in), _local_one_over_h (one_over_h_in),
    _local_cell_volume (cell_volume_in), _local_Q_hat (Q_hat_in),
    _local_U (U_in), _local_V (V_in), _local_dPdF (dPdF_in), _local_df (df_in),
    _local_dq (dq_in)
  {
  }
  void Execute (int index)
  {
    // full array typedefs
    //typedef int (&refArray)[SIZE][3][8];
    typedef float (&fullArray1)[SIZE][3][8][16];
    typedef float (&fullArray2)[SIZE][16];
    typedef float (&fullArray3)[SIZE][16];
    typedef float (&fullArray4)[SIZE][16];
    typedef float (&fullArray5)[SIZE][16];
    typedef float (&fullArray6)[SIZE][16];
    typedef float (&fullArray7)[SIZE][3][16];
    typedef float (&fullArray8)[SIZE][9][16];
    typedef float (&fullArray9)[SIZE][9][16];
    typedef float (&fullArray10)[SIZE][12][16];
    typedef float (&fullArray11)[SIZE][3][8][16];
    typedef float (&fullArray12)[SIZE][16];

    // full array extractions
    //refArray _rA = reinterpret_cast < refArray >(*_A);
    fullArray1 _rdu = reinterpret_cast < fullArray1 > (*_local_du);
    fullArray2 _rdp = reinterpret_cast < fullArray2 > (*_local_dp);
    fullArray3 _ralpha_squared_over_kappa =
      reinterpret_cast < fullArray3 > (*_local_alpha_squared_over_kappa);
    fullArray4 _ralpha = reinterpret_cast < fullArray4 > (*_local_alpha);
    fullArray5 _rone_over_h =
      reinterpret_cast < fullArray5 > (*_local_one_over_h);
    fullArray6 _rcell_volume =
      reinterpret_cast < fullArray6 > (*_local_cell_volume);
    fullArray7 _rQ_hat = reinterpret_cast < fullArray7 > (*_local_Q_hat);
    fullArray8 _rU = reinterpret_cast < fullArray8 > (*_local_U);
    fullArray9 _rV = reinterpret_cast < fullArray9 > (*_local_V);
    fullArray10 _rdPdF = reinterpret_cast < fullArray10 > (*_local_dPdF);
    fullArray11 _rdf = reinterpret_cast < fullArray11 > (*_local_df);
    fullArray12 _rdq = reinterpret_cast < fullArray12 > (*_local_dq);

    const int ChunkSize = 4;
    // chunk typedef
    //typedef int (&refArray1)[3][8];
    typedef float (&refArray1)[3][8][16];
    typedef float (&refArray2)[16];
    typedef float (&refArray3)[16];
    typedef float (&refArray4)[16];
    typedef float (&refArray5)[16];
    typedef float (&refArray6)[16];
    typedef float (&refArray7)[3][16];
    typedef float (&refArray8)[9][16];
    typedef float (&refArray9)[9][16];
    typedef float (&refArray10)[12][16];
    typedef float (&refArray11)[3][8][16];
    typedef float (&refArray12)[16];

    for (int chunk_offset = 0; chunk_offset < 16; chunk_offset += ChunkSize)
      {
        //refArray1 Ak = reinterpret_cast < refArray1 > (_rA[entry][0][chunk_offset]);
        refArray1 duk =
          reinterpret_cast < refArray1 > (_rdu[index][0][0][chunk_offset]);
        refArray2 dpk =
          reinterpret_cast < refArray2 > (_rdp[index][chunk_offset]);
        refArray3 alpha_squared_over_kappak =
          reinterpret_cast < refArray3 >
          (_ralpha_squared_over_kappa[index][chunk_offset]);
        refArray4 alphak =
          reinterpret_cast < refArray4 > (_ralpha[index][chunk_offset]);
        refArray5 one_over_hk =
          reinterpret_cast < refArray5 > (_rone_over_h[index][chunk_offset]);
        refArray6 cell_volumek =
          reinterpret_cast < refArray6 > (_rcell_volume[index][chunk_offset]);
        refArray7 Q_hatk =
          reinterpret_cast < refArray7 > (_rQ_hat[index][0][chunk_offset]);
        refArray8 Uk =
          reinterpret_cast < refArray8 > (_rU[index][0][chunk_offset]);
        refArray9 Vk =
          reinterpret_cast < refArray9 > (_rV[index][0][chunk_offset]);
        refArray10 dPdFk =
          reinterpret_cast < refArray10 > (_rdPdF[index][0][chunk_offset]);
        refArray11 dfk =
          reinterpret_cast < refArray11 > (_rdf[index][0][0][chunk_offset]);
        refArray12 dqk =
          reinterpret_cast < refArray12 > (_rdq[index][chunk_offset]);

        Add_Force_Differential < float32x4_t, float[16], int[16] > (duk, dpk,
                                                                    alpha_squared_over_kappak,
                                                                    alphak,
                                                                    one_over_hk,
                                                                    cell_volumek,
                                                                    Q_hatk, Uk,
                                                                    Vk, dPdFk,
                                                                    dfk, dqk);
      }

  }
};

#endif

#ifdef ENABLE_MIC_INSTRUCTION_SET

template < int SIZE > class Add_Force_Differential_MIC_None
{
private:
  // Generate Variables Here
  float *_local_du;
  float *_local_dp;
  float *_local_alpha_squared_over_kappa;
  float *_local_alpha;
  float *_local_one_over_h;
  float *_local_cell_volume;
  float *_local_Q_hat;
  float *_local_U;
  float *_local_V;
  float *_local_dPdF;
  float *_local_df;
  float *_local_dq;

public:
    explicit Add_Force_Differential_MIC_None (float *du_in, float *dp_in,
                                              float
                                              *alpha_squared_over_kappa_in,
                                              float *alpha_in,
                                              float *one_over_h_in,
                                              float *cell_volume_in,
                                              float *Q_hat_in, float *U_in,
                                              float *V_in, float *dPdF_in,
                                              float *df_in,
                                              float *dq_in):_local_du (du_in),
    _local_dp (dp_in),
    _local_alpha_squared_over_kappa (alpha_squared_over_kappa_in),
    _local_alpha (alpha_in), _local_one_over_h (one_over_h_in),
    _local_cell_volume (cell_volume_in), _local_Q_hat (Q_hat_in),
    _local_U (U_in), _local_V (V_in), _local_dPdF (dPdF_in), _local_df (df_in),
    _local_dq (dq_in)
  {
  }
  void Execute (int index)
  {
    // full array typedefs
    //typedef int (&refArray)[SIZE][3][8];
    typedef float (&fullArray1)[SIZE][3][8][16];
    typedef float (&fullArray2)[SIZE][16];
    typedef float (&fullArray3)[SIZE][16];
    typedef float (&fullArray4)[SIZE][16];
    typedef float (&fullArray5)[SIZE][16];
    typedef float (&fullArray6)[SIZE][16];
    typedef float (&fullArray7)[SIZE][3][16];
    typedef float (&fullArray8)[SIZE][9][16];
    typedef float (&fullArray9)[SIZE][9][16];
    typedef float (&fullArray10)[SIZE][12][16];
    typedef float (&fullArray11)[SIZE][3][8][16];
    typedef float (&fullArray12)[SIZE][16];

    // full array extractions
    //refArray _rA = reinterpret_cast < refArray >(*_A);
    fullArray1 _rdu = reinterpret_cast < fullArray1 > (*_local_du);
    fullArray2 _rdp = reinterpret_cast < fullArray2 > (*_local_dp);
    fullArray3 _ralpha_squared_over_kappa =
      reinterpret_cast < fullArray3 > (*_local_alpha_squared_over_kappa);
    fullArray4 _ralpha = reinterpret_cast < fullArray4 > (*_local_alpha);
    fullArray5 _rone_over_h =
      reinterpret_cast < fullArray5 > (*_local_one_over_h);
    fullArray6 _rcell_volume =
      reinterpret_cast < fullArray6 > (*_local_cell_volume);
    fullArray7 _rQ_hat = reinterpret_cast < fullArray7 > (*_local_Q_hat);
    fullArray8 _rU = reinterpret_cast < fullArray8 > (*_local_U);
    fullArray9 _rV = reinterpret_cast < fullArray9 > (*_local_V);
    fullArray10 _rdPdF = reinterpret_cast < fullArray10 > (*_local_dPdF);
    fullArray11 _rdf = reinterpret_cast < fullArray11 > (*_local_df);
    fullArray12 _rdq = reinterpret_cast < fullArray12 > (*_local_dq);

    const int ChunkSize = 16;
    // chunk typedef
    //typedef int (&refArray1)[3][8];
    typedef float (&refArray1)[3][8][16];
    typedef float (&refArray2)[16];
    typedef float (&refArray3)[16];
    typedef float (&refArray4)[16];
    typedef float (&refArray5)[16];
    typedef float (&refArray6)[16];
    typedef float (&refArray7)[3][16];
    typedef float (&refArray8)[9][16];
    typedef float (&refArray9)[9][16];
    typedef float (&refArray10)[12][16];
    typedef float (&refArray11)[3][8][16];
    typedef float (&refArray12)[16];

    for (int chunk_offset = 0; chunk_offset < 16; chunk_offset += ChunkSize)
      {
        //refArray1 Ak = reinterpret_cast < refArray1 > (_rA[entry][0][chunk_offset]);
        refArray1 duk =
          reinterpret_cast < refArray1 > (_rdu[index][0][0][chunk_offset]);
        refArray2 dpk =
          reinterpret_cast < refArray2 > (_rdp[index][chunk_offset]);
        refArray3 alpha_squared_over_kappak =
          reinterpret_cast < refArray3 >
          (_ralpha_squared_over_kappa[index][chunk_offset]);
        refArray4 alphak =
          reinterpret_cast < refArray4 > (_ralpha[index][chunk_offset]);
        refArray5 one_over_hk =
          reinterpret_cast < refArray5 > (_rone_over_h[index][chunk_offset]);
        refArray6 cell_volumek =
          reinterpret_cast < refArray6 > (_rcell_volume[index][chunk_offset]);
        refArray7 Q_hatk =
          reinterpret_cast < refArray7 > (_rQ_hat[index][0][chunk_offset]);
        refArray8 Uk =
          reinterpret_cast < refArray8 > (_rU[index][0][chunk_offset]);
        refArray9 Vk =
          reinterpret_cast < refArray9 > (_rV[index][0][chunk_offset]);
        refArray10 dPdFk =
          reinterpret_cast < refArray10 > (_rdPdF[index][0][chunk_offset]);
        refArray11 dfk =
          reinterpret_cast < refArray11 > (_rdf[index][0][0][chunk_offset]);
        refArray12 dqk =
          reinterpret_cast < refArray12 > (_rdq[index][chunk_offset]);

        Add_Force_Differential < __m512, float[16], int[16] > (duk, dpk,
                                                               alpha_squared_over_kappak,
                                                               alphak,
                                                               one_over_hk,
                                                               cell_volumek,
                                                               Q_hatk, Uk, Vk,
                                                               dPdFk, dfk, dqk);
      }

  }
};

#endif

int
main (int argc, char *argv[])
{
  typedef float T;

  int seed = 1;
  int threads = 1;
  int threads_max = 1;
  int passes = 1;
  const int data_size = 9610;
  //const int data_size = 1000000;
  if (argc >= 2)
    threads = atoi (argv[1]);
  if (argc >= 3)
    threads_max = atoi (argv[2]);
  if (argc >= 4)
    passes = atoi (argv[3]);
  srand (seed);

  pthread_queue = new PTHREAD_QUEUE (threads_max);

  std::
    cout << "Preparing to Run " << data_size << " of all kernels with " <<
    threads << " threads." << std::endl;



  {
    std::cout << "Running Thread Test for Add_Force_Differential " << std::endl;

//=======================================================
//
//        DEFINE ALL VARIABLES USED BY KERNEL
//
//=======================================================
    std::cout << "\nAllocating all data: ";
    std::cout.flush ();

    start_timer ();
    typedef T (&du_type)[data_size][3][8][16];
    du_type du =
      reinterpret_cast < du_type >
      (*((T *) (_mm_malloc (data_size * 3 * 8 * 16 * sizeof (T), 64))));
    typedef T (&dp_type)[data_size][16];
    dp_type dp =
      reinterpret_cast < dp_type >
      (*((T *) (_mm_malloc (data_size * 16 * sizeof (T), 64))));
    typedef T (&alpha_squared_over_kappa_type)[data_size][16];
    alpha_squared_over_kappa_type alpha_squared_over_kappa =
      reinterpret_cast < alpha_squared_over_kappa_type >
      (*((T *) (_mm_malloc (data_size * 16 * sizeof (T), 64))));
    typedef T (&alpha_type)[data_size][16];
    alpha_type alpha =
      reinterpret_cast < alpha_type >
      (*((T *) (_mm_malloc (data_size * 16 * sizeof (T), 64))));
    typedef T (&one_over_h_type)[data_size][16];
    one_over_h_type one_over_h =
      reinterpret_cast < one_over_h_type >
      (*((T *) (_mm_malloc (data_size * 16 * sizeof (T), 64))));
    typedef T (&cell_volume_type)[data_size][16];
    cell_volume_type cell_volume =
      reinterpret_cast < cell_volume_type >
      (*((T *) (_mm_malloc (data_size * 16 * sizeof (T), 64))));
    typedef T (&Q_hat_type)[data_size][3][16];
    Q_hat_type Q_hat =
      reinterpret_cast < Q_hat_type >
      (*((T *) (_mm_malloc (data_size * 3 * 16 * sizeof (T), 64))));
    typedef T (&U_type)[data_size][9][16];
    U_type U =
      reinterpret_cast < U_type >
      (*((T *) (_mm_malloc (data_size * 9 * 16 * sizeof (T), 64))));
    typedef T (&V_type)[data_size][9][16];
    V_type V =
      reinterpret_cast < V_type >
      (*((T *) (_mm_malloc (data_size * 9 * 16 * sizeof (T), 64))));
    typedef T (&dPdF_type)[data_size][12][16];
    dPdF_type dPdF =
      reinterpret_cast < dPdF_type >
      (*((T *) (_mm_malloc (data_size * 12 * 16 * sizeof (T), 64))));
    typedef T (&df_type)[data_size][3][8][16];
    df_type df =
      reinterpret_cast < df_type >
      (*((T *) (_mm_malloc (data_size * 3 * 8 * 16 * sizeof (T), 64))));
    typedef T (&dq_type)[data_size][16];
    dq_type dq =
      reinterpret_cast < dq_type >
      (*((T *) (_mm_malloc (data_size * 16 * sizeof (T), 64))));


    for (int __a = 0; __a < data_size; __a++)
      for (int __b = 0; __b < 3; __b++)
        for (int __c = 0; __c < 8; __c++)
          for (int __d = 0; __d < 16; __d++)
            {
              du[__a][__b][__c][__d] = Get_Random < float >();
    } for (int __a = 0; __a < data_size; __a++)
      for (int __b = 0; __b < 16; __b++)
        {
          dp[__a][__b] = Get_Random < float >();
    } for (int __a = 0; __a < data_size; __a++)
      for (int __b = 0; __b < 16; __b++)
        {
          alpha_squared_over_kappa[__a][__b] = Get_Random < float >();
    } for (int __a = 0; __a < data_size; __a++)
      for (int __b = 0; __b < 16; __b++)
        {
          alpha[__a][__b] = Get_Random < float >();
    } for (int __a = 0; __a < data_size; __a++)
      for (int __b = 0; __b < 16; __b++)
        {
          one_over_h[__a][__b] = Get_Random < float >();
    } for (int __a = 0; __a < data_size; __a++)
      for (int __b = 0; __b < 16; __b++)
        {
          cell_volume[__a][__b] = Get_Random < float >();
    } for (int __a = 0; __a < data_size; __a++)
      for (int __b = 0; __b < 3; __b++)
        for (int __c = 0; __c < 16; __c++)
          {
            Q_hat[__a][__b][__c] = Get_Random < float >();
    } for (int __a = 0; __a < data_size; __a++)
      for (int __b = 0; __b < 9; __b++)
        for (int __c = 0; __c < 16; __c++)
          {
            U[__a][__b][__c] = Get_Random < float >();
    } for (int __a = 0; __a < data_size; __a++)
      for (int __b = 0; __b < 9; __b++)
        for (int __c = 0; __c < 16; __c++)
          {
            V[__a][__b][__c] = Get_Random < float >();
    } for (int __a = 0; __a < data_size; __a++)
      for (int __b = 0; __b < 12; __b++)
        for (int __c = 0; __c < 16; __c++)
          {
            dPdF[__a][__b][__c] = Get_Random < float >();
    } for (int __a = 0; __a < data_size; __a++)
      for (int __b = 0; __b < 3; __b++)
        for (int __c = 0; __c < 8; __c++)
          for (int __d = 0; __d < 16; __d++)
            {
              df[__a][__b][__c][__d] = Get_Random < float >();
    } for (int __a = 0; __a < data_size; __a++)
      for (int __b = 0; __b < 16; __b++)
        {
          dq[__a][__b] = Get_Random < float >();
        }
    stop_timer ();

    std::cout << get_time () << "s\n\n" << std::endl;


//=======================================================
//
//             COMPUTE SCALAR RESULTS
//
//=======================================================

    {
      std::cout << "	Running " << data_size << " of SCALAR :  " << std::endl;


      Add_Force_Differential_SCALAR_None < data_size > op ((float *) &du,
                                                           (float *) &dp,
                                                           (float *)
                                                           &alpha_squared_over_kappa,
                                                           (float *) &alpha,
                                                           (float *)
                                                           &one_over_h,
                                                           (float *)
                                                           &cell_volume,
                                                           (float *) &Q_hat,
                                                           (float *) &U,
                                                           (float *) &V,
                                                           (float *) &dPdF,
                                                           (float *) &df,
                                                           (float *) &dq);

      for (int t = threads; t <= threads_max;
           t += std::max < int >(((threads_max - threads) / 30), 1))
        {
          std::cout << "Running Test with " << t << " threads." << std::endl;
          MT_Streaming_Kernels::Kernel_Serial_Base_Helper <
            Add_Force_Differential_SCALAR_None < data_size > >helper (op,
                                                                      data_size,
                                                                      t);

          double min_time = 10000000;
          double max_time = -1;
          double avg_time = 0;

          for (int i = 0; i < passes; i++)
            {
              start_timer ();
              helper.Run_Parallel ();
              stop_timer ();
              std::cout << get_time () << "s" << std::endl;
              min_time = std::min < double >(min_time, get_time ());
              max_time = std::max < double >(max_time, get_time ());
              avg_time += get_time ();
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
      std::cout << "	Running " << data_size << " of SSE :  " << std::endl;


      Add_Force_Differential_SSE_None < data_size > op ((float *) &du,
                                                        (float *) &dp,
                                                        (float *)
                                                        &alpha_squared_over_kappa,
                                                        (float *) &alpha,
                                                        (float *) &one_over_h,
                                                        (float *) &cell_volume,
                                                        (float *) &Q_hat,
                                                        (float *) &U,
                                                        (float *) &V,
                                                        (float *) &dPdF,
                                                        (float *) &df,
                                                        (float *) &dq);

      for (int t = threads; t <= threads_max;
           t += std::max < int >(((threads_max - threads) / 30), 1))
        {
          std::cout << "Running Test with " << t << " threads." << std::endl;
          MT_Streaming_Kernels::Kernel_Serial_Base_Helper <
            Add_Force_Differential_SSE_None < data_size > >helper (op,
                                                                   data_size,
                                                                   t);

          double min_time = 10000000;
          double max_time = -1;
          double avg_time = 0;

          for (int i = 0; i < passes; i++)
            {
              start_timer ();
              helper.Run_Parallel ();
              stop_timer ();
              std::cout << get_time () << "s" << std::endl;
              min_time = std::min < double >(min_time, get_time ());
              max_time = std::max < double >(max_time, get_time ());
              avg_time += get_time ();
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
      std::cout << "	Running " << data_size << " of AVX :  " << std::endl;


    typedef T (*vertex_ptr_type)[data_size][3][8][16];
    typedef T (*cell_ptr_type)  [data_size][16];
    typedef T (*DMat3_ptr_type) [data_size][3][16];
    typedef T (*Mat3_ptr_type)  [data_size][9][16];
    typedef T (*dPdF_ptr_type)  [data_size][12][16];

    vertex_ptr_type du_ptr=reinterpret_cast<float (*) [data_size][3][8][16] >(&du[0][0][0][0]);
    cell_ptr_type dp_ptr=reinterpret_cast<float (*) [data_size][16] >(&dp[0][0]);
    cell_ptr_type alpha_squared_over_kappa_ptr=reinterpret_cast<float (*) [data_size][16] >(&alpha_squared_over_kappa[0][0]);
    cell_ptr_type alpha_ptr=reinterpret_cast<float (*) [data_size][16] >(&alpha[0][0]);
    cell_ptr_type one_over_h_ptr=reinterpret_cast<float (*) [data_size][16] >(&one_over_h[0][0]);
    cell_ptr_type cell_volume_ptr=reinterpret_cast<float (*) [data_size][16] >(&cell_volume[0][0]);
    DMat3_ptr_type Q_hat_ptr=reinterpret_cast<float (*) [data_size][3][16] >(&Q_hat[0][0][0]);
    Mat3_ptr_type U_ptr=reinterpret_cast<float (*) [data_size][9][16] >(&U[0][0][0]);
    Mat3_ptr_type V_ptr=reinterpret_cast<float (*) [data_size][9][16] >(&V[0][0][0]);
    dPdF_ptr_type dPdF_ptr=reinterpret_cast<float (*) [data_size][12][16] >(&dPdF[0][0][0]);
    vertex_ptr_type df_ptr=reinterpret_cast<float (*) [data_size][3][8][16] >(&df[0][0][0][0]);
    cell_ptr_type dq_ptr=reinterpret_cast<float (*) [data_size][16] >(&dq[0][0]);


    vertex_ptr_type du2_ptr=reinterpret_cast<float (*) [data_size][3][8][16] >(&du[0][0][0][8]);
    cell_ptr_type dp2_ptr=reinterpret_cast<float (*) [data_size][16] >(&dp[0][8]);
    cell_ptr_type alpha_squared_over_kappa2_ptr=reinterpret_cast<float (*) [data_size][16] >(&alpha_squared_over_kappa[0][8]);
    cell_ptr_type alpha2_ptr=reinterpret_cast<float (*) [data_size][16] >(&alpha[0][8]);
    cell_ptr_type one_over_h2_ptr=reinterpret_cast<float (*) [data_size][16] >(&one_over_h[0][8]);
    cell_ptr_type cell_volume2_ptr=reinterpret_cast<float (*) [data_size][16] >(&cell_volume[0][8]);
    DMat3_ptr_type Q_hat2_ptr=reinterpret_cast<float (*) [data_size][3][16] >(&Q_hat[0][0][8]);
    Mat3_ptr_type U2_ptr=reinterpret_cast<float (*) [data_size][9][16] >(&U[0][0][8]);
    Mat3_ptr_type V2_ptr=reinterpret_cast<float (*) [data_size][9][16] >(&V[0][0][8]);
    dPdF_ptr_type dPdF2_ptr=reinterpret_cast<float (*) [data_size][12][16] >(&dPdF[0][0][8]);
    vertex_ptr_type df2_ptr=reinterpret_cast<float (*) [data_size][3][8][16] >(&df[0][0][0][8]);
    cell_ptr_type dq2_ptr=reinterpret_cast<float (*) [data_size][16] >(&dq[0][8]);

/*
      Add_Force_Differential_AVX_None < data_size > op ((float *) &du,
                                                        (float *) &dp,
                                                        (float *)
                                                        &alpha_squared_over_kappa,
                                                        (float *) &alpha,
                                                        (float *) &one_over_h,
                                                        (float *) &cell_volume,
                                                        (float *) &Q_hat,
                                                        (float *) &U,
                                                        (float *) &V,
                                                        (float *) &dPdF,
                                                        (float *) &df,
                                                        (float *) &dq);

      for (int t = threads; t <= threads_max;
           t += std::max < int >(((threads_max - threads) / 30), 1))
        {
          std::cout << "Running Test with " << t << " threads." << std::endl;
          MT_Streaming_Kernels::Kernel_Serial_Base_Helper <
            Add_Force_Differential_AVX_None < data_size > >helper (op,
                                                                   data_size,
                                                                   t);

          double min_time = 10000000;
          double max_time = -1;
          double avg_time = 0;
*/
          for (int i = 0; i < passes; i++)
            {
              start_timer ();
#if 0
              helper.Run_Parallel_OMP ();
#else
#pragma omp parallel for num_threads(threads)
              for( int j=0; j<data_size; j++){
                  Add_Force_Differential < __m256, float[16], int[16] > (
                      reinterpret_cast< const float (&)[3][8][16] >((*du_ptr)[j]),
                      reinterpret_cast< const float (&)[16] >((*dp_ptr)[j]),
                      reinterpret_cast< const float (&)[16] >((*alpha_squared_over_kappa_ptr)[j]),
                      reinterpret_cast< const float (&)[16] >((*alpha_ptr)[j]),
                      reinterpret_cast< const float (&)[16] >((*one_over_h_ptr)[j]),
                      reinterpret_cast< const float (&)[16] >((*cell_volume_ptr)[j]),
                      reinterpret_cast< const float (&)[3][16] >((*Q_hat_ptr)[j]),
                      reinterpret_cast< const float (&)[9][16] >((*U_ptr)[j]),
                      reinterpret_cast< const float (&)[9][16] >((*V_ptr)[j]),
                      reinterpret_cast< const float (&)[12][16] >((*dPdF_ptr)[j]),
                      reinterpret_cast< float (&)[3][8][16] >((*df_ptr)[j]),
                      reinterpret_cast< float (&)[16] >((*dq_ptr)[j]));   

                  Add_Force_Differential < __m256, float[16], int[16] > (
                      reinterpret_cast< const float (&)[3][8][16] >((*du2_ptr)[j]),
                      reinterpret_cast< const float (&)[16] >((*dp2_ptr)[j]),
                      reinterpret_cast< const float (&)[16] >((*alpha_squared_over_kappa2_ptr)[j]),
                      reinterpret_cast< const float (&)[16] >((*alpha2_ptr)[j]),
                      reinterpret_cast< const float (&)[16] >((*one_over_h2_ptr)[j]),
                      reinterpret_cast< const float (&)[16] >((*cell_volume2_ptr)[j]),
                      reinterpret_cast< const float (&)[3][16] >((*Q_hat2_ptr)[j]),
                      reinterpret_cast< const float (&)[9][16] >((*U2_ptr)[j]),
                      reinterpret_cast< const float (&)[9][16] >((*V2_ptr)[j]),
                      reinterpret_cast< const float (&)[12][16] >((*dPdF2_ptr)[j]),
                      reinterpret_cast< float (&)[3][8][16] >((*df2_ptr)[j]),
                      reinterpret_cast< float (&)[16] >((*dq2_ptr)[j]));   
              }
#endif
              stop_timer ();
              std::cout << get_time () << "s" << std::endl;

            }
//          avg_time = avg_time / passes;
//          std::cout << "Min pass time: " << min_time << std::endl;
//          std::cout << "Max pass time: " << max_time << std::endl;
//          std::cout << "Avg pass time: " << avg_time << std::endl;
    
    }
#endif

//=======================================================
//
//             COMPUTE NEON RESULTS
//
//=======================================================
#ifdef ENABLE_NEON_INSTRUCTION_SET
    {
      std::cout << "	Running " << data_size << " of NEON :  " << std::endl;


      Add_Force_Differential_NEON_None < data_size > op ((float *) &du,
                                                         (float *) &dp,
                                                         (float *)
                                                         &alpha_squared_over_kappa,
                                                         (float *) &alpha,
                                                         (float *) &one_over_h,
                                                         (float *) &cell_volume,
                                                         (float *) &Q_hat,
                                                         (float *) &U,
                                                         (float *) &V,
                                                         (float *) &dPdF,
                                                         (float *) &df,
                                                         (float *) &dq);

      for (int t = threads; t < threads_max;
           t += std::max < int >(((threads_max - threads) / 30), 1))
        {
          std::cout << "Running Test with " << t << " threads." << std::endl;
          MT_Streaming_Kernels::Kernel_Serial_Base_Helper <
            Add_Force_Differential_NEON_None < data_size > >helper (op,
                                                                    data_size,
                                                                    t);

          double min_time = 10000000;
          double max_time = -1;
          double avg_time = 0;

          for (int i = 0; i < passes; i++)
            {
              start_timer ();
              helper.Run_Parallel ();
              stop_timer ();
              std::cout << get_time () << "s" << std::endl;
              min_time = std::min < double >(min_time, get_time ());
              max_time = std::max < double >(max_time, get_time ());
              avg_time += get_time ();
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
      std::cout << "	Running " << data_size << " of NEON :  " << std::endl;


      Add_Force_Differential_NEON_None < data_size > op ((float *) &du,
                                                         (float *) &dp,
                                                         (float *)
                                                         &alpha_squared_over_kappa,
                                                         (float *) &alpha,
                                                         (float *) &one_over_h,
                                                         (float *) &cell_volume,
                                                         (float *) &Q_hat,
                                                         (float *) &U,
                                                         (float *) &V,
                                                         (float *) &dPdF,
                                                         (float *) &df,
                                                         (float *) &dq);

      for (int t = threads; t <= threads_max;
           t += std::max < int >(((threads_max - threads) / 30), 1))
        {
          std::cout << "Running Test with " << t << " threads." << std::endl;
          MT_Streaming_Kernels::Kernel_Serial_Base_Helper <
            Add_Force_Differential_NEON_None < data_size > >helper (op,
                                                                    data_size,
                                                                    t);

          double min_time = 10000000;
          double max_time = -1;
          double avg_time = 0;

          for (int i = 0; i < passes; i++)
            {
              start_timer ();
              helper.Run_Parallel ();
              stop_timer ();
              std::cout << get_time () << "s" << std::endl;
              min_time = std::min < double >(min_time, get_time ());
              max_time = std::max < double >(max_time, get_time ());
              avg_time += get_time ();
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

    typedef T (*vertex_ptr_type)[data_size][3][8][16];
    typedef T (*cell_ptr_type)  [data_size][16];
    typedef T (*DMat3_ptr_type) [data_size][3][16];
    typedef T (*Mat3_ptr_type)  [data_size][9][16];
    typedef T (*dPdF_ptr_type)  [data_size][12][16];

    vertex_ptr_type du_ptr=reinterpret_cast<float (*) [data_size][3][8][16] >(&du[0][0][0][0]);
    cell_ptr_type dp_ptr=reinterpret_cast<float (*) [data_size][16] >(&dp[0][0]);
    cell_ptr_type alpha_squared_over_kappa_ptr=reinterpret_cast<float (*) [data_size][16] >(&alpha_squared_over_kappa[0][0]);
    cell_ptr_type alpha_ptr=reinterpret_cast<float (*) [data_size][16] >(&alpha[0][0]);
    cell_ptr_type one_over_h_ptr=reinterpret_cast<float (*) [data_size][16] >(&one_over_h[0][0]);
    cell_ptr_type cell_volume_ptr=reinterpret_cast<float (*) [data_size][16] >(&cell_volume[0][0]);
    DMat3_ptr_type Q_hat_ptr=reinterpret_cast<float (*) [data_size][3][16] >(&Q_hat[0][0][0]);
    Mat3_ptr_type U_ptr=reinterpret_cast<float (*) [data_size][9][16] >(&U[0][0][0]);
    Mat3_ptr_type V_ptr=reinterpret_cast<float (*) [data_size][9][16] >(&V[0][0][0]);
    dPdF_ptr_type dPdF_ptr=reinterpret_cast<float (*) [data_size][12][16] >(&dPdF[0][0][0]);
    vertex_ptr_type df_ptr=reinterpret_cast<float (*) [data_size][3][8][16] >(&df[0][0][0][0]);
    cell_ptr_type dq_ptr=reinterpret_cast<float (*) [data_size][16] >(&dq[0][0]);


      std::cout << "	Running " << data_size << " of MIC :  " << std::endl;


      // Add_Force_Differential_MIC_None < data_size > op ((float *) &du,
      //                                                   (float *) &dp,
      //                                                   (float *)
      //                                                   &alpha_squared_over_kappa,
      //                                                   (float *) &alpha,
      //                                                   (float *) &one_over_h,
      //                                                   (float *) &cell_volume,
      //                                                   (float *) &Q_hat,
      //                                                   (float *) &U,
      //                                                   (float *) &V,
      //                                                   (float *) &dPdF,
      //                                                   (float *) &df,
      //                                                   (float *) &dq);

      // for (int t = threads; t <= threads_max;
      //      t += std::max < int >(((threads_max - threads) / 30), 1))
      //   {
      //     std::cout << "Running Test with " << t << " threads." << std::endl;
      //     MT_Streaming_Kernels::Kernel_Serial_Base_Helper <
      //       Add_Force_Differential_MIC_None < data_size > >helper (op,
      //                                                              data_size,
      //                                                              t);

      //     double min_time = 10000000;
      //     double max_time = -1;
      //     double avg_time = 0;

          for (int i = 0; i < passes; i++)
            {
              start_timer ();
#pragma omp parallel for
              for( int j=0; j<data_size; j++){
                  Add_Force_Differential < __m512, float[16], int[16] > (
                      reinterpret_cast< const float (&)[3][8][16] >((*du_ptr)[j]),
                      reinterpret_cast< const float (&)[16] >((*dp_ptr)[j]),
                      reinterpret_cast< const float (&)[16] >((*alpha_squared_over_kappa_ptr)[j]),
                      reinterpret_cast< const float (&)[16] >((*alpha_ptr)[j]),
                      reinterpret_cast< const float (&)[16] >((*one_over_h_ptr)[j]),
                      reinterpret_cast< const float (&)[16] >((*cell_volume_ptr)[j]),
                      reinterpret_cast< const float (&)[3][16] >((*Q_hat_ptr)[j]),
                      reinterpret_cast< const float (&)[9][16] >((*U_ptr)[j]),
                      reinterpret_cast< const float (&)[9][16] >((*V_ptr)[j]),
                      reinterpret_cast< const float (&)[12][16] >((*dPdF_ptr)[j]),
                      reinterpret_cast< float (&)[3][8][16] >((*df_ptr)[j]),
                      reinterpret_cast< float (&)[16] >((*dq_ptr)[j]));   
              }
              stop_timer ();
              std::cout << get_time () << "s" << std::endl;
            }
        
    }
#endif

//=======================================================
//
//        FREE MEMORY USED BY ALL VARIABLES
//
//=======================================================
    std::cout << "\nFreeing all data: " << std::endl;
    std::cout.flush ();

    _mm_free (reinterpret_cast < void *>(du));
    _mm_free (reinterpret_cast < void *>(dp));
    _mm_free (reinterpret_cast < void *>(alpha_squared_over_kappa));
    _mm_free (reinterpret_cast < void *>(alpha));
    _mm_free (reinterpret_cast < void *>(one_over_h));
    _mm_free (reinterpret_cast < void *>(cell_volume));
    _mm_free (reinterpret_cast < void *>(Q_hat));
    _mm_free (reinterpret_cast < void *>(U));
    _mm_free (reinterpret_cast < void *>(V));
    _mm_free (reinterpret_cast < void *>(dPdF));
    _mm_free (reinterpret_cast < void *>(df));
    _mm_free (reinterpret_cast < void *>(dq));


  }


  return 0;

}

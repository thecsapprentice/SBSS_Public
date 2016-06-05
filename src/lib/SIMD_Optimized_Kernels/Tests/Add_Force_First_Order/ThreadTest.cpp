
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <sys/time.h>
#include "KernelCommon.h"
struct NEOHOOKEAN_TAG;
struct COROTATED_TAG;

#include "Add_Force_First_Order.h"

#include <Thread_Queueing/PTHREAD_QUEUE.h>
#include <Kernel_Serial_Base_Helper.h>


template < class T > T Get_Random (const T a = (T) - 1., const T b = (T) 1.)
{
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


template < int SIZE > class Add_Force_First_Order_SCALAR_NEOHOOKEAN
{
private:
  // Generate Variables Here
  float *_local_u;
  float *_local_p;
  float *_local_mu;
  float *_local_alpha;
  float *_local_alpha_sqr_over_kappa;
  float *_local_kappa;
  float *_local_one_over_h;
  float *_local_cell_volume;
  float *_local_U;
  float *_local_V;
  float *_local_Sigma;
  float *_local_Q_Hat;
  float *_local_f;
  float *_local_q;

public:
    explicit Add_Force_First_Order_SCALAR_NEOHOOKEAN (float *u_in, float *p_in,
                                                      float *mu_in,
                                                      float *alpha_in,
                                                      float
                                                      *alpha_sqr_over_kappa_in,
                                                      float *kappa_in,
                                                      float *one_over_h_in,
                                                      float *cell_volume_in,
                                                      float *U_in, float *V_in,
                                                      float *Sigma_in,
                                                      float *Q_Hat_in,
                                                      float *f_in,
                                                      float
                                                      *q_in):_local_u (u_in),
    _local_p (p_in), _local_mu (mu_in), _local_alpha (alpha_in),
    _local_alpha_sqr_over_kappa (alpha_sqr_over_kappa_in),
    _local_kappa (kappa_in), _local_one_over_h (one_over_h_in),
    _local_cell_volume (cell_volume_in), _local_U (U_in), _local_V (V_in),
    _local_Sigma (Sigma_in), _local_Q_Hat (Q_Hat_in), _local_f (f_in),
    _local_q (q_in)
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
    typedef float (&fullArray7)[SIZE][16];
    typedef float (&fullArray8)[SIZE][16];
    typedef float (&fullArray9)[SIZE][9][16];
    typedef float (&fullArray10)[SIZE][9][16];
    typedef float (&fullArray11)[SIZE][3][16];
    typedef float (&fullArray12)[SIZE][3][16];
    typedef float (&fullArray13)[SIZE][3][8][16];
    typedef float (&fullArray14)[SIZE][16];

    // full array extractions
    //refArray _rA = reinterpret_cast < refArray >(*_A);
    fullArray1 _ru = reinterpret_cast < fullArray1 > (*_local_u);
    fullArray2 _rp = reinterpret_cast < fullArray2 > (*_local_p);
    fullArray3 _rmu = reinterpret_cast < fullArray3 > (*_local_mu);
    fullArray4 _ralpha = reinterpret_cast < fullArray4 > (*_local_alpha);
    fullArray5 _ralpha_sqr_over_kappa =
      reinterpret_cast < fullArray5 > (*_local_alpha_sqr_over_kappa);
    fullArray6 _rkappa = reinterpret_cast < fullArray6 > (*_local_kappa);
    fullArray7 _rone_over_h =
      reinterpret_cast < fullArray7 > (*_local_one_over_h);
    fullArray8 _rcell_volume =
      reinterpret_cast < fullArray8 > (*_local_cell_volume);
    fullArray9 _rU = reinterpret_cast < fullArray9 > (*_local_U);
    fullArray10 _rV = reinterpret_cast < fullArray10 > (*_local_V);
    fullArray11 _rSigma = reinterpret_cast < fullArray11 > (*_local_Sigma);
    fullArray12 _rQ_Hat = reinterpret_cast < fullArray12 > (*_local_Q_Hat);
    fullArray13 _rf = reinterpret_cast < fullArray13 > (*_local_f);
    fullArray14 _rq = reinterpret_cast < fullArray14 > (*_local_q);

    const int ChunkSize = 1;
    // chunk typedef
    //typedef int (&refArray1)[3][8];
    typedef float (&refArray1)[3][8][16];
    typedef float (&refArray2)[16];
    typedef float (&refArray3)[16];
    typedef float (&refArray4)[16];
    typedef float (&refArray5)[16];
    typedef float (&refArray6)[16];
    typedef float (&refArray7)[16];
    typedef float (&refArray8)[16];
    typedef float (&refArray9)[9][16];
    typedef float (&refArray10)[9][16];
    typedef float (&refArray11)[3][16];
    typedef float (&refArray12)[3][16];
    typedef float (&refArray13)[3][8][16];
    typedef float (&refArray14)[16];

    for (int chunk_offset = 0; chunk_offset < 16; chunk_offset += ChunkSize)
      {
        //refArray1 Ak = reinterpret_cast < refArray1 > (_rA[entry][0][chunk_offset]);
        refArray1 uk =
          reinterpret_cast < refArray1 > (_ru[index][0][0][chunk_offset]);
        refArray2 pk =
          reinterpret_cast < refArray2 > (_rp[index][chunk_offset]);
        refArray3 muk =
          reinterpret_cast < refArray3 > (_rmu[index][chunk_offset]);
        refArray4 alphak =
          reinterpret_cast < refArray4 > (_ralpha[index][chunk_offset]);
        refArray5 alpha_sqr_over_kappak =
          reinterpret_cast < refArray5 >
          (_ralpha_sqr_over_kappa[index][chunk_offset]);
        refArray6 kappak =
          reinterpret_cast < refArray6 > (_rkappa[index][chunk_offset]);
        refArray7 one_over_hk =
          reinterpret_cast < refArray7 > (_rone_over_h[index][chunk_offset]);
        refArray8 cell_volumek =
          reinterpret_cast < refArray8 > (_rcell_volume[index][chunk_offset]);
        refArray9 Uk =
          reinterpret_cast < refArray9 > (_rU[index][0][chunk_offset]);
        refArray10 Vk =
          reinterpret_cast < refArray10 > (_rV[index][0][chunk_offset]);
        refArray11 Sigmak =
          reinterpret_cast < refArray11 > (_rSigma[index][0][chunk_offset]);
        refArray12 Q_Hatk =
          reinterpret_cast < refArray12 > (_rQ_Hat[index][0][chunk_offset]);
        refArray13 fk =
          reinterpret_cast < refArray13 > (_rf[index][0][0][chunk_offset]);
        refArray14 qk =
          reinterpret_cast < refArray14 > (_rq[index][chunk_offset]);

        Add_Force_First_Order < NEOHOOKEAN_TAG, float, float[16],
          int[16] >::Run (uk, pk, muk, alphak, alpha_sqr_over_kappak, kappak,
                          one_over_hk, cell_volumek, Uk, Vk, Sigmak, Q_Hatk, fk,
                          qk);
      }

  }
};

template < int SIZE > class Add_Force_First_Order_SCALAR_COROTATED
{
private:
  // Generate Variables Here
  float *_local_u;
  float *_local_p;
  float *_local_mu;
  float *_local_alpha;
  float *_local_alpha_sqr_over_kappa;
  float *_local_kappa;
  float *_local_one_over_h;
  float *_local_cell_volume;
  float *_local_U;
  float *_local_V;
  float *_local_Sigma;
  float *_local_Q_Hat;
  float *_local_f;
  float *_local_q;

public:
    explicit Add_Force_First_Order_SCALAR_COROTATED (float *u_in, float *p_in,
                                                     float *mu_in,
                                                     float *alpha_in,
                                                     float
                                                     *alpha_sqr_over_kappa_in,
                                                     float *kappa_in,
                                                     float *one_over_h_in,
                                                     float *cell_volume_in,
                                                     float *U_in, float *V_in,
                                                     float *Sigma_in,
                                                     float *Q_Hat_in,
                                                     float *f_in,
                                                     float
                                                     *q_in):_local_u (u_in),
    _local_p (p_in), _local_mu (mu_in), _local_alpha (alpha_in),
    _local_alpha_sqr_over_kappa (alpha_sqr_over_kappa_in),
    _local_kappa (kappa_in), _local_one_over_h (one_over_h_in),
    _local_cell_volume (cell_volume_in), _local_U (U_in), _local_V (V_in),
    _local_Sigma (Sigma_in), _local_Q_Hat (Q_Hat_in), _local_f (f_in),
    _local_q (q_in)
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
    typedef float (&fullArray7)[SIZE][16];
    typedef float (&fullArray8)[SIZE][16];
    typedef float (&fullArray9)[SIZE][9][16];
    typedef float (&fullArray10)[SIZE][9][16];
    typedef float (&fullArray11)[SIZE][3][16];
    typedef float (&fullArray12)[SIZE][3][16];
    typedef float (&fullArray13)[SIZE][3][8][16];
    typedef float (&fullArray14)[SIZE][16];

    // full array extractions
    //refArray _rA = reinterpret_cast < refArray >(*_A);
    fullArray1 _ru = reinterpret_cast < fullArray1 > (*_local_u);
    fullArray2 _rp = reinterpret_cast < fullArray2 > (*_local_p);
    fullArray3 _rmu = reinterpret_cast < fullArray3 > (*_local_mu);
    fullArray4 _ralpha = reinterpret_cast < fullArray4 > (*_local_alpha);
    fullArray5 _ralpha_sqr_over_kappa =
      reinterpret_cast < fullArray5 > (*_local_alpha_sqr_over_kappa);
    fullArray6 _rkappa = reinterpret_cast < fullArray6 > (*_local_kappa);
    fullArray7 _rone_over_h =
      reinterpret_cast < fullArray7 > (*_local_one_over_h);
    fullArray8 _rcell_volume =
      reinterpret_cast < fullArray8 > (*_local_cell_volume);
    fullArray9 _rU = reinterpret_cast < fullArray9 > (*_local_U);
    fullArray10 _rV = reinterpret_cast < fullArray10 > (*_local_V);
    fullArray11 _rSigma = reinterpret_cast < fullArray11 > (*_local_Sigma);
    fullArray12 _rQ_Hat = reinterpret_cast < fullArray12 > (*_local_Q_Hat);
    fullArray13 _rf = reinterpret_cast < fullArray13 > (*_local_f);
    fullArray14 _rq = reinterpret_cast < fullArray14 > (*_local_q);

    const int ChunkSize = 1;
    // chunk typedef
    //typedef int (&refArray1)[3][8];
    typedef float (&refArray1)[3][8][16];
    typedef float (&refArray2)[16];
    typedef float (&refArray3)[16];
    typedef float (&refArray4)[16];
    typedef float (&refArray5)[16];
    typedef float (&refArray6)[16];
    typedef float (&refArray7)[16];
    typedef float (&refArray8)[16];
    typedef float (&refArray9)[9][16];
    typedef float (&refArray10)[9][16];
    typedef float (&refArray11)[3][16];
    typedef float (&refArray12)[3][16];
    typedef float (&refArray13)[3][8][16];
    typedef float (&refArray14)[16];

    for (int chunk_offset = 0; chunk_offset < 16; chunk_offset += ChunkSize)
      {
        //refArray1 Ak = reinterpret_cast < refArray1 > (_rA[entry][0][chunk_offset]);
        refArray1 uk =
          reinterpret_cast < refArray1 > (_ru[index][0][0][chunk_offset]);
        refArray2 pk =
          reinterpret_cast < refArray2 > (_rp[index][chunk_offset]);
        refArray3 muk =
          reinterpret_cast < refArray3 > (_rmu[index][chunk_offset]);
        refArray4 alphak =
          reinterpret_cast < refArray4 > (_ralpha[index][chunk_offset]);
        refArray5 alpha_sqr_over_kappak =
          reinterpret_cast < refArray5 >
          (_ralpha_sqr_over_kappa[index][chunk_offset]);
        refArray6 kappak =
          reinterpret_cast < refArray6 > (_rkappa[index][chunk_offset]);
        refArray7 one_over_hk =
          reinterpret_cast < refArray7 > (_rone_over_h[index][chunk_offset]);
        refArray8 cell_volumek =
          reinterpret_cast < refArray8 > (_rcell_volume[index][chunk_offset]);
        refArray9 Uk =
          reinterpret_cast < refArray9 > (_rU[index][0][chunk_offset]);
        refArray10 Vk =
          reinterpret_cast < refArray10 > (_rV[index][0][chunk_offset]);
        refArray11 Sigmak =
          reinterpret_cast < refArray11 > (_rSigma[index][0][chunk_offset]);
        refArray12 Q_Hatk =
          reinterpret_cast < refArray12 > (_rQ_Hat[index][0][chunk_offset]);
        refArray13 fk =
          reinterpret_cast < refArray13 > (_rf[index][0][0][chunk_offset]);
        refArray14 qk =
          reinterpret_cast < refArray14 > (_rq[index][chunk_offset]);

        Add_Force_First_Order < COROTATED_TAG, float, float[16],
          int[16] >::Run (uk, pk, muk, alphak, alpha_sqr_over_kappak, kappak,
                          one_over_hk, cell_volumek, Uk, Vk, Sigmak, Q_Hatk, fk,
                          qk);
      }

  }
};


#ifdef ENABLE_SSE_INSTRUCTION_SET

template < int SIZE > class Add_Force_First_Order_SSE_NEOHOOKEAN
{
private:
  // Generate Variables Here
  float *_local_u;
  float *_local_p;
  float *_local_mu;
  float *_local_alpha;
  float *_local_alpha_sqr_over_kappa;
  float *_local_kappa;
  float *_local_one_over_h;
  float *_local_cell_volume;
  float *_local_U;
  float *_local_V;
  float *_local_Sigma;
  float *_local_Q_Hat;
  float *_local_f;
  float *_local_q;

public:
    explicit Add_Force_First_Order_SSE_NEOHOOKEAN (float *u_in, float *p_in,
                                                   float *mu_in,
                                                   float *alpha_in,
                                                   float
                                                   *alpha_sqr_over_kappa_in,
                                                   float *kappa_in,
                                                   float *one_over_h_in,
                                                   float *cell_volume_in,
                                                   float *U_in, float *V_in,
                                                   float *Sigma_in,
                                                   float *Q_Hat_in, float *f_in,
                                                   float *q_in):_local_u (u_in),
    _local_p (p_in), _local_mu (mu_in), _local_alpha (alpha_in),
    _local_alpha_sqr_over_kappa (alpha_sqr_over_kappa_in),
    _local_kappa (kappa_in), _local_one_over_h (one_over_h_in),
    _local_cell_volume (cell_volume_in), _local_U (U_in), _local_V (V_in),
    _local_Sigma (Sigma_in), _local_Q_Hat (Q_Hat_in), _local_f (f_in),
    _local_q (q_in)
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
    typedef float (&fullArray7)[SIZE][16];
    typedef float (&fullArray8)[SIZE][16];
    typedef float (&fullArray9)[SIZE][9][16];
    typedef float (&fullArray10)[SIZE][9][16];
    typedef float (&fullArray11)[SIZE][3][16];
    typedef float (&fullArray12)[SIZE][3][16];
    typedef float (&fullArray13)[SIZE][3][8][16];
    typedef float (&fullArray14)[SIZE][16];

    // full array extractions
    //refArray _rA = reinterpret_cast < refArray >(*_A);
    fullArray1 _ru = reinterpret_cast < fullArray1 > (*_local_u);
    fullArray2 _rp = reinterpret_cast < fullArray2 > (*_local_p);
    fullArray3 _rmu = reinterpret_cast < fullArray3 > (*_local_mu);
    fullArray4 _ralpha = reinterpret_cast < fullArray4 > (*_local_alpha);
    fullArray5 _ralpha_sqr_over_kappa =
      reinterpret_cast < fullArray5 > (*_local_alpha_sqr_over_kappa);
    fullArray6 _rkappa = reinterpret_cast < fullArray6 > (*_local_kappa);
    fullArray7 _rone_over_h =
      reinterpret_cast < fullArray7 > (*_local_one_over_h);
    fullArray8 _rcell_volume =
      reinterpret_cast < fullArray8 > (*_local_cell_volume);
    fullArray9 _rU = reinterpret_cast < fullArray9 > (*_local_U);
    fullArray10 _rV = reinterpret_cast < fullArray10 > (*_local_V);
    fullArray11 _rSigma = reinterpret_cast < fullArray11 > (*_local_Sigma);
    fullArray12 _rQ_Hat = reinterpret_cast < fullArray12 > (*_local_Q_Hat);
    fullArray13 _rf = reinterpret_cast < fullArray13 > (*_local_f);
    fullArray14 _rq = reinterpret_cast < fullArray14 > (*_local_q);

    const int ChunkSize = 4;
    // chunk typedef
    //typedef int (&refArray1)[3][8];
    typedef float (&refArray1)[3][8][16];
    typedef float (&refArray2)[16];
    typedef float (&refArray3)[16];
    typedef float (&refArray4)[16];
    typedef float (&refArray5)[16];
    typedef float (&refArray6)[16];
    typedef float (&refArray7)[16];
    typedef float (&refArray8)[16];
    typedef float (&refArray9)[9][16];
    typedef float (&refArray10)[9][16];
    typedef float (&refArray11)[3][16];
    typedef float (&refArray12)[3][16];
    typedef float (&refArray13)[3][8][16];
    typedef float (&refArray14)[16];

    for (int chunk_offset = 0; chunk_offset < 16; chunk_offset += ChunkSize)
      {
        //refArray1 Ak = reinterpret_cast < refArray1 > (_rA[entry][0][chunk_offset]);
        refArray1 uk =
          reinterpret_cast < refArray1 > (_ru[index][0][0][chunk_offset]);
        refArray2 pk =
          reinterpret_cast < refArray2 > (_rp[index][chunk_offset]);
        refArray3 muk =
          reinterpret_cast < refArray3 > (_rmu[index][chunk_offset]);
        refArray4 alphak =
          reinterpret_cast < refArray4 > (_ralpha[index][chunk_offset]);
        refArray5 alpha_sqr_over_kappak =
          reinterpret_cast < refArray5 >
          (_ralpha_sqr_over_kappa[index][chunk_offset]);
        refArray6 kappak =
          reinterpret_cast < refArray6 > (_rkappa[index][chunk_offset]);
        refArray7 one_over_hk =
          reinterpret_cast < refArray7 > (_rone_over_h[index][chunk_offset]);
        refArray8 cell_volumek =
          reinterpret_cast < refArray8 > (_rcell_volume[index][chunk_offset]);
        refArray9 Uk =
          reinterpret_cast < refArray9 > (_rU[index][0][chunk_offset]);
        refArray10 Vk =
          reinterpret_cast < refArray10 > (_rV[index][0][chunk_offset]);
        refArray11 Sigmak =
          reinterpret_cast < refArray11 > (_rSigma[index][0][chunk_offset]);
        refArray12 Q_Hatk =
          reinterpret_cast < refArray12 > (_rQ_Hat[index][0][chunk_offset]);
        refArray13 fk =
          reinterpret_cast < refArray13 > (_rf[index][0][0][chunk_offset]);
        refArray14 qk =
          reinterpret_cast < refArray14 > (_rq[index][chunk_offset]);

        Add_Force_First_Order < NEOHOOKEAN_TAG, __m128, float[16],
          int[16] >::Run (uk, pk, muk, alphak, alpha_sqr_over_kappak, kappak,
                          one_over_hk, cell_volumek, Uk, Vk, Sigmak, Q_Hatk, fk,
                          qk);
      }

  }
};

template < int SIZE > class Add_Force_First_Order_SSE_COROTATED
{
private:
  // Generate Variables Here
  float *_local_u;
  float *_local_p;
  float *_local_mu;
  float *_local_alpha;
  float *_local_alpha_sqr_over_kappa;
  float *_local_kappa;
  float *_local_one_over_h;
  float *_local_cell_volume;
  float *_local_U;
  float *_local_V;
  float *_local_Sigma;
  float *_local_Q_Hat;
  float *_local_f;
  float *_local_q;

public:
    explicit Add_Force_First_Order_SSE_COROTATED (float *u_in, float *p_in,
                                                  float *mu_in, float *alpha_in,
                                                  float
                                                  *alpha_sqr_over_kappa_in,
                                                  float *kappa_in,
                                                  float *one_over_h_in,
                                                  float *cell_volume_in,
                                                  float *U_in, float *V_in,
                                                  float *Sigma_in,
                                                  float *Q_Hat_in, float *f_in,
                                                  float *q_in):_local_u (u_in),
    _local_p (p_in), _local_mu (mu_in), _local_alpha (alpha_in),
    _local_alpha_sqr_over_kappa (alpha_sqr_over_kappa_in),
    _local_kappa (kappa_in), _local_one_over_h (one_over_h_in),
    _local_cell_volume (cell_volume_in), _local_U (U_in), _local_V (V_in),
    _local_Sigma (Sigma_in), _local_Q_Hat (Q_Hat_in), _local_f (f_in),
    _local_q (q_in)
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
    typedef float (&fullArray7)[SIZE][16];
    typedef float (&fullArray8)[SIZE][16];
    typedef float (&fullArray9)[SIZE][9][16];
    typedef float (&fullArray10)[SIZE][9][16];
    typedef float (&fullArray11)[SIZE][3][16];
    typedef float (&fullArray12)[SIZE][3][16];
    typedef float (&fullArray13)[SIZE][3][8][16];
    typedef float (&fullArray14)[SIZE][16];

    // full array extractions
    //refArray _rA = reinterpret_cast < refArray >(*_A);
    fullArray1 _ru = reinterpret_cast < fullArray1 > (*_local_u);
    fullArray2 _rp = reinterpret_cast < fullArray2 > (*_local_p);
    fullArray3 _rmu = reinterpret_cast < fullArray3 > (*_local_mu);
    fullArray4 _ralpha = reinterpret_cast < fullArray4 > (*_local_alpha);
    fullArray5 _ralpha_sqr_over_kappa =
      reinterpret_cast < fullArray5 > (*_local_alpha_sqr_over_kappa);
    fullArray6 _rkappa = reinterpret_cast < fullArray6 > (*_local_kappa);
    fullArray7 _rone_over_h =
      reinterpret_cast < fullArray7 > (*_local_one_over_h);
    fullArray8 _rcell_volume =
      reinterpret_cast < fullArray8 > (*_local_cell_volume);
    fullArray9 _rU = reinterpret_cast < fullArray9 > (*_local_U);
    fullArray10 _rV = reinterpret_cast < fullArray10 > (*_local_V);
    fullArray11 _rSigma = reinterpret_cast < fullArray11 > (*_local_Sigma);
    fullArray12 _rQ_Hat = reinterpret_cast < fullArray12 > (*_local_Q_Hat);
    fullArray13 _rf = reinterpret_cast < fullArray13 > (*_local_f);
    fullArray14 _rq = reinterpret_cast < fullArray14 > (*_local_q);

    const int ChunkSize = 4;
    // chunk typedef
    //typedef int (&refArray1)[3][8];
    typedef float (&refArray1)[3][8][16];
    typedef float (&refArray2)[16];
    typedef float (&refArray3)[16];
    typedef float (&refArray4)[16];
    typedef float (&refArray5)[16];
    typedef float (&refArray6)[16];
    typedef float (&refArray7)[16];
    typedef float (&refArray8)[16];
    typedef float (&refArray9)[9][16];
    typedef float (&refArray10)[9][16];
    typedef float (&refArray11)[3][16];
    typedef float (&refArray12)[3][16];
    typedef float (&refArray13)[3][8][16];
    typedef float (&refArray14)[16];

    for (int chunk_offset = 0; chunk_offset < 16; chunk_offset += ChunkSize)
      {
        //refArray1 Ak = reinterpret_cast < refArray1 > (_rA[entry][0][chunk_offset]);
        refArray1 uk =
          reinterpret_cast < refArray1 > (_ru[index][0][0][chunk_offset]);
        refArray2 pk =
          reinterpret_cast < refArray2 > (_rp[index][chunk_offset]);
        refArray3 muk =
          reinterpret_cast < refArray3 > (_rmu[index][chunk_offset]);
        refArray4 alphak =
          reinterpret_cast < refArray4 > (_ralpha[index][chunk_offset]);
        refArray5 alpha_sqr_over_kappak =
          reinterpret_cast < refArray5 >
          (_ralpha_sqr_over_kappa[index][chunk_offset]);
        refArray6 kappak =
          reinterpret_cast < refArray6 > (_rkappa[index][chunk_offset]);
        refArray7 one_over_hk =
          reinterpret_cast < refArray7 > (_rone_over_h[index][chunk_offset]);
        refArray8 cell_volumek =
          reinterpret_cast < refArray8 > (_rcell_volume[index][chunk_offset]);
        refArray9 Uk =
          reinterpret_cast < refArray9 > (_rU[index][0][chunk_offset]);
        refArray10 Vk =
          reinterpret_cast < refArray10 > (_rV[index][0][chunk_offset]);
        refArray11 Sigmak =
          reinterpret_cast < refArray11 > (_rSigma[index][0][chunk_offset]);
        refArray12 Q_Hatk =
          reinterpret_cast < refArray12 > (_rQ_Hat[index][0][chunk_offset]);
        refArray13 fk =
          reinterpret_cast < refArray13 > (_rf[index][0][0][chunk_offset]);
        refArray14 qk =
          reinterpret_cast < refArray14 > (_rq[index][chunk_offset]);

        Add_Force_First_Order < COROTATED_TAG, __m128, float[16],
          int[16] >::Run (uk, pk, muk, alphak, alpha_sqr_over_kappak, kappak,
                          one_over_hk, cell_volumek, Uk, Vk, Sigmak, Q_Hatk, fk,
                          qk);
      }

  }
};

#endif

#ifdef ENABLE_AVX_INSTRUCTION_SET

template < int SIZE > class Add_Force_First_Order_AVX_NEOHOOKEAN
{
private:
  // Generate Variables Here
  float *_local_u;
  float *_local_p;
  float *_local_mu;
  float *_local_alpha;
  float *_local_alpha_sqr_over_kappa;
  float *_local_kappa;
  float *_local_one_over_h;
  float *_local_cell_volume;
  float *_local_U;
  float *_local_V;
  float *_local_Sigma;
  float *_local_Q_Hat;
  float *_local_f;
  float *_local_q;

public:
    explicit Add_Force_First_Order_AVX_NEOHOOKEAN (float *u_in, float *p_in,
                                                   float *mu_in,
                                                   float *alpha_in,
                                                   float
                                                   *alpha_sqr_over_kappa_in,
                                                   float *kappa_in,
                                                   float *one_over_h_in,
                                                   float *cell_volume_in,
                                                   float *U_in, float *V_in,
                                                   float *Sigma_in,
                                                   float *Q_Hat_in, float *f_in,
                                                   float *q_in):_local_u (u_in),
    _local_p (p_in), _local_mu (mu_in), _local_alpha (alpha_in),
    _local_alpha_sqr_over_kappa (alpha_sqr_over_kappa_in),
    _local_kappa (kappa_in), _local_one_over_h (one_over_h_in),
    _local_cell_volume (cell_volume_in), _local_U (U_in), _local_V (V_in),
    _local_Sigma (Sigma_in), _local_Q_Hat (Q_Hat_in), _local_f (f_in),
    _local_q (q_in)
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
    typedef float (&fullArray7)[SIZE][16];
    typedef float (&fullArray8)[SIZE][16];
    typedef float (&fullArray9)[SIZE][9][16];
    typedef float (&fullArray10)[SIZE][9][16];
    typedef float (&fullArray11)[SIZE][3][16];
    typedef float (&fullArray12)[SIZE][3][16];
    typedef float (&fullArray13)[SIZE][3][8][16];
    typedef float (&fullArray14)[SIZE][16];

    // full array extractions
    //refArray _rA = reinterpret_cast < refArray >(*_A);
    fullArray1 _ru = reinterpret_cast < fullArray1 > (*_local_u);
    fullArray2 _rp = reinterpret_cast < fullArray2 > (*_local_p);
    fullArray3 _rmu = reinterpret_cast < fullArray3 > (*_local_mu);
    fullArray4 _ralpha = reinterpret_cast < fullArray4 > (*_local_alpha);
    fullArray5 _ralpha_sqr_over_kappa =
      reinterpret_cast < fullArray5 > (*_local_alpha_sqr_over_kappa);
    fullArray6 _rkappa = reinterpret_cast < fullArray6 > (*_local_kappa);
    fullArray7 _rone_over_h =
      reinterpret_cast < fullArray7 > (*_local_one_over_h);
    fullArray8 _rcell_volume =
      reinterpret_cast < fullArray8 > (*_local_cell_volume);
    fullArray9 _rU = reinterpret_cast < fullArray9 > (*_local_U);
    fullArray10 _rV = reinterpret_cast < fullArray10 > (*_local_V);
    fullArray11 _rSigma = reinterpret_cast < fullArray11 > (*_local_Sigma);
    fullArray12 _rQ_Hat = reinterpret_cast < fullArray12 > (*_local_Q_Hat);
    fullArray13 _rf = reinterpret_cast < fullArray13 > (*_local_f);
    fullArray14 _rq = reinterpret_cast < fullArray14 > (*_local_q);

    const int ChunkSize = 8;
    // chunk typedef
    //typedef int (&refArray1)[3][8];
    typedef float (&refArray1)[3][8][16];
    typedef float (&refArray2)[16];
    typedef float (&refArray3)[16];
    typedef float (&refArray4)[16];
    typedef float (&refArray5)[16];
    typedef float (&refArray6)[16];
    typedef float (&refArray7)[16];
    typedef float (&refArray8)[16];
    typedef float (&refArray9)[9][16];
    typedef float (&refArray10)[9][16];
    typedef float (&refArray11)[3][16];
    typedef float (&refArray12)[3][16];
    typedef float (&refArray13)[3][8][16];
    typedef float (&refArray14)[16];

    for (int chunk_offset = 0; chunk_offset < 16; chunk_offset += ChunkSize)
      {
        //refArray1 Ak = reinterpret_cast < refArray1 > (_rA[entry][0][chunk_offset]);
        refArray1 uk =
          reinterpret_cast < refArray1 > (_ru[index][0][0][chunk_offset]);
        refArray2 pk =
          reinterpret_cast < refArray2 > (_rp[index][chunk_offset]);
        refArray3 muk =
          reinterpret_cast < refArray3 > (_rmu[index][chunk_offset]);
        refArray4 alphak =
          reinterpret_cast < refArray4 > (_ralpha[index][chunk_offset]);
        refArray5 alpha_sqr_over_kappak =
          reinterpret_cast < refArray5 >
          (_ralpha_sqr_over_kappa[index][chunk_offset]);
        refArray6 kappak =
          reinterpret_cast < refArray6 > (_rkappa[index][chunk_offset]);
        refArray7 one_over_hk =
          reinterpret_cast < refArray7 > (_rone_over_h[index][chunk_offset]);
        refArray8 cell_volumek =
          reinterpret_cast < refArray8 > (_rcell_volume[index][chunk_offset]);
        refArray9 Uk =
          reinterpret_cast < refArray9 > (_rU[index][0][chunk_offset]);
        refArray10 Vk =
          reinterpret_cast < refArray10 > (_rV[index][0][chunk_offset]);
        refArray11 Sigmak =
          reinterpret_cast < refArray11 > (_rSigma[index][0][chunk_offset]);
        refArray12 Q_Hatk =
          reinterpret_cast < refArray12 > (_rQ_Hat[index][0][chunk_offset]);
        refArray13 fk =
          reinterpret_cast < refArray13 > (_rf[index][0][0][chunk_offset]);
        refArray14 qk =
          reinterpret_cast < refArray14 > (_rq[index][chunk_offset]);

        Add_Force_First_Order < NEOHOOKEAN_TAG, __m256, float[16],
          int[16] >::Run (uk, pk, muk, alphak, alpha_sqr_over_kappak, kappak,
                          one_over_hk, cell_volumek, Uk, Vk, Sigmak, Q_Hatk, fk,
                          qk);
      }

  }
};

template < int SIZE > class Add_Force_First_Order_AVX_COROTATED
{
private:
  // Generate Variables Here
  float *_local_u;
  float *_local_p;
  float *_local_mu;
  float *_local_alpha;
  float *_local_alpha_sqr_over_kappa;
  float *_local_kappa;
  float *_local_one_over_h;
  float *_local_cell_volume;
  float *_local_U;
  float *_local_V;
  float *_local_Sigma;
  float *_local_Q_Hat;
  float *_local_f;
  float *_local_q;

public:
    explicit Add_Force_First_Order_AVX_COROTATED (float *u_in, float *p_in,
                                                  float *mu_in, float *alpha_in,
                                                  float
                                                  *alpha_sqr_over_kappa_in,
                                                  float *kappa_in,
                                                  float *one_over_h_in,
                                                  float *cell_volume_in,
                                                  float *U_in, float *V_in,
                                                  float *Sigma_in,
                                                  float *Q_Hat_in, float *f_in,
                                                  float *q_in):_local_u (u_in),
    _local_p (p_in), _local_mu (mu_in), _local_alpha (alpha_in),
    _local_alpha_sqr_over_kappa (alpha_sqr_over_kappa_in),
    _local_kappa (kappa_in), _local_one_over_h (one_over_h_in),
    _local_cell_volume (cell_volume_in), _local_U (U_in), _local_V (V_in),
    _local_Sigma (Sigma_in), _local_Q_Hat (Q_Hat_in), _local_f (f_in),
    _local_q (q_in)
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
    typedef float (&fullArray7)[SIZE][16];
    typedef float (&fullArray8)[SIZE][16];
    typedef float (&fullArray9)[SIZE][9][16];
    typedef float (&fullArray10)[SIZE][9][16];
    typedef float (&fullArray11)[SIZE][3][16];
    typedef float (&fullArray12)[SIZE][3][16];
    typedef float (&fullArray13)[SIZE][3][8][16];
    typedef float (&fullArray14)[SIZE][16];

    // full array extractions
    //refArray _rA = reinterpret_cast < refArray >(*_A);
    fullArray1 _ru = reinterpret_cast < fullArray1 > (*_local_u);
    fullArray2 _rp = reinterpret_cast < fullArray2 > (*_local_p);
    fullArray3 _rmu = reinterpret_cast < fullArray3 > (*_local_mu);
    fullArray4 _ralpha = reinterpret_cast < fullArray4 > (*_local_alpha);
    fullArray5 _ralpha_sqr_over_kappa =
      reinterpret_cast < fullArray5 > (*_local_alpha_sqr_over_kappa);
    fullArray6 _rkappa = reinterpret_cast < fullArray6 > (*_local_kappa);
    fullArray7 _rone_over_h =
      reinterpret_cast < fullArray7 > (*_local_one_over_h);
    fullArray8 _rcell_volume =
      reinterpret_cast < fullArray8 > (*_local_cell_volume);
    fullArray9 _rU = reinterpret_cast < fullArray9 > (*_local_U);
    fullArray10 _rV = reinterpret_cast < fullArray10 > (*_local_V);
    fullArray11 _rSigma = reinterpret_cast < fullArray11 > (*_local_Sigma);
    fullArray12 _rQ_Hat = reinterpret_cast < fullArray12 > (*_local_Q_Hat);
    fullArray13 _rf = reinterpret_cast < fullArray13 > (*_local_f);
    fullArray14 _rq = reinterpret_cast < fullArray14 > (*_local_q);

    const int ChunkSize = 8;
    // chunk typedef
    //typedef int (&refArray1)[3][8];
    typedef float (&refArray1)[3][8][16];
    typedef float (&refArray2)[16];
    typedef float (&refArray3)[16];
    typedef float (&refArray4)[16];
    typedef float (&refArray5)[16];
    typedef float (&refArray6)[16];
    typedef float (&refArray7)[16];
    typedef float (&refArray8)[16];
    typedef float (&refArray9)[9][16];
    typedef float (&refArray10)[9][16];
    typedef float (&refArray11)[3][16];
    typedef float (&refArray12)[3][16];
    typedef float (&refArray13)[3][8][16];
    typedef float (&refArray14)[16];

    for (int chunk_offset = 0; chunk_offset < 16; chunk_offset += ChunkSize)
      {
        //refArray1 Ak = reinterpret_cast < refArray1 > (_rA[entry][0][chunk_offset]);
        refArray1 uk =
          reinterpret_cast < refArray1 > (_ru[index][0][0][chunk_offset]);
        refArray2 pk =
          reinterpret_cast < refArray2 > (_rp[index][chunk_offset]);
        refArray3 muk =
          reinterpret_cast < refArray3 > (_rmu[index][chunk_offset]);
        refArray4 alphak =
          reinterpret_cast < refArray4 > (_ralpha[index][chunk_offset]);
        refArray5 alpha_sqr_over_kappak =
          reinterpret_cast < refArray5 >
          (_ralpha_sqr_over_kappa[index][chunk_offset]);
        refArray6 kappak =
          reinterpret_cast < refArray6 > (_rkappa[index][chunk_offset]);
        refArray7 one_over_hk =
          reinterpret_cast < refArray7 > (_rone_over_h[index][chunk_offset]);
        refArray8 cell_volumek =
          reinterpret_cast < refArray8 > (_rcell_volume[index][chunk_offset]);
        refArray9 Uk =
          reinterpret_cast < refArray9 > (_rU[index][0][chunk_offset]);
        refArray10 Vk =
          reinterpret_cast < refArray10 > (_rV[index][0][chunk_offset]);
        refArray11 Sigmak =
          reinterpret_cast < refArray11 > (_rSigma[index][0][chunk_offset]);
        refArray12 Q_Hatk =
          reinterpret_cast < refArray12 > (_rQ_Hat[index][0][chunk_offset]);
        refArray13 fk =
          reinterpret_cast < refArray13 > (_rf[index][0][0][chunk_offset]);
        refArray14 qk =
          reinterpret_cast < refArray14 > (_rq[index][chunk_offset]);

        Add_Force_First_Order < COROTATED_TAG, __m256, float[16],
          int[16] >::Run (uk, pk, muk, alphak, alpha_sqr_over_kappak, kappak,
                          one_over_hk, cell_volumek, Uk, Vk, Sigmak, Q_Hatk, fk,
                          qk);
      }

  }
};

#endif

#ifdef ENABLE_NEON_INSTRUCTION_SET

template < int SIZE > class Add_Force_First_Order_NEON_NEOHOOKEAN
{
private:
  // Generate Variables Here
  float *_local_u;
  float *_local_p;
  float *_local_mu;
  float *_local_alpha;
  float *_local_alpha_sqr_over_kappa;
  float *_local_kappa;
  float *_local_one_over_h;
  float *_local_cell_volume;
  float *_local_U;
  float *_local_V;
  float *_local_Sigma;
  float *_local_Q_Hat;
  float *_local_f;
  float *_local_q;

public:
    explicit Add_Force_First_Order_NEON_NEOHOOKEAN (float *u_in, float *p_in,
                                                    float *mu_in,
                                                    float *alpha_in,
                                                    float
                                                    *alpha_sqr_over_kappa_in,
                                                    float *kappa_in,
                                                    float *one_over_h_in,
                                                    float *cell_volume_in,
                                                    float *U_in, float *V_in,
                                                    float *Sigma_in,
                                                    float *Q_Hat_in,
                                                    float *f_in,
                                                    float
                                                    *q_in):_local_u (u_in),
    _local_p (p_in), _local_mu (mu_in), _local_alpha (alpha_in),
    _local_alpha_sqr_over_kappa (alpha_sqr_over_kappa_in),
    _local_kappa (kappa_in), _local_one_over_h (one_over_h_in),
    _local_cell_volume (cell_volume_in), _local_U (U_in), _local_V (V_in),
    _local_Sigma (Sigma_in), _local_Q_Hat (Q_Hat_in), _local_f (f_in),
    _local_q (q_in)
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
    typedef float (&fullArray7)[SIZE][16];
    typedef float (&fullArray8)[SIZE][16];
    typedef float (&fullArray9)[SIZE][9][16];
    typedef float (&fullArray10)[SIZE][9][16];
    typedef float (&fullArray11)[SIZE][3][16];
    typedef float (&fullArray12)[SIZE][3][16];
    typedef float (&fullArray13)[SIZE][3][8][16];
    typedef float (&fullArray14)[SIZE][16];

    // full array extractions
    //refArray _rA = reinterpret_cast < refArray >(*_A);
    fullArray1 _ru = reinterpret_cast < fullArray1 > (*_local_u);
    fullArray2 _rp = reinterpret_cast < fullArray2 > (*_local_p);
    fullArray3 _rmu = reinterpret_cast < fullArray3 > (*_local_mu);
    fullArray4 _ralpha = reinterpret_cast < fullArray4 > (*_local_alpha);
    fullArray5 _ralpha_sqr_over_kappa =
      reinterpret_cast < fullArray5 > (*_local_alpha_sqr_over_kappa);
    fullArray6 _rkappa = reinterpret_cast < fullArray6 > (*_local_kappa);
    fullArray7 _rone_over_h =
      reinterpret_cast < fullArray7 > (*_local_one_over_h);
    fullArray8 _rcell_volume =
      reinterpret_cast < fullArray8 > (*_local_cell_volume);
    fullArray9 _rU = reinterpret_cast < fullArray9 > (*_local_U);
    fullArray10 _rV = reinterpret_cast < fullArray10 > (*_local_V);
    fullArray11 _rSigma = reinterpret_cast < fullArray11 > (*_local_Sigma);
    fullArray12 _rQ_Hat = reinterpret_cast < fullArray12 > (*_local_Q_Hat);
    fullArray13 _rf = reinterpret_cast < fullArray13 > (*_local_f);
    fullArray14 _rq = reinterpret_cast < fullArray14 > (*_local_q);

    const int ChunkSize = 4;
    // chunk typedef
    //typedef int (&refArray1)[3][8];
    typedef float (&refArray1)[3][8][16];
    typedef float (&refArray2)[16];
    typedef float (&refArray3)[16];
    typedef float (&refArray4)[16];
    typedef float (&refArray5)[16];
    typedef float (&refArray6)[16];
    typedef float (&refArray7)[16];
    typedef float (&refArray8)[16];
    typedef float (&refArray9)[9][16];
    typedef float (&refArray10)[9][16];
    typedef float (&refArray11)[3][16];
    typedef float (&refArray12)[3][16];
    typedef float (&refArray13)[3][8][16];
    typedef float (&refArray14)[16];

    for (int chunk_offset = 0; chunk_offset < 16; chunk_offset += ChunkSize)
      {
        //refArray1 Ak = reinterpret_cast < refArray1 > (_rA[entry][0][chunk_offset]);
        refArray1 uk =
          reinterpret_cast < refArray1 > (_ru[index][0][0][chunk_offset]);
        refArray2 pk =
          reinterpret_cast < refArray2 > (_rp[index][chunk_offset]);
        refArray3 muk =
          reinterpret_cast < refArray3 > (_rmu[index][chunk_offset]);
        refArray4 alphak =
          reinterpret_cast < refArray4 > (_ralpha[index][chunk_offset]);
        refArray5 alpha_sqr_over_kappak =
          reinterpret_cast < refArray5 >
          (_ralpha_sqr_over_kappa[index][chunk_offset]);
        refArray6 kappak =
          reinterpret_cast < refArray6 > (_rkappa[index][chunk_offset]);
        refArray7 one_over_hk =
          reinterpret_cast < refArray7 > (_rone_over_h[index][chunk_offset]);
        refArray8 cell_volumek =
          reinterpret_cast < refArray8 > (_rcell_volume[index][chunk_offset]);
        refArray9 Uk =
          reinterpret_cast < refArray9 > (_rU[index][0][chunk_offset]);
        refArray10 Vk =
          reinterpret_cast < refArray10 > (_rV[index][0][chunk_offset]);
        refArray11 Sigmak =
          reinterpret_cast < refArray11 > (_rSigma[index][0][chunk_offset]);
        refArray12 Q_Hatk =
          reinterpret_cast < refArray12 > (_rQ_Hat[index][0][chunk_offset]);
        refArray13 fk =
          reinterpret_cast < refArray13 > (_rf[index][0][0][chunk_offset]);
        refArray14 qk =
          reinterpret_cast < refArray14 > (_rq[index][chunk_offset]);

        Add_Force_First_Order < NEOHOOKEAN_TAG, float32x4_t, float[16],
          int[16] >::Run (uk, pk, muk, alphak, alpha_sqr_over_kappak, kappak,
                          one_over_hk, cell_volumek, Uk, Vk, Sigmak, Q_Hatk, fk,
                          qk);
      }

  }
};

template < int SIZE > class Add_Force_First_Order_NEON_COROTATED
{
private:
  // Generate Variables Here
  float *_local_u;
  float *_local_p;
  float *_local_mu;
  float *_local_alpha;
  float *_local_alpha_sqr_over_kappa;
  float *_local_kappa;
  float *_local_one_over_h;
  float *_local_cell_volume;
  float *_local_U;
  float *_local_V;
  float *_local_Sigma;
  float *_local_Q_Hat;
  float *_local_f;
  float *_local_q;

public:
    explicit Add_Force_First_Order_NEON_COROTATED (float *u_in, float *p_in,
                                                   float *mu_in,
                                                   float *alpha_in,
                                                   float
                                                   *alpha_sqr_over_kappa_in,
                                                   float *kappa_in,
                                                   float *one_over_h_in,
                                                   float *cell_volume_in,
                                                   float *U_in, float *V_in,
                                                   float *Sigma_in,
                                                   float *Q_Hat_in, float *f_in,
                                                   float *q_in):_local_u (u_in),
    _local_p (p_in), _local_mu (mu_in), _local_alpha (alpha_in),
    _local_alpha_sqr_over_kappa (alpha_sqr_over_kappa_in),
    _local_kappa (kappa_in), _local_one_over_h (one_over_h_in),
    _local_cell_volume (cell_volume_in), _local_U (U_in), _local_V (V_in),
    _local_Sigma (Sigma_in), _local_Q_Hat (Q_Hat_in), _local_f (f_in),
    _local_q (q_in)
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
    typedef float (&fullArray7)[SIZE][16];
    typedef float (&fullArray8)[SIZE][16];
    typedef float (&fullArray9)[SIZE][9][16];
    typedef float (&fullArray10)[SIZE][9][16];
    typedef float (&fullArray11)[SIZE][3][16];
    typedef float (&fullArray12)[SIZE][3][16];
    typedef float (&fullArray13)[SIZE][3][8][16];
    typedef float (&fullArray14)[SIZE][16];

    // full array extractions
    //refArray _rA = reinterpret_cast < refArray >(*_A);
    fullArray1 _ru = reinterpret_cast < fullArray1 > (*_local_u);
    fullArray2 _rp = reinterpret_cast < fullArray2 > (*_local_p);
    fullArray3 _rmu = reinterpret_cast < fullArray3 > (*_local_mu);
    fullArray4 _ralpha = reinterpret_cast < fullArray4 > (*_local_alpha);
    fullArray5 _ralpha_sqr_over_kappa =
      reinterpret_cast < fullArray5 > (*_local_alpha_sqr_over_kappa);
    fullArray6 _rkappa = reinterpret_cast < fullArray6 > (*_local_kappa);
    fullArray7 _rone_over_h =
      reinterpret_cast < fullArray7 > (*_local_one_over_h);
    fullArray8 _rcell_volume =
      reinterpret_cast < fullArray8 > (*_local_cell_volume);
    fullArray9 _rU = reinterpret_cast < fullArray9 > (*_local_U);
    fullArray10 _rV = reinterpret_cast < fullArray10 > (*_local_V);
    fullArray11 _rSigma = reinterpret_cast < fullArray11 > (*_local_Sigma);
    fullArray12 _rQ_Hat = reinterpret_cast < fullArray12 > (*_local_Q_Hat);
    fullArray13 _rf = reinterpret_cast < fullArray13 > (*_local_f);
    fullArray14 _rq = reinterpret_cast < fullArray14 > (*_local_q);

    const int ChunkSize = 4;
    // chunk typedef
    //typedef int (&refArray1)[3][8];
    typedef float (&refArray1)[3][8][16];
    typedef float (&refArray2)[16];
    typedef float (&refArray3)[16];
    typedef float (&refArray4)[16];
    typedef float (&refArray5)[16];
    typedef float (&refArray6)[16];
    typedef float (&refArray7)[16];
    typedef float (&refArray8)[16];
    typedef float (&refArray9)[9][16];
    typedef float (&refArray10)[9][16];
    typedef float (&refArray11)[3][16];
    typedef float (&refArray12)[3][16];
    typedef float (&refArray13)[3][8][16];
    typedef float (&refArray14)[16];

    for (int chunk_offset = 0; chunk_offset < 16; chunk_offset += ChunkSize)
      {
        //refArray1 Ak = reinterpret_cast < refArray1 > (_rA[entry][0][chunk_offset]);
        refArray1 uk =
          reinterpret_cast < refArray1 > (_ru[index][0][0][chunk_offset]);
        refArray2 pk =
          reinterpret_cast < refArray2 > (_rp[index][chunk_offset]);
        refArray3 muk =
          reinterpret_cast < refArray3 > (_rmu[index][chunk_offset]);
        refArray4 alphak =
          reinterpret_cast < refArray4 > (_ralpha[index][chunk_offset]);
        refArray5 alpha_sqr_over_kappak =
          reinterpret_cast < refArray5 >
          (_ralpha_sqr_over_kappa[index][chunk_offset]);
        refArray6 kappak =
          reinterpret_cast < refArray6 > (_rkappa[index][chunk_offset]);
        refArray7 one_over_hk =
          reinterpret_cast < refArray7 > (_rone_over_h[index][chunk_offset]);
        refArray8 cell_volumek =
          reinterpret_cast < refArray8 > (_rcell_volume[index][chunk_offset]);
        refArray9 Uk =
          reinterpret_cast < refArray9 > (_rU[index][0][chunk_offset]);
        refArray10 Vk =
          reinterpret_cast < refArray10 > (_rV[index][0][chunk_offset]);
        refArray11 Sigmak =
          reinterpret_cast < refArray11 > (_rSigma[index][0][chunk_offset]);
        refArray12 Q_Hatk =
          reinterpret_cast < refArray12 > (_rQ_Hat[index][0][chunk_offset]);
        refArray13 fk =
          reinterpret_cast < refArray13 > (_rf[index][0][0][chunk_offset]);
        refArray14 qk =
          reinterpret_cast < refArray14 > (_rq[index][chunk_offset]);

        Add_Force_First_Order < COROTATED_TAG, float32x4_t, float[16],
          int[16] >::Run (uk, pk, muk, alphak, alpha_sqr_over_kappak, kappak,
                          one_over_hk, cell_volumek, Uk, Vk, Sigmak, Q_Hatk, fk,
                          qk);
      }

  }
};

#endif

#ifdef ENABLE_MIC_INSTRUCTION_SET

template < int SIZE > class Add_Force_First_Order_MIC_NEOHOOKEAN
{
private:
  // Generate Variables Here
  float *_local_u;
  float *_local_p;
  float *_local_mu;
  float *_local_alpha;
  float *_local_alpha_sqr_over_kappa;
  float *_local_kappa;
  float *_local_one_over_h;
  float *_local_cell_volume;
  float *_local_U;
  float *_local_V;
  float *_local_Sigma;
  float *_local_Q_Hat;
  float *_local_f;
  float *_local_q;

public:
    explicit Add_Force_First_Order_MIC_NEOHOOKEAN (float *u_in, float *p_in,
                                                   float *mu_in,
                                                   float *alpha_in,
                                                   float
                                                   *alpha_sqr_over_kappa_in,
                                                   float *kappa_in,
                                                   float *one_over_h_in,
                                                   float *cell_volume_in,
                                                   float *U_in, float *V_in,
                                                   float *Sigma_in,
                                                   float *Q_Hat_in, float *f_in,
                                                   float *q_in):_local_u (u_in),
    _local_p (p_in), _local_mu (mu_in), _local_alpha (alpha_in),
    _local_alpha_sqr_over_kappa (alpha_sqr_over_kappa_in),
    _local_kappa (kappa_in), _local_one_over_h (one_over_h_in),
    _local_cell_volume (cell_volume_in), _local_U (U_in), _local_V (V_in),
    _local_Sigma (Sigma_in), _local_Q_Hat (Q_Hat_in), _local_f (f_in),
    _local_q (q_in)
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
    typedef float (&fullArray7)[SIZE][16];
    typedef float (&fullArray8)[SIZE][16];
    typedef float (&fullArray9)[SIZE][9][16];
    typedef float (&fullArray10)[SIZE][9][16];
    typedef float (&fullArray11)[SIZE][3][16];
    typedef float (&fullArray12)[SIZE][3][16];
    typedef float (&fullArray13)[SIZE][3][8][16];
    typedef float (&fullArray14)[SIZE][16];

    // full array extractions
    //refArray _rA = reinterpret_cast < refArray >(*_A);
    fullArray1 _ru = reinterpret_cast < fullArray1 > (*_local_u);
    fullArray2 _rp = reinterpret_cast < fullArray2 > (*_local_p);
    fullArray3 _rmu = reinterpret_cast < fullArray3 > (*_local_mu);
    fullArray4 _ralpha = reinterpret_cast < fullArray4 > (*_local_alpha);
    fullArray5 _ralpha_sqr_over_kappa =
      reinterpret_cast < fullArray5 > (*_local_alpha_sqr_over_kappa);
    fullArray6 _rkappa = reinterpret_cast < fullArray6 > (*_local_kappa);
    fullArray7 _rone_over_h =
      reinterpret_cast < fullArray7 > (*_local_one_over_h);
    fullArray8 _rcell_volume =
      reinterpret_cast < fullArray8 > (*_local_cell_volume);
    fullArray9 _rU = reinterpret_cast < fullArray9 > (*_local_U);
    fullArray10 _rV = reinterpret_cast < fullArray10 > (*_local_V);
    fullArray11 _rSigma = reinterpret_cast < fullArray11 > (*_local_Sigma);
    fullArray12 _rQ_Hat = reinterpret_cast < fullArray12 > (*_local_Q_Hat);
    fullArray13 _rf = reinterpret_cast < fullArray13 > (*_local_f);
    fullArray14 _rq = reinterpret_cast < fullArray14 > (*_local_q);

    const int ChunkSize = 16;
    // chunk typedef
    //typedef int (&refArray1)[3][8];
    typedef float (&refArray1)[3][8][16];
    typedef float (&refArray2)[16];
    typedef float (&refArray3)[16];
    typedef float (&refArray4)[16];
    typedef float (&refArray5)[16];
    typedef float (&refArray6)[16];
    typedef float (&refArray7)[16];
    typedef float (&refArray8)[16];
    typedef float (&refArray9)[9][16];
    typedef float (&refArray10)[9][16];
    typedef float (&refArray11)[3][16];
    typedef float (&refArray12)[3][16];
    typedef float (&refArray13)[3][8][16];
    typedef float (&refArray14)[16];

    for (int chunk_offset = 0; chunk_offset < 16; chunk_offset += ChunkSize)
      {
        //refArray1 Ak = reinterpret_cast < refArray1 > (_rA[entry][0][chunk_offset]);
        refArray1 uk =
          reinterpret_cast < refArray1 > (_ru[index][0][0][chunk_offset]);
        refArray2 pk =
          reinterpret_cast < refArray2 > (_rp[index][chunk_offset]);
        refArray3 muk =
          reinterpret_cast < refArray3 > (_rmu[index][chunk_offset]);
        refArray4 alphak =
          reinterpret_cast < refArray4 > (_ralpha[index][chunk_offset]);
        refArray5 alpha_sqr_over_kappak =
          reinterpret_cast < refArray5 >
          (_ralpha_sqr_over_kappa[index][chunk_offset]);
        refArray6 kappak =
          reinterpret_cast < refArray6 > (_rkappa[index][chunk_offset]);
        refArray7 one_over_hk =
          reinterpret_cast < refArray7 > (_rone_over_h[index][chunk_offset]);
        refArray8 cell_volumek =
          reinterpret_cast < refArray8 > (_rcell_volume[index][chunk_offset]);
        refArray9 Uk =
          reinterpret_cast < refArray9 > (_rU[index][0][chunk_offset]);
        refArray10 Vk =
          reinterpret_cast < refArray10 > (_rV[index][0][chunk_offset]);
        refArray11 Sigmak =
          reinterpret_cast < refArray11 > (_rSigma[index][0][chunk_offset]);
        refArray12 Q_Hatk =
          reinterpret_cast < refArray12 > (_rQ_Hat[index][0][chunk_offset]);
        refArray13 fk =
          reinterpret_cast < refArray13 > (_rf[index][0][0][chunk_offset]);
        refArray14 qk =
          reinterpret_cast < refArray14 > (_rq[index][chunk_offset]);

        Add_Force_First_Order < NEOHOOKEAN_TAG, __m512, float[16],
          int[16] >::Run (uk, pk, muk, alphak, alpha_sqr_over_kappak, kappak,
                          one_over_hk, cell_volumek, Uk, Vk, Sigmak, Q_Hatk, fk,
                          qk);
      }

  }
};

template < int SIZE > class Add_Force_First_Order_MIC_COROTATED
{
private:
  // Generate Variables Here
  float *_local_u;
  float *_local_p;
  float *_local_mu;
  float *_local_alpha;
  float *_local_alpha_sqr_over_kappa;
  float *_local_kappa;
  float *_local_one_over_h;
  float *_local_cell_volume;
  float *_local_U;
  float *_local_V;
  float *_local_Sigma;
  float *_local_Q_Hat;
  float *_local_f;
  float *_local_q;

public:
    explicit Add_Force_First_Order_MIC_COROTATED (float *u_in, float *p_in,
                                                  float *mu_in, float *alpha_in,
                                                  float
                                                  *alpha_sqr_over_kappa_in,
                                                  float *kappa_in,
                                                  float *one_over_h_in,
                                                  float *cell_volume_in,
                                                  float *U_in, float *V_in,
                                                  float *Sigma_in,
                                                  float *Q_Hat_in, float *f_in,
                                                  float *q_in):_local_u (u_in),
    _local_p (p_in), _local_mu (mu_in), _local_alpha (alpha_in),
    _local_alpha_sqr_over_kappa (alpha_sqr_over_kappa_in),
    _local_kappa (kappa_in), _local_one_over_h (one_over_h_in),
    _local_cell_volume (cell_volume_in), _local_U (U_in), _local_V (V_in),
    _local_Sigma (Sigma_in), _local_Q_Hat (Q_Hat_in), _local_f (f_in),
    _local_q (q_in)
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
    typedef float (&fullArray7)[SIZE][16];
    typedef float (&fullArray8)[SIZE][16];
    typedef float (&fullArray9)[SIZE][9][16];
    typedef float (&fullArray10)[SIZE][9][16];
    typedef float (&fullArray11)[SIZE][3][16];
    typedef float (&fullArray12)[SIZE][3][16];
    typedef float (&fullArray13)[SIZE][3][8][16];
    typedef float (&fullArray14)[SIZE][16];

    // full array extractions
    //refArray _rA = reinterpret_cast < refArray >(*_A);
    fullArray1 _ru = reinterpret_cast < fullArray1 > (*_local_u);
    fullArray2 _rp = reinterpret_cast < fullArray2 > (*_local_p);
    fullArray3 _rmu = reinterpret_cast < fullArray3 > (*_local_mu);
    fullArray4 _ralpha = reinterpret_cast < fullArray4 > (*_local_alpha);
    fullArray5 _ralpha_sqr_over_kappa =
      reinterpret_cast < fullArray5 > (*_local_alpha_sqr_over_kappa);
    fullArray6 _rkappa = reinterpret_cast < fullArray6 > (*_local_kappa);
    fullArray7 _rone_over_h =
      reinterpret_cast < fullArray7 > (*_local_one_over_h);
    fullArray8 _rcell_volume =
      reinterpret_cast < fullArray8 > (*_local_cell_volume);
    fullArray9 _rU = reinterpret_cast < fullArray9 > (*_local_U);
    fullArray10 _rV = reinterpret_cast < fullArray10 > (*_local_V);
    fullArray11 _rSigma = reinterpret_cast < fullArray11 > (*_local_Sigma);
    fullArray12 _rQ_Hat = reinterpret_cast < fullArray12 > (*_local_Q_Hat);
    fullArray13 _rf = reinterpret_cast < fullArray13 > (*_local_f);
    fullArray14 _rq = reinterpret_cast < fullArray14 > (*_local_q);

    const int ChunkSize = 16;
    // chunk typedef
    //typedef int (&refArray1)[3][8];
    typedef float (&refArray1)[3][8][16];
    typedef float (&refArray2)[16];
    typedef float (&refArray3)[16];
    typedef float (&refArray4)[16];
    typedef float (&refArray5)[16];
    typedef float (&refArray6)[16];
    typedef float (&refArray7)[16];
    typedef float (&refArray8)[16];
    typedef float (&refArray9)[9][16];
    typedef float (&refArray10)[9][16];
    typedef float (&refArray11)[3][16];
    typedef float (&refArray12)[3][16];
    typedef float (&refArray13)[3][8][16];
    typedef float (&refArray14)[16];

    for (int chunk_offset = 0; chunk_offset < 16; chunk_offset += ChunkSize)
      {
        //refArray1 Ak = reinterpret_cast < refArray1 > (_rA[entry][0][chunk_offset]);
        refArray1 uk =
          reinterpret_cast < refArray1 > (_ru[index][0][0][chunk_offset]);
        refArray2 pk =
          reinterpret_cast < refArray2 > (_rp[index][chunk_offset]);
        refArray3 muk =
          reinterpret_cast < refArray3 > (_rmu[index][chunk_offset]);
        refArray4 alphak =
          reinterpret_cast < refArray4 > (_ralpha[index][chunk_offset]);
        refArray5 alpha_sqr_over_kappak =
          reinterpret_cast < refArray5 >
          (_ralpha_sqr_over_kappa[index][chunk_offset]);
        refArray6 kappak =
          reinterpret_cast < refArray6 > (_rkappa[index][chunk_offset]);
        refArray7 one_over_hk =
          reinterpret_cast < refArray7 > (_rone_over_h[index][chunk_offset]);
        refArray8 cell_volumek =
          reinterpret_cast < refArray8 > (_rcell_volume[index][chunk_offset]);
        refArray9 Uk =
          reinterpret_cast < refArray9 > (_rU[index][0][chunk_offset]);
        refArray10 Vk =
          reinterpret_cast < refArray10 > (_rV[index][0][chunk_offset]);
        refArray11 Sigmak =
          reinterpret_cast < refArray11 > (_rSigma[index][0][chunk_offset]);
        refArray12 Q_Hatk =
          reinterpret_cast < refArray12 > (_rQ_Hat[index][0][chunk_offset]);
        refArray13 fk =
          reinterpret_cast < refArray13 > (_rf[index][0][0][chunk_offset]);
        refArray14 qk =
          reinterpret_cast < refArray14 > (_rq[index][chunk_offset]);

        Add_Force_First_Order < COROTATED_TAG, __m512, float[16],
          int[16] >::Run (uk, pk, muk, alphak, alpha_sqr_over_kappak, kappak,
                          one_over_hk, cell_volumek, Uk, Vk, Sigmak, Q_Hatk, fk,
                          qk);
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
  const int data_size = 1000000;
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
    std::cout << "Running Thread Test for Add_Force_First_Order " << std::endl;

//=======================================================
//
//        DEFINE ALL VARIABLES USED BY KERNEL
//
//=======================================================
    std::cout << "\nAllocating all data: ";
    std::cout.flush ();

    start_timer ();
    typedef T (&u_type)[data_size][3][8][16];
    u_type u =
      reinterpret_cast < u_type >
      (*((T *) (_mm_malloc (data_size * 3 * 8 * 16 * sizeof (T), 64))));
    typedef T (&p_type)[data_size][16];
    p_type p =
      reinterpret_cast < p_type >
      (*((T *) (_mm_malloc (data_size * 16 * sizeof (T), 64))));
    typedef T (&mu_type)[data_size][16];
    mu_type mu =
      reinterpret_cast < mu_type >
      (*((T *) (_mm_malloc (data_size * 16 * sizeof (T), 64))));
    typedef T (&alpha_type)[data_size][16];
    alpha_type alpha =
      reinterpret_cast < alpha_type >
      (*((T *) (_mm_malloc (data_size * 16 * sizeof (T), 64))));
    typedef T (&alpha_sqr_over_kappa_type)[data_size][16];
    alpha_sqr_over_kappa_type alpha_sqr_over_kappa =
      reinterpret_cast < alpha_sqr_over_kappa_type >
      (*((T *) (_mm_malloc (data_size * 16 * sizeof (T), 64))));
    typedef T (&kappa_type)[data_size][16];
    kappa_type kappa =
      reinterpret_cast < kappa_type >
      (*((T *) (_mm_malloc (data_size * 16 * sizeof (T), 64))));
    typedef T (&one_over_h_type)[data_size][16];
    one_over_h_type one_over_h =
      reinterpret_cast < one_over_h_type >
      (*((T *) (_mm_malloc (data_size * 16 * sizeof (T), 64))));
    typedef T (&cell_volume_type)[data_size][16];
    cell_volume_type cell_volume =
      reinterpret_cast < cell_volume_type >
      (*((T *) (_mm_malloc (data_size * 16 * sizeof (T), 64))));
    typedef T (&U_type)[data_size][9][16];
    U_type U =
      reinterpret_cast < U_type >
      (*((T *) (_mm_malloc (data_size * 9 * 16 * sizeof (T), 64))));
    typedef T (&V_type)[data_size][9][16];
    V_type V =
      reinterpret_cast < V_type >
      (*((T *) (_mm_malloc (data_size * 9 * 16 * sizeof (T), 64))));
    typedef T (&Sigma_type)[data_size][3][16];
    Sigma_type Sigma =
      reinterpret_cast < Sigma_type >
      (*((T *) (_mm_malloc (data_size * 3 * 16 * sizeof (T), 64))));
    typedef T (&Q_Hat_type)[data_size][3][16];
    Q_Hat_type Q_Hat =
      reinterpret_cast < Q_Hat_type >
      (*((T *) (_mm_malloc (data_size * 3 * 16 * sizeof (T), 64))));
    typedef T (&f_type)[data_size][3][8][16];
    f_type f =
      reinterpret_cast < f_type >
      (*((T *) (_mm_malloc (data_size * 3 * 8 * 16 * sizeof (T), 64))));
    typedef T (&q_type)[data_size][16];
    q_type q =
      reinterpret_cast < q_type >
      (*((T *) (_mm_malloc (data_size * 16 * sizeof (T), 64))));


    for (int __a = 0; __a < data_size; __a++)
      for (int __b = 0; __b < 3; __b++)
        for (int __c = 0; __c < 8; __c++)
          for (int __d = 0; __d < 16; __d++)
            {
              u[__a][__b][__c][__d] = Get_Random < float >();
    } for (int __a = 0; __a < data_size; __a++)
      for (int __b = 0; __b < 16; __b++)
        {
          p[__a][__b] = Get_Random < float >();
    } for (int __a = 0; __a < data_size; __a++)
      for (int __b = 0; __b < 16; __b++)
        {
          mu[__a][__b] = Get_Random < float >();
    } for (int __a = 0; __a < data_size; __a++)
      for (int __b = 0; __b < 16; __b++)
        {
          alpha[__a][__b] = Get_Random < float >();
    } for (int __a = 0; __a < data_size; __a++)
      for (int __b = 0; __b < 16; __b++)
        {
          alpha_sqr_over_kappa[__a][__b] = Get_Random < float >();
    } for (int __a = 0; __a < data_size; __a++)
      for (int __b = 0; __b < 16; __b++)
        {
          kappa[__a][__b] = Get_Random < float >();
    } for (int __a = 0; __a < data_size; __a++)
      for (int __b = 0; __b < 16; __b++)
        {
          one_over_h[__a][__b] = Get_Random < float >();
    } for (int __a = 0; __a < data_size; __a++)
      for (int __b = 0; __b < 16; __b++)
        {
          cell_volume[__a][__b] = Get_Random < float >();
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
      for (int __b = 0; __b < 3; __b++)
        for (int __c = 0; __c < 16; __c++)
          {
            Sigma[__a][__b][__c] = Get_Random < float >();
    } for (int __a = 0; __a < data_size; __a++)
      for (int __b = 0; __b < 3; __b++)
        for (int __c = 0; __c < 16; __c++)
          {
            Q_Hat[__a][__b][__c] = Get_Random < float >();
    } for (int __a = 0; __a < data_size; __a++)
      for (int __b = 0; __b < 3; __b++)
        for (int __c = 0; __c < 8; __c++)
          for (int __d = 0; __d < 16; __d++)
            {
              f[__a][__b][__c][__d] = Get_Random < float >();
    } for (int __a = 0; __a < data_size; __a++)
      for (int __b = 0; __b < 16; __b++)
        {
          q[__a][__b] = Get_Random < float >();
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


      Add_Force_First_Order_SCALAR_NEOHOOKEAN < data_size > op ((float *) &u,
                                                                (float *) &p,
                                                                (float *) &mu,
                                                                (float *)
                                                                &alpha,
                                                                (float *)
                                                                &alpha_sqr_over_kappa,
                                                                (float *)
                                                                &kappa,
                                                                (float *)
                                                                &one_over_h,
                                                                (float *)
                                                                &cell_volume,
                                                                (float *) &U,
                                                                (float *) &V,
                                                                (float *)
                                                                &Sigma,
                                                                (float *)
                                                                &Q_Hat,
                                                                (float *) &f,
                                                                (float *) &q);

      for (int t = threads; t <= threads_max;
           t += std::max < int >(((threads_max - threads) / 30), 1))
        {
          std::cout << "Running Test with " << t << " threads." << std::endl;
          MT_Streaming_Kernels::Kernel_Serial_Base_Helper <
            Add_Force_First_Order_SCALAR_NEOHOOKEAN < data_size > >helper (op,
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


      Add_Force_First_Order_SSE_NEOHOOKEAN < data_size > op ((float *) &u,
                                                             (float *) &p,
                                                             (float *) &mu,
                                                             (float *) &alpha,
                                                             (float *)
                                                             &alpha_sqr_over_kappa,
                                                             (float *) &kappa,
                                                             (float *)
                                                             &one_over_h,
                                                             (float *)
                                                             &cell_volume,
                                                             (float *) &U,
                                                             (float *) &V,
                                                             (float *) &Sigma,
                                                             (float *) &Q_Hat,
                                                             (float *) &f,
                                                             (float *) &q);

      for (int t = threads; t <= threads_max;
           t += std::max < int >(((threads_max - threads) / 30), 1))
        {
          std::cout << "Running Test with " << t << " threads." << std::endl;
          MT_Streaming_Kernels::Kernel_Serial_Base_Helper <
            Add_Force_First_Order_SSE_NEOHOOKEAN < data_size > >helper (op,
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


      Add_Force_First_Order_AVX_NEOHOOKEAN < data_size > op ((float *) &u,
                                                             (float *) &p,
                                                             (float *) &mu,
                                                             (float *) &alpha,
                                                             (float *)
                                                             &alpha_sqr_over_kappa,
                                                             (float *) &kappa,
                                                             (float *)
                                                             &one_over_h,
                                                             (float *)
                                                             &cell_volume,
                                                             (float *) &U,
                                                             (float *) &V,
                                                             (float *) &Sigma,
                                                             (float *) &Q_Hat,
                                                             (float *) &f,
                                                             (float *) &q);

      for (int t = threads; t <= threads_max;
           t += std::max < int >(((threads_max - threads) / 30), 1))
        {
          std::cout << "Running Test with " << t << " threads." << std::endl;
          MT_Streaming_Kernels::Kernel_Serial_Base_Helper <
            Add_Force_First_Order_AVX_NEOHOOKEAN < data_size > >helper (op,
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


      Add_Force_First_Order_NEON_NEOHOOKEAN < data_size > op ((float *) &u,
                                                              (float *) &p,
                                                              (float *) &mu,
                                                              (float *) &alpha,
                                                              (float *)
                                                              &alpha_sqr_over_kappa,
                                                              (float *) &kappa,
                                                              (float *)
                                                              &one_over_h,
                                                              (float *)
                                                              &cell_volume,
                                                              (float *) &U,
                                                              (float *) &V,
                                                              (float *) &Sigma,
                                                              (float *) &Q_Hat,
                                                              (float *) &f,
                                                              (float *) &q);

      for (int t = threads; t <= threads_max;
           t += std::max < int >(((threads_max - threads) / 30), 1))
        {
          std::cout << "Running Test with " << t << " threads." << std::endl;
          MT_Streaming_Kernels::Kernel_Serial_Base_Helper <
            Add_Force_First_Order_NEON_NEOHOOKEAN < data_size > >helper (op,
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
      std::cout << "	Running " << data_size << " of MIC :  " << std::endl;


      Add_Force_First_Order_MIC_NEOHOOKEAN < data_size > op ((float *) &u,
                                                             (float *) &p,
                                                             (float *) &mu,
                                                             (float *) &alpha,
                                                             (float *)
                                                             &alpha_sqr_over_kappa,
                                                             (float *) &kappa,
                                                             (float *)
                                                             &one_over_h,
                                                             (float *)
                                                             &cell_volume,
                                                             (float *) &U,
                                                             (float *) &V,
                                                             (float *) &Sigma,
                                                             (float *) &Q_Hat,
                                                             (float *) &f,
                                                             (float *) &q);

      for (int t = threads; t <= threads_max;
           t += std::max < int >(((threads_max - threads) / 30), 1))
        {
          std::cout << "Running Test with " << t << " threads." << std::endl;
          MT_Streaming_Kernels::Kernel_Serial_Base_Helper <
            Add_Force_First_Order_MIC_NEOHOOKEAN < data_size > >helper (op,
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
//        FREE MEMORY USED BY ALL VARIABLES
//
//=======================================================
    std::cout << "\nFreeing all data: " << std::endl;
    std::cout.flush ();

    _mm_free (reinterpret_cast < void *>(u));
    _mm_free (reinterpret_cast < void *>(p));
    _mm_free (reinterpret_cast < void *>(mu));
    _mm_free (reinterpret_cast < void *>(alpha));
    _mm_free (reinterpret_cast < void *>(alpha_sqr_over_kappa));
    _mm_free (reinterpret_cast < void *>(kappa));
    _mm_free (reinterpret_cast < void *>(one_over_h));
    _mm_free (reinterpret_cast < void *>(cell_volume));
    _mm_free (reinterpret_cast < void *>(U));
    _mm_free (reinterpret_cast < void *>(V));
    _mm_free (reinterpret_cast < void *>(Sigma));
    _mm_free (reinterpret_cast < void *>(Q_Hat));
    _mm_free (reinterpret_cast < void *>(f));
    _mm_free (reinterpret_cast < void *>(q));


  }


  {
    std::cout << "Running Thread Test for Add_Force_First_Order " << std::endl;

//=======================================================
//
//        DEFINE ALL VARIABLES USED BY KERNEL
//
//=======================================================
    std::cout << "\nAllocating all data: ";
    std::cout.flush ();

    start_timer ();
    typedef T (&u_type)[data_size][3][8][16];
    u_type u =
      reinterpret_cast < u_type >
      (*((T *) (_mm_malloc (data_size * 3 * 8 * 16 * sizeof (T), 64))));
    typedef T (&p_type)[data_size][16];
    p_type p =
      reinterpret_cast < p_type >
      (*((T *) (_mm_malloc (data_size * 16 * sizeof (T), 64))));
    typedef T (&mu_type)[data_size][16];
    mu_type mu =
      reinterpret_cast < mu_type >
      (*((T *) (_mm_malloc (data_size * 16 * sizeof (T), 64))));
    typedef T (&alpha_type)[data_size][16];
    alpha_type alpha =
      reinterpret_cast < alpha_type >
      (*((T *) (_mm_malloc (data_size * 16 * sizeof (T), 64))));
    typedef T (&alpha_sqr_over_kappa_type)[data_size][16];
    alpha_sqr_over_kappa_type alpha_sqr_over_kappa =
      reinterpret_cast < alpha_sqr_over_kappa_type >
      (*((T *) (_mm_malloc (data_size * 16 * sizeof (T), 64))));
    typedef T (&kappa_type)[data_size][16];
    kappa_type kappa =
      reinterpret_cast < kappa_type >
      (*((T *) (_mm_malloc (data_size * 16 * sizeof (T), 64))));
    typedef T (&one_over_h_type)[data_size][16];
    one_over_h_type one_over_h =
      reinterpret_cast < one_over_h_type >
      (*((T *) (_mm_malloc (data_size * 16 * sizeof (T), 64))));
    typedef T (&cell_volume_type)[data_size][16];
    cell_volume_type cell_volume =
      reinterpret_cast < cell_volume_type >
      (*((T *) (_mm_malloc (data_size * 16 * sizeof (T), 64))));
    typedef T (&U_type)[data_size][9][16];
    U_type U =
      reinterpret_cast < U_type >
      (*((T *) (_mm_malloc (data_size * 9 * 16 * sizeof (T), 64))));
    typedef T (&V_type)[data_size][9][16];
    V_type V =
      reinterpret_cast < V_type >
      (*((T *) (_mm_malloc (data_size * 9 * 16 * sizeof (T), 64))));
    typedef T (&Sigma_type)[data_size][3][16];
    Sigma_type Sigma =
      reinterpret_cast < Sigma_type >
      (*((T *) (_mm_malloc (data_size * 3 * 16 * sizeof (T), 64))));
    typedef T (&Q_Hat_type)[data_size][3][16];
    Q_Hat_type Q_Hat =
      reinterpret_cast < Q_Hat_type >
      (*((T *) (_mm_malloc (data_size * 3 * 16 * sizeof (T), 64))));
    typedef T (&f_type)[data_size][3][8][16];
    f_type f =
      reinterpret_cast < f_type >
      (*((T *) (_mm_malloc (data_size * 3 * 8 * 16 * sizeof (T), 64))));
    typedef T (&q_type)[data_size][16];
    q_type q =
      reinterpret_cast < q_type >
      (*((T *) (_mm_malloc (data_size * 16 * sizeof (T), 64))));


    for (int __a = 0; __a < data_size; __a++)
      for (int __b = 0; __b < 3; __b++)
        for (int __c = 0; __c < 8; __c++)
          for (int __d = 0; __d < 16; __d++)
            {
              u[__a][__b][__c][__d] = Get_Random < float >();
    } for (int __a = 0; __a < data_size; __a++)
      for (int __b = 0; __b < 16; __b++)
        {
          p[__a][__b] = Get_Random < float >();
    } for (int __a = 0; __a < data_size; __a++)
      for (int __b = 0; __b < 16; __b++)
        {
          mu[__a][__b] = Get_Random < float >();
    } for (int __a = 0; __a < data_size; __a++)
      for (int __b = 0; __b < 16; __b++)
        {
          alpha[__a][__b] = Get_Random < float >();
    } for (int __a = 0; __a < data_size; __a++)
      for (int __b = 0; __b < 16; __b++)
        {
          alpha_sqr_over_kappa[__a][__b] = Get_Random < float >();
    } for (int __a = 0; __a < data_size; __a++)
      for (int __b = 0; __b < 16; __b++)
        {
          kappa[__a][__b] = Get_Random < float >();
    } for (int __a = 0; __a < data_size; __a++)
      for (int __b = 0; __b < 16; __b++)
        {
          one_over_h[__a][__b] = Get_Random < float >();
    } for (int __a = 0; __a < data_size; __a++)
      for (int __b = 0; __b < 16; __b++)
        {
          cell_volume[__a][__b] = Get_Random < float >();
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
      for (int __b = 0; __b < 3; __b++)
        for (int __c = 0; __c < 16; __c++)
          {
            Sigma[__a][__b][__c] = Get_Random < float >();
    } for (int __a = 0; __a < data_size; __a++)
      for (int __b = 0; __b < 3; __b++)
        for (int __c = 0; __c < 16; __c++)
          {
            Q_Hat[__a][__b][__c] = Get_Random < float >();
    } for (int __a = 0; __a < data_size; __a++)
      for (int __b = 0; __b < 3; __b++)
        for (int __c = 0; __c < 8; __c++)
          for (int __d = 0; __d < 16; __d++)
            {
              f[__a][__b][__c][__d] = Get_Random < float >();
    } for (int __a = 0; __a < data_size; __a++)
      for (int __b = 0; __b < 16; __b++)
        {
          q[__a][__b] = Get_Random < float >();
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


      Add_Force_First_Order_SCALAR_COROTATED < data_size > op ((float *) &u,
                                                               (float *) &p,
                                                               (float *) &mu,
                                                               (float *) &alpha,
                                                               (float *)
                                                               &alpha_sqr_over_kappa,
                                                               (float *) &kappa,
                                                               (float *)
                                                               &one_over_h,
                                                               (float *)
                                                               &cell_volume,
                                                               (float *) &U,
                                                               (float *) &V,
                                                               (float *) &Sigma,
                                                               (float *) &Q_Hat,
                                                               (float *) &f,
                                                               (float *) &q);

      for (int t = threads; t <= threads_max;
           t += std::max < int >(((threads_max - threads) / 30), 1))
        {
          std::cout << "Running Test with " << t << " threads." << std::endl;
          MT_Streaming_Kernels::Kernel_Serial_Base_Helper <
            Add_Force_First_Order_SCALAR_COROTATED < data_size > >helper (op,
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


      Add_Force_First_Order_SSE_COROTATED < data_size > op ((float *) &u,
                                                            (float *) &p,
                                                            (float *) &mu,
                                                            (float *) &alpha,
                                                            (float *)
                                                            &alpha_sqr_over_kappa,
                                                            (float *) &kappa,
                                                            (float *)
                                                            &one_over_h,
                                                            (float *)
                                                            &cell_volume,
                                                            (float *) &U,
                                                            (float *) &V,
                                                            (float *) &Sigma,
                                                            (float *) &Q_Hat,
                                                            (float *) &f,
                                                            (float *) &q);

      for (int t = threads; t <= threads_max;
           t += std::max < int >(((threads_max - threads) / 30), 1))
        {
          std::cout << "Running Test with " << t << " threads." << std::endl;
          MT_Streaming_Kernels::Kernel_Serial_Base_Helper <
            Add_Force_First_Order_SSE_COROTATED < data_size > >helper (op,
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


      Add_Force_First_Order_AVX_COROTATED < data_size > op ((float *) &u,
                                                            (float *) &p,
                                                            (float *) &mu,
                                                            (float *) &alpha,
                                                            (float *)
                                                            &alpha_sqr_over_kappa,
                                                            (float *) &kappa,
                                                            (float *)
                                                            &one_over_h,
                                                            (float *)
                                                            &cell_volume,
                                                            (float *) &U,
                                                            (float *) &V,
                                                            (float *) &Sigma,
                                                            (float *) &Q_Hat,
                                                            (float *) &f,
                                                            (float *) &q);

      for (int t = threads; t <= threads_max;
           t += std::max < int >(((threads_max - threads) / 30), 1))
        {
          std::cout << "Running Test with " << t << " threads." << std::endl;
          MT_Streaming_Kernels::Kernel_Serial_Base_Helper <
            Add_Force_First_Order_AVX_COROTATED < data_size > >helper (op,
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


      Add_Force_First_Order_NEON_COROTATED < data_size > op ((float *) &u,
                                                             (float *) &p,
                                                             (float *) &mu,
                                                             (float *) &alpha,
                                                             (float *)
                                                             &alpha_sqr_over_kappa,
                                                             (float *) &kappa,
                                                             (float *)
                                                             &one_over_h,
                                                             (float *)
                                                             &cell_volume,
                                                             (float *) &U,
                                                             (float *) &V,
                                                             (float *) &Sigma,
                                                             (float *) &Q_Hat,
                                                             (float *) &f,
                                                             (float *) &q);

      for (int t = threads; t <= threads_max;
           t += std::max < int >(((threads_max - threads) / 30), 1))
        {
          std::cout << "Running Test with " << t << " threads." << std::endl;
          MT_Streaming_Kernels::Kernel_Serial_Base_Helper <
            Add_Force_First_Order_NEON_COROTATED < data_size > >helper (op,
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
      std::cout << "	Running " << data_size << " of MIC :  " << std::endl;


      Add_Force_First_Order_MIC_COROTATED < data_size > op ((float *) &u,
                                                            (float *) &p,
                                                            (float *) &mu,
                                                            (float *) &alpha,
                                                            (float *)
                                                            &alpha_sqr_over_kappa,
                                                            (float *) &kappa,
                                                            (float *)
                                                            &one_over_h,
                                                            (float *)
                                                            &cell_volume,
                                                            (float *) &U,
                                                            (float *) &V,
                                                            (float *) &Sigma,
                                                            (float *) &Q_Hat,
                                                            (float *) &f,
                                                            (float *) &q);

      for (int t = threads; t <= threads_max;
           t += std::max < int >(((threads_max - threads) / 30), 1))
        {
          std::cout << "Running Test with " << t << " threads." << std::endl;
          MT_Streaming_Kernels::Kernel_Serial_Base_Helper <
            Add_Force_First_Order_MIC_COROTATED < data_size > >helper (op,
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
//        FREE MEMORY USED BY ALL VARIABLES
//
//=======================================================
    std::cout << "\nFreeing all data: " << std::endl;
    std::cout.flush ();

    _mm_free (reinterpret_cast < void *>(u));
    _mm_free (reinterpret_cast < void *>(p));
    _mm_free (reinterpret_cast < void *>(mu));
    _mm_free (reinterpret_cast < void *>(alpha));
    _mm_free (reinterpret_cast < void *>(alpha_sqr_over_kappa));
    _mm_free (reinterpret_cast < void *>(kappa));
    _mm_free (reinterpret_cast < void *>(one_over_h));
    _mm_free (reinterpret_cast < void *>(cell_volume));
    _mm_free (reinterpret_cast < void *>(U));
    _mm_free (reinterpret_cast < void *>(V));
    _mm_free (reinterpret_cast < void *>(Sigma));
    _mm_free (reinterpret_cast < void *>(Q_Hat));
    _mm_free (reinterpret_cast < void *>(f));
    _mm_free (reinterpret_cast < void *>(q));


  }


  return 0;

}

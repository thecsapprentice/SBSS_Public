
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <sys/time.h>
#include "KernelCommon.h"
struct NEOHOOKEAN_TAG;
struct COROTATED_TAG;

#include "Augmented_Rotated_Stress_Derivative.h"

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


template < int SIZE >
  class Augmented_Rotated_Stress_Derivative_SCALAR_NEOHOOKEAN
{
private:
  // Generate Variables Here
  float *_local_dPdF;
  float *_local_Sigma;
  float *_local_p;
  float *_local_mu;
  float *_local_alpha;

public:
    explicit Augmented_Rotated_Stress_Derivative_SCALAR_NEOHOOKEAN (float
                                                                    *dPdF_in,
                                                                    float
                                                                    *Sigma_in,
                                                                    float *p_in,
                                                                    float
                                                                    *mu_in,
                                                                    float
                                                                    *alpha_in):_local_dPdF
    (dPdF_in), _local_Sigma (Sigma_in), _local_p (p_in), _local_mu (mu_in),
    _local_alpha (alpha_in)
  {
  }
  void Execute (int index)
  {
    // full array typedefs
    //typedef int (&refArray)[SIZE][3][8];
    typedef float (&fullArray1)[SIZE][12][16];
    typedef float (&fullArray2)[SIZE][3][16];
    typedef float (&fullArray3)[SIZE][16];
    typedef float (&fullArray4)[SIZE][16];
    typedef float (&fullArray5)[SIZE][16];

    // full array extractions
    //refArray _rA = reinterpret_cast < refArray >(*_A);
    fullArray1 _rdPdF = reinterpret_cast < fullArray1 > (*_local_dPdF);
    fullArray2 _rSigma = reinterpret_cast < fullArray2 > (*_local_Sigma);
    fullArray3 _rp = reinterpret_cast < fullArray3 > (*_local_p);
    fullArray4 _rmu = reinterpret_cast < fullArray4 > (*_local_mu);
    fullArray5 _ralpha = reinterpret_cast < fullArray5 > (*_local_alpha);

    const int ChunkSize = 1;
    // chunk typedef
    //typedef int (&refArray1)[3][8];
    typedef float (&refArray1)[12][16];
    typedef float (&refArray2)[3][16];
    typedef float (&refArray3)[16];
    typedef float (&refArray4)[16];
    typedef float (&refArray5)[16];

    for (int chunk_offset = 0; chunk_offset < 16; chunk_offset += ChunkSize)
      {
        //refArray1 Ak = reinterpret_cast < refArray1 > (_rA[entry][0][chunk_offset]);
        refArray1 dPdFk =
          reinterpret_cast < refArray1 > (_rdPdF[index][0][chunk_offset]);
        refArray2 Sigmak =
          reinterpret_cast < refArray2 > (_rSigma[index][0][chunk_offset]);
        refArray3 pk =
          reinterpret_cast < refArray3 > (_rp[index][chunk_offset]);
        refArray4 muk =
          reinterpret_cast < refArray4 > (_rmu[index][chunk_offset]);
        refArray5 alphak =
          reinterpret_cast < refArray5 > (_ralpha[index][chunk_offset]);

        Augmented_Rotated_Stress_Derivative < NEOHOOKEAN_TAG, float, float[16],
          int[16] >::Run (dPdFk, Sigmak, pk, muk, alphak);
      }

  }
};

template < int SIZE > class Augmented_Rotated_Stress_Derivative_SCALAR_COROTATED
{
private:
  // Generate Variables Here
  float *_local_dPdF;
  float *_local_Sigma;
  float *_local_p;
  float *_local_mu;
  float *_local_alpha;

public:
    explicit Augmented_Rotated_Stress_Derivative_SCALAR_COROTATED (float
                                                                   *dPdF_in,
                                                                   float
                                                                   *Sigma_in,
                                                                   float *p_in,
                                                                   float *mu_in,
                                                                   float
                                                                   *alpha_in):_local_dPdF
    (dPdF_in), _local_Sigma (Sigma_in), _local_p (p_in), _local_mu (mu_in),
    _local_alpha (alpha_in)
  {
  }
  void Execute (int index)
  {
    // full array typedefs
    //typedef int (&refArray)[SIZE][3][8];
    typedef float (&fullArray1)[SIZE][12][16];
    typedef float (&fullArray2)[SIZE][3][16];
    typedef float (&fullArray3)[SIZE][16];
    typedef float (&fullArray4)[SIZE][16];
    typedef float (&fullArray5)[SIZE][16];

    // full array extractions
    //refArray _rA = reinterpret_cast < refArray >(*_A);
    fullArray1 _rdPdF = reinterpret_cast < fullArray1 > (*_local_dPdF);
    fullArray2 _rSigma = reinterpret_cast < fullArray2 > (*_local_Sigma);
    fullArray3 _rp = reinterpret_cast < fullArray3 > (*_local_p);
    fullArray4 _rmu = reinterpret_cast < fullArray4 > (*_local_mu);
    fullArray5 _ralpha = reinterpret_cast < fullArray5 > (*_local_alpha);

    const int ChunkSize = 1;
    // chunk typedef
    //typedef int (&refArray1)[3][8];
    typedef float (&refArray1)[12][16];
    typedef float (&refArray2)[3][16];
    typedef float (&refArray3)[16];
    typedef float (&refArray4)[16];
    typedef float (&refArray5)[16];

    for (int chunk_offset = 0; chunk_offset < 16; chunk_offset += ChunkSize)
      {
        //refArray1 Ak = reinterpret_cast < refArray1 > (_rA[entry][0][chunk_offset]);
        refArray1 dPdFk =
          reinterpret_cast < refArray1 > (_rdPdF[index][0][chunk_offset]);
        refArray2 Sigmak =
          reinterpret_cast < refArray2 > (_rSigma[index][0][chunk_offset]);
        refArray3 pk =
          reinterpret_cast < refArray3 > (_rp[index][chunk_offset]);
        refArray4 muk =
          reinterpret_cast < refArray4 > (_rmu[index][chunk_offset]);
        refArray5 alphak =
          reinterpret_cast < refArray5 > (_ralpha[index][chunk_offset]);

        Augmented_Rotated_Stress_Derivative < COROTATED_TAG, float, float[16],
          int[16] >::Run (dPdFk, Sigmak, pk, muk, alphak);
      }

  }
};


#ifdef ENABLE_SSE_INSTRUCTION_SET

template < int SIZE > class Augmented_Rotated_Stress_Derivative_SSE_NEOHOOKEAN
{
private:
  // Generate Variables Here
  float *_local_dPdF;
  float *_local_Sigma;
  float *_local_p;
  float *_local_mu;
  float *_local_alpha;

public:
    explicit Augmented_Rotated_Stress_Derivative_SSE_NEOHOOKEAN (float *dPdF_in,
                                                                 float
                                                                 *Sigma_in,
                                                                 float *p_in,
                                                                 float *mu_in,
                                                                 float
                                                                 *alpha_in):_local_dPdF
    (dPdF_in), _local_Sigma (Sigma_in), _local_p (p_in), _local_mu (mu_in),
    _local_alpha (alpha_in)
  {
  }
  void Execute (int index)
  {
    // full array typedefs
    //typedef int (&refArray)[SIZE][3][8];
    typedef float (&fullArray1)[SIZE][12][16];
    typedef float (&fullArray2)[SIZE][3][16];
    typedef float (&fullArray3)[SIZE][16];
    typedef float (&fullArray4)[SIZE][16];
    typedef float (&fullArray5)[SIZE][16];

    // full array extractions
    //refArray _rA = reinterpret_cast < refArray >(*_A);
    fullArray1 _rdPdF = reinterpret_cast < fullArray1 > (*_local_dPdF);
    fullArray2 _rSigma = reinterpret_cast < fullArray2 > (*_local_Sigma);
    fullArray3 _rp = reinterpret_cast < fullArray3 > (*_local_p);
    fullArray4 _rmu = reinterpret_cast < fullArray4 > (*_local_mu);
    fullArray5 _ralpha = reinterpret_cast < fullArray5 > (*_local_alpha);

    const int ChunkSize = 4;
    // chunk typedef
    //typedef int (&refArray1)[3][8];
    typedef float (&refArray1)[12][16];
    typedef float (&refArray2)[3][16];
    typedef float (&refArray3)[16];
    typedef float (&refArray4)[16];
    typedef float (&refArray5)[16];

    for (int chunk_offset = 0; chunk_offset < 16; chunk_offset += ChunkSize)
      {
        //refArray1 Ak = reinterpret_cast < refArray1 > (_rA[entry][0][chunk_offset]);
        refArray1 dPdFk =
          reinterpret_cast < refArray1 > (_rdPdF[index][0][chunk_offset]);
        refArray2 Sigmak =
          reinterpret_cast < refArray2 > (_rSigma[index][0][chunk_offset]);
        refArray3 pk =
          reinterpret_cast < refArray3 > (_rp[index][chunk_offset]);
        refArray4 muk =
          reinterpret_cast < refArray4 > (_rmu[index][chunk_offset]);
        refArray5 alphak =
          reinterpret_cast < refArray5 > (_ralpha[index][chunk_offset]);

        Augmented_Rotated_Stress_Derivative < NEOHOOKEAN_TAG, __m128, float[16],
          int[16] >::Run (dPdFk, Sigmak, pk, muk, alphak);
      }

  }
};

template < int SIZE > class Augmented_Rotated_Stress_Derivative_SSE_COROTATED
{
private:
  // Generate Variables Here
  float *_local_dPdF;
  float *_local_Sigma;
  float *_local_p;
  float *_local_mu;
  float *_local_alpha;

public:
    explicit Augmented_Rotated_Stress_Derivative_SSE_COROTATED (float *dPdF_in,
                                                                float *Sigma_in,
                                                                float *p_in,
                                                                float *mu_in,
                                                                float
                                                                *alpha_in):_local_dPdF
    (dPdF_in), _local_Sigma (Sigma_in), _local_p (p_in), _local_mu (mu_in),
    _local_alpha (alpha_in)
  {
  }
  void Execute (int index)
  {
    // full array typedefs
    //typedef int (&refArray)[SIZE][3][8];
    typedef float (&fullArray1)[SIZE][12][16];
    typedef float (&fullArray2)[SIZE][3][16];
    typedef float (&fullArray3)[SIZE][16];
    typedef float (&fullArray4)[SIZE][16];
    typedef float (&fullArray5)[SIZE][16];

    // full array extractions
    //refArray _rA = reinterpret_cast < refArray >(*_A);
    fullArray1 _rdPdF = reinterpret_cast < fullArray1 > (*_local_dPdF);
    fullArray2 _rSigma = reinterpret_cast < fullArray2 > (*_local_Sigma);
    fullArray3 _rp = reinterpret_cast < fullArray3 > (*_local_p);
    fullArray4 _rmu = reinterpret_cast < fullArray4 > (*_local_mu);
    fullArray5 _ralpha = reinterpret_cast < fullArray5 > (*_local_alpha);

    const int ChunkSize = 4;
    // chunk typedef
    //typedef int (&refArray1)[3][8];
    typedef float (&refArray1)[12][16];
    typedef float (&refArray2)[3][16];
    typedef float (&refArray3)[16];
    typedef float (&refArray4)[16];
    typedef float (&refArray5)[16];

    for (int chunk_offset = 0; chunk_offset < 16; chunk_offset += ChunkSize)
      {
        //refArray1 Ak = reinterpret_cast < refArray1 > (_rA[entry][0][chunk_offset]);
        refArray1 dPdFk =
          reinterpret_cast < refArray1 > (_rdPdF[index][0][chunk_offset]);
        refArray2 Sigmak =
          reinterpret_cast < refArray2 > (_rSigma[index][0][chunk_offset]);
        refArray3 pk =
          reinterpret_cast < refArray3 > (_rp[index][chunk_offset]);
        refArray4 muk =
          reinterpret_cast < refArray4 > (_rmu[index][chunk_offset]);
        refArray5 alphak =
          reinterpret_cast < refArray5 > (_ralpha[index][chunk_offset]);

        Augmented_Rotated_Stress_Derivative < COROTATED_TAG, __m128, float[16],
          int[16] >::Run (dPdFk, Sigmak, pk, muk, alphak);
      }

  }
};

#endif

#ifdef ENABLE_AVX_INSTRUCTION_SET

template < int SIZE > class Augmented_Rotated_Stress_Derivative_AVX_NEOHOOKEAN
{
private:
  // Generate Variables Here
  float *_local_dPdF;
  float *_local_Sigma;
  float *_local_p;
  float *_local_mu;
  float *_local_alpha;

public:
    explicit Augmented_Rotated_Stress_Derivative_AVX_NEOHOOKEAN (float *dPdF_in,
                                                                 float
                                                                 *Sigma_in,
                                                                 float *p_in,
                                                                 float *mu_in,
                                                                 float
                                                                 *alpha_in):_local_dPdF
    (dPdF_in), _local_Sigma (Sigma_in), _local_p (p_in), _local_mu (mu_in),
    _local_alpha (alpha_in)
  {
  }
  void Execute (int index)
  {
    // full array typedefs
    //typedef int (&refArray)[SIZE][3][8];
    typedef float (&fullArray1)[SIZE][12][16];
    typedef float (&fullArray2)[SIZE][3][16];
    typedef float (&fullArray3)[SIZE][16];
    typedef float (&fullArray4)[SIZE][16];
    typedef float (&fullArray5)[SIZE][16];

    // full array extractions
    //refArray _rA = reinterpret_cast < refArray >(*_A);
    fullArray1 _rdPdF = reinterpret_cast < fullArray1 > (*_local_dPdF);
    fullArray2 _rSigma = reinterpret_cast < fullArray2 > (*_local_Sigma);
    fullArray3 _rp = reinterpret_cast < fullArray3 > (*_local_p);
    fullArray4 _rmu = reinterpret_cast < fullArray4 > (*_local_mu);
    fullArray5 _ralpha = reinterpret_cast < fullArray5 > (*_local_alpha);

    const int ChunkSize = 8;
    // chunk typedef
    //typedef int (&refArray1)[3][8];
    typedef float (&refArray1)[12][16];
    typedef float (&refArray2)[3][16];
    typedef float (&refArray3)[16];
    typedef float (&refArray4)[16];
    typedef float (&refArray5)[16];

    for (int chunk_offset = 0; chunk_offset < 16; chunk_offset += ChunkSize)
      {
        //refArray1 Ak = reinterpret_cast < refArray1 > (_rA[entry][0][chunk_offset]);
        refArray1 dPdFk =
          reinterpret_cast < refArray1 > (_rdPdF[index][0][chunk_offset]);
        refArray2 Sigmak =
          reinterpret_cast < refArray2 > (_rSigma[index][0][chunk_offset]);
        refArray3 pk =
          reinterpret_cast < refArray3 > (_rp[index][chunk_offset]);
        refArray4 muk =
          reinterpret_cast < refArray4 > (_rmu[index][chunk_offset]);
        refArray5 alphak =
          reinterpret_cast < refArray5 > (_ralpha[index][chunk_offset]);

        Augmented_Rotated_Stress_Derivative < NEOHOOKEAN_TAG, __m256, float[16],
          int[16] >::Run (dPdFk, Sigmak, pk, muk, alphak);
      }

  }
};

template < int SIZE > class Augmented_Rotated_Stress_Derivative_AVX_COROTATED
{
private:
  // Generate Variables Here
  float *_local_dPdF;
  float *_local_Sigma;
  float *_local_p;
  float *_local_mu;
  float *_local_alpha;

public:
    explicit Augmented_Rotated_Stress_Derivative_AVX_COROTATED (float *dPdF_in,
                                                                float *Sigma_in,
                                                                float *p_in,
                                                                float *mu_in,
                                                                float
                                                                *alpha_in):_local_dPdF
    (dPdF_in), _local_Sigma (Sigma_in), _local_p (p_in), _local_mu (mu_in),
    _local_alpha (alpha_in)
  {
  }
  void Execute (int index)
  {
    // full array typedefs
    //typedef int (&refArray)[SIZE][3][8];
    typedef float (&fullArray1)[SIZE][12][16];
    typedef float (&fullArray2)[SIZE][3][16];
    typedef float (&fullArray3)[SIZE][16];
    typedef float (&fullArray4)[SIZE][16];
    typedef float (&fullArray5)[SIZE][16];

    // full array extractions
    //refArray _rA = reinterpret_cast < refArray >(*_A);
    fullArray1 _rdPdF = reinterpret_cast < fullArray1 > (*_local_dPdF);
    fullArray2 _rSigma = reinterpret_cast < fullArray2 > (*_local_Sigma);
    fullArray3 _rp = reinterpret_cast < fullArray3 > (*_local_p);
    fullArray4 _rmu = reinterpret_cast < fullArray4 > (*_local_mu);
    fullArray5 _ralpha = reinterpret_cast < fullArray5 > (*_local_alpha);

    const int ChunkSize = 8;
    // chunk typedef
    //typedef int (&refArray1)[3][8];
    typedef float (&refArray1)[12][16];
    typedef float (&refArray2)[3][16];
    typedef float (&refArray3)[16];
    typedef float (&refArray4)[16];
    typedef float (&refArray5)[16];

    for (int chunk_offset = 0; chunk_offset < 16; chunk_offset += ChunkSize)
      {
        //refArray1 Ak = reinterpret_cast < refArray1 > (_rA[entry][0][chunk_offset]);
        refArray1 dPdFk =
          reinterpret_cast < refArray1 > (_rdPdF[index][0][chunk_offset]);
        refArray2 Sigmak =
          reinterpret_cast < refArray2 > (_rSigma[index][0][chunk_offset]);
        refArray3 pk =
          reinterpret_cast < refArray3 > (_rp[index][chunk_offset]);
        refArray4 muk =
          reinterpret_cast < refArray4 > (_rmu[index][chunk_offset]);
        refArray5 alphak =
          reinterpret_cast < refArray5 > (_ralpha[index][chunk_offset]);

        Augmented_Rotated_Stress_Derivative < COROTATED_TAG, __m256, float[16],
          int[16] >::Run (dPdFk, Sigmak, pk, muk, alphak);
      }

  }
};

#endif

#ifdef ENABLE_NEON_INSTRUCTION_SET

template < int SIZE > class Augmented_Rotated_Stress_Derivative_NEON_NEOHOOKEAN
{
private:
  // Generate Variables Here
  float *_local_dPdF;
  float *_local_Sigma;
  float *_local_p;
  float *_local_mu;
  float *_local_alpha;

public:
    explicit Augmented_Rotated_Stress_Derivative_NEON_NEOHOOKEAN (float
                                                                  *dPdF_in,
                                                                  float
                                                                  *Sigma_in,
                                                                  float *p_in,
                                                                  float *mu_in,
                                                                  float
                                                                  *alpha_in):_local_dPdF
    (dPdF_in), _local_Sigma (Sigma_in), _local_p (p_in), _local_mu (mu_in),
    _local_alpha (alpha_in)
  {
  }
  void Execute (int index)
  {
    // full array typedefs
    //typedef int (&refArray)[SIZE][3][8];
    typedef float (&fullArray1)[SIZE][12][16];
    typedef float (&fullArray2)[SIZE][3][16];
    typedef float (&fullArray3)[SIZE][16];
    typedef float (&fullArray4)[SIZE][16];
    typedef float (&fullArray5)[SIZE][16];

    // full array extractions
    //refArray _rA = reinterpret_cast < refArray >(*_A);
    fullArray1 _rdPdF = reinterpret_cast < fullArray1 > (*_local_dPdF);
    fullArray2 _rSigma = reinterpret_cast < fullArray2 > (*_local_Sigma);
    fullArray3 _rp = reinterpret_cast < fullArray3 > (*_local_p);
    fullArray4 _rmu = reinterpret_cast < fullArray4 > (*_local_mu);
    fullArray5 _ralpha = reinterpret_cast < fullArray5 > (*_local_alpha);

    const int ChunkSize = 4;
    // chunk typedef
    //typedef int (&refArray1)[3][8];
    typedef float (&refArray1)[12][16];
    typedef float (&refArray2)[3][16];
    typedef float (&refArray3)[16];
    typedef float (&refArray4)[16];
    typedef float (&refArray5)[16];

    for (int chunk_offset = 0; chunk_offset < 16; chunk_offset += ChunkSize)
      {
        //refArray1 Ak = reinterpret_cast < refArray1 > (_rA[entry][0][chunk_offset]);
        refArray1 dPdFk =
          reinterpret_cast < refArray1 > (_rdPdF[index][0][chunk_offset]);
        refArray2 Sigmak =
          reinterpret_cast < refArray2 > (_rSigma[index][0][chunk_offset]);
        refArray3 pk =
          reinterpret_cast < refArray3 > (_rp[index][chunk_offset]);
        refArray4 muk =
          reinterpret_cast < refArray4 > (_rmu[index][chunk_offset]);
        refArray5 alphak =
          reinterpret_cast < refArray5 > (_ralpha[index][chunk_offset]);

        Augmented_Rotated_Stress_Derivative < NEOHOOKEAN_TAG, float32x4_t,
          float[16], int[16] >::Run (dPdFk, Sigmak, pk, muk, alphak);
      }

  }
};

template < int SIZE > class Augmented_Rotated_Stress_Derivative_NEON_COROTATED
{
private:
  // Generate Variables Here
  float *_local_dPdF;
  float *_local_Sigma;
  float *_local_p;
  float *_local_mu;
  float *_local_alpha;

public:
    explicit Augmented_Rotated_Stress_Derivative_NEON_COROTATED (float *dPdF_in,
                                                                 float
                                                                 *Sigma_in,
                                                                 float *p_in,
                                                                 float *mu_in,
                                                                 float
                                                                 *alpha_in):_local_dPdF
    (dPdF_in), _local_Sigma (Sigma_in), _local_p (p_in), _local_mu (mu_in),
    _local_alpha (alpha_in)
  {
  }
  void Execute (int index)
  {
    // full array typedefs
    //typedef int (&refArray)[SIZE][3][8];
    typedef float (&fullArray1)[SIZE][12][16];
    typedef float (&fullArray2)[SIZE][3][16];
    typedef float (&fullArray3)[SIZE][16];
    typedef float (&fullArray4)[SIZE][16];
    typedef float (&fullArray5)[SIZE][16];

    // full array extractions
    //refArray _rA = reinterpret_cast < refArray >(*_A);
    fullArray1 _rdPdF = reinterpret_cast < fullArray1 > (*_local_dPdF);
    fullArray2 _rSigma = reinterpret_cast < fullArray2 > (*_local_Sigma);
    fullArray3 _rp = reinterpret_cast < fullArray3 > (*_local_p);
    fullArray4 _rmu = reinterpret_cast < fullArray4 > (*_local_mu);
    fullArray5 _ralpha = reinterpret_cast < fullArray5 > (*_local_alpha);

    const int ChunkSize = 4;
    // chunk typedef
    //typedef int (&refArray1)[3][8];
    typedef float (&refArray1)[12][16];
    typedef float (&refArray2)[3][16];
    typedef float (&refArray3)[16];
    typedef float (&refArray4)[16];
    typedef float (&refArray5)[16];

    for (int chunk_offset = 0; chunk_offset < 16; chunk_offset += ChunkSize)
      {
        //refArray1 Ak = reinterpret_cast < refArray1 > (_rA[entry][0][chunk_offset]);
        refArray1 dPdFk =
          reinterpret_cast < refArray1 > (_rdPdF[index][0][chunk_offset]);
        refArray2 Sigmak =
          reinterpret_cast < refArray2 > (_rSigma[index][0][chunk_offset]);
        refArray3 pk =
          reinterpret_cast < refArray3 > (_rp[index][chunk_offset]);
        refArray4 muk =
          reinterpret_cast < refArray4 > (_rmu[index][chunk_offset]);
        refArray5 alphak =
          reinterpret_cast < refArray5 > (_ralpha[index][chunk_offset]);

        Augmented_Rotated_Stress_Derivative < COROTATED_TAG, float32x4_t,
          float[16], int[16] >::Run (dPdFk, Sigmak, pk, muk, alphak);
      }

  }
};

#endif

#ifdef ENABLE_MIC_INSTRUCTION_SET

template < int SIZE > class Augmented_Rotated_Stress_Derivative_MIC_NEOHOOKEAN
{
private:
  // Generate Variables Here
  float *_local_dPdF;
  float *_local_Sigma;
  float *_local_p;
  float *_local_mu;
  float *_local_alpha;

public:
    explicit Augmented_Rotated_Stress_Derivative_MIC_NEOHOOKEAN (float *dPdF_in,
                                                                 float
                                                                 *Sigma_in,
                                                                 float *p_in,
                                                                 float *mu_in,
                                                                 float
                                                                 *alpha_in):_local_dPdF
    (dPdF_in), _local_Sigma (Sigma_in), _local_p (p_in), _local_mu (mu_in),
    _local_alpha (alpha_in)
  {
  }
  void Execute (int index)
  {
    // full array typedefs
    //typedef int (&refArray)[SIZE][3][8];
    typedef float (&fullArray1)[SIZE][12][16];
    typedef float (&fullArray2)[SIZE][3][16];
    typedef float (&fullArray3)[SIZE][16];
    typedef float (&fullArray4)[SIZE][16];
    typedef float (&fullArray5)[SIZE][16];

    // full array extractions
    //refArray _rA = reinterpret_cast < refArray >(*_A);
    fullArray1 _rdPdF = reinterpret_cast < fullArray1 > (*_local_dPdF);
    fullArray2 _rSigma = reinterpret_cast < fullArray2 > (*_local_Sigma);
    fullArray3 _rp = reinterpret_cast < fullArray3 > (*_local_p);
    fullArray4 _rmu = reinterpret_cast < fullArray4 > (*_local_mu);
    fullArray5 _ralpha = reinterpret_cast < fullArray5 > (*_local_alpha);

    const int ChunkSize = 16;
    // chunk typedef
    //typedef int (&refArray1)[3][8];
    typedef float (&refArray1)[12][16];
    typedef float (&refArray2)[3][16];
    typedef float (&refArray3)[16];
    typedef float (&refArray4)[16];
    typedef float (&refArray5)[16];

    for (int chunk_offset = 0; chunk_offset < 16; chunk_offset += ChunkSize)
      {
        //refArray1 Ak = reinterpret_cast < refArray1 > (_rA[entry][0][chunk_offset]);
        refArray1 dPdFk =
          reinterpret_cast < refArray1 > (_rdPdF[index][0][chunk_offset]);
        refArray2 Sigmak =
          reinterpret_cast < refArray2 > (_rSigma[index][0][chunk_offset]);
        refArray3 pk =
          reinterpret_cast < refArray3 > (_rp[index][chunk_offset]);
        refArray4 muk =
          reinterpret_cast < refArray4 > (_rmu[index][chunk_offset]);
        refArray5 alphak =
          reinterpret_cast < refArray5 > (_ralpha[index][chunk_offset]);

        Augmented_Rotated_Stress_Derivative < NEOHOOKEAN_TAG, __m512, float[16],
          int[16] >::Run (dPdFk, Sigmak, pk, muk, alphak);
      }

  }
};

template < int SIZE > class Augmented_Rotated_Stress_Derivative_MIC_COROTATED
{
private:
  // Generate Variables Here
  float *_local_dPdF;
  float *_local_Sigma;
  float *_local_p;
  float *_local_mu;
  float *_local_alpha;

public:
    explicit Augmented_Rotated_Stress_Derivative_MIC_COROTATED (float *dPdF_in,
                                                                float *Sigma_in,
                                                                float *p_in,
                                                                float *mu_in,
                                                                float
                                                                *alpha_in):_local_dPdF
    (dPdF_in), _local_Sigma (Sigma_in), _local_p (p_in), _local_mu (mu_in),
    _local_alpha (alpha_in)
  {
  }
  void Execute (int index)
  {
    // full array typedefs
    //typedef int (&refArray)[SIZE][3][8];
    typedef float (&fullArray1)[SIZE][12][16];
    typedef float (&fullArray2)[SIZE][3][16];
    typedef float (&fullArray3)[SIZE][16];
    typedef float (&fullArray4)[SIZE][16];
    typedef float (&fullArray5)[SIZE][16];

    // full array extractions
    //refArray _rA = reinterpret_cast < refArray >(*_A);
    fullArray1 _rdPdF = reinterpret_cast < fullArray1 > (*_local_dPdF);
    fullArray2 _rSigma = reinterpret_cast < fullArray2 > (*_local_Sigma);
    fullArray3 _rp = reinterpret_cast < fullArray3 > (*_local_p);
    fullArray4 _rmu = reinterpret_cast < fullArray4 > (*_local_mu);
    fullArray5 _ralpha = reinterpret_cast < fullArray5 > (*_local_alpha);

    const int ChunkSize = 16;
    // chunk typedef
    //typedef int (&refArray1)[3][8];
    typedef float (&refArray1)[12][16];
    typedef float (&refArray2)[3][16];
    typedef float (&refArray3)[16];
    typedef float (&refArray4)[16];
    typedef float (&refArray5)[16];

    for (int chunk_offset = 0; chunk_offset < 16; chunk_offset += ChunkSize)
      {
        //refArray1 Ak = reinterpret_cast < refArray1 > (_rA[entry][0][chunk_offset]);
        refArray1 dPdFk =
          reinterpret_cast < refArray1 > (_rdPdF[index][0][chunk_offset]);
        refArray2 Sigmak =
          reinterpret_cast < refArray2 > (_rSigma[index][0][chunk_offset]);
        refArray3 pk =
          reinterpret_cast < refArray3 > (_rp[index][chunk_offset]);
        refArray4 muk =
          reinterpret_cast < refArray4 > (_rmu[index][chunk_offset]);
        refArray5 alphak =
          reinterpret_cast < refArray5 > (_ralpha[index][chunk_offset]);

        Augmented_Rotated_Stress_Derivative < COROTATED_TAG, __m512, float[16],
          int[16] >::Run (dPdFk, Sigmak, pk, muk, alphak);
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
    std::
      cout << "Running Thread Test for Augmented_Rotated_Stress_Derivative " <<
      std::endl;

//=======================================================
//
//        DEFINE ALL VARIABLES USED BY KERNEL
//
//=======================================================
    std::cout << "\nAllocating all data: ";
    std::cout.flush ();

    start_timer ();
    typedef T (&dPdF_type)[data_size][12][16];
    dPdF_type dPdF =
      reinterpret_cast < dPdF_type >
      (*((T *) (_mm_malloc (data_size * 12 * 16 * sizeof (T), 64))));
    typedef T (&Sigma_type)[data_size][3][16];
    Sigma_type Sigma =
      reinterpret_cast < Sigma_type >
      (*((T *) (_mm_malloc (data_size * 3 * 16 * sizeof (T), 64))));
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


    for (int __a = 0; __a < data_size; __a++)
      for (int __b = 0; __b < 12; __b++)
        for (int __c = 0; __c < 16; __c++)
          {
            dPdF[__a][__b][__c] = Get_Random < float >();
    } for (int __a = 0; __a < data_size; __a++)
      for (int __b = 0; __b < 3; __b++)
        for (int __c = 0; __c < 16; __c++)
          {
            Sigma[__a][__b][__c] = Get_Random < float >();
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


      Augmented_Rotated_Stress_Derivative_SCALAR_NEOHOOKEAN < data_size >
        op ((float *) &dPdF, (float *) &Sigma, (float *) &p, (float *) &mu,
            (float *) &alpha);

      for (int t = threads; t <= threads_max;
           t += std::max < int >(((threads_max - threads) / 30), 1))
        {
          std::cout << "Running Test with " << t << " threads." << std::endl;
          MT_Streaming_Kernels::Kernel_Serial_Base_Helper <
            Augmented_Rotated_Stress_Derivative_SCALAR_NEOHOOKEAN < data_size >
            >helper (op, data_size, t);

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


      Augmented_Rotated_Stress_Derivative_SSE_NEOHOOKEAN < data_size >
        op ((float *) &dPdF, (float *) &Sigma, (float *) &p, (float *) &mu,
            (float *) &alpha);

      for (int t = threads; t <= threads_max;
           t += std::max < int >(((threads_max - threads) / 30), 1))
        {
          std::cout << "Running Test with " << t << " threads." << std::endl;
          MT_Streaming_Kernels::Kernel_Serial_Base_Helper <
            Augmented_Rotated_Stress_Derivative_SSE_NEOHOOKEAN < data_size >
            >helper (op, data_size, t);

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


      Augmented_Rotated_Stress_Derivative_AVX_NEOHOOKEAN < data_size >
        op ((float *) &dPdF, (float *) &Sigma, (float *) &p, (float *) &mu,
            (float *) &alpha);

      for (int t = threads; t <= threads_max;
           t += std::max < int >(((threads_max - threads) / 30), 1))
        {
          std::cout << "Running Test with " << t << " threads." << std::endl;
          MT_Streaming_Kernels::Kernel_Serial_Base_Helper <
            Augmented_Rotated_Stress_Derivative_AVX_NEOHOOKEAN < data_size >
            >helper (op, data_size, t);

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


      Augmented_Rotated_Stress_Derivative_NEON_NEOHOOKEAN < data_size >
        op ((float *) &dPdF, (float *) &Sigma, (float *) &p, (float *) &mu,
            (float *) &alpha);

      for (int t = threads; t <= threads_max;
           t += std::max < int >(((threads_max - threads) / 30), 1))
        {
          std::cout << "Running Test with " << t << " threads." << std::endl;
          MT_Streaming_Kernels::Kernel_Serial_Base_Helper <
            Augmented_Rotated_Stress_Derivative_NEON_NEOHOOKEAN < data_size >
            >helper (op, data_size, t);

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


      Augmented_Rotated_Stress_Derivative_MIC_NEOHOOKEAN < data_size >
        op ((float *) &dPdF, (float *) &Sigma, (float *) &p, (float *) &mu,
            (float *) &alpha);

      for (int t = threads; t <= threads_max;
           t += std::max < int >(((threads_max - threads) / 30), 1))
        {
          std::cout << "Running Test with " << t << " threads." << std::endl;
          MT_Streaming_Kernels::Kernel_Serial_Base_Helper <
            Augmented_Rotated_Stress_Derivative_MIC_NEOHOOKEAN < data_size >
            >helper (op, data_size, t);

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

    _mm_free (reinterpret_cast < void *>(dPdF));
    _mm_free (reinterpret_cast < void *>(Sigma));
    _mm_free (reinterpret_cast < void *>(p));
    _mm_free (reinterpret_cast < void *>(mu));
    _mm_free (reinterpret_cast < void *>(alpha));


  }


  {
    std::
      cout << "Running Thread Test for Augmented_Rotated_Stress_Derivative " <<
      std::endl;

//=======================================================
//
//        DEFINE ALL VARIABLES USED BY KERNEL
//
//=======================================================
    std::cout << "\nAllocating all data: ";
    std::cout.flush ();

    start_timer ();
    typedef T (&dPdF_type)[data_size][12][16];
    dPdF_type dPdF =
      reinterpret_cast < dPdF_type >
      (*((T *) (_mm_malloc (data_size * 12 * 16 * sizeof (T), 64))));
    typedef T (&Sigma_type)[data_size][3][16];
    Sigma_type Sigma =
      reinterpret_cast < Sigma_type >
      (*((T *) (_mm_malloc (data_size * 3 * 16 * sizeof (T), 64))));
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


    for (int __a = 0; __a < data_size; __a++)
      for (int __b = 0; __b < 12; __b++)
        for (int __c = 0; __c < 16; __c++)
          {
            dPdF[__a][__b][__c] = Get_Random < float >();
    } for (int __a = 0; __a < data_size; __a++)
      for (int __b = 0; __b < 3; __b++)
        for (int __c = 0; __c < 16; __c++)
          {
            Sigma[__a][__b][__c] = Get_Random < float >();
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


      Augmented_Rotated_Stress_Derivative_SCALAR_COROTATED < data_size >
        op ((float *) &dPdF, (float *) &Sigma, (float *) &p, (float *) &mu,
            (float *) &alpha);

      for (int t = threads; t <= threads_max;
           t += std::max < int >(((threads_max - threads) / 30), 1))
        {
          std::cout << "Running Test with " << t << " threads." << std::endl;
          MT_Streaming_Kernels::Kernel_Serial_Base_Helper <
            Augmented_Rotated_Stress_Derivative_SCALAR_COROTATED < data_size >
            >helper (op, data_size, t);

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


      Augmented_Rotated_Stress_Derivative_SSE_COROTATED < data_size >
        op ((float *) &dPdF, (float *) &Sigma, (float *) &p, (float *) &mu,
            (float *) &alpha);

      for (int t = threads; t <= threads_max;
           t += std::max < int >(((threads_max - threads) / 30), 1))
        {
          std::cout << "Running Test with " << t << " threads." << std::endl;
          MT_Streaming_Kernels::Kernel_Serial_Base_Helper <
            Augmented_Rotated_Stress_Derivative_SSE_COROTATED < data_size >
            >helper (op, data_size, t);

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


      Augmented_Rotated_Stress_Derivative_AVX_COROTATED < data_size >
        op ((float *) &dPdF, (float *) &Sigma, (float *) &p, (float *) &mu,
            (float *) &alpha);

      for (int t = threads; t <= threads_max;
           t += std::max < int >(((threads_max - threads) / 30), 1))
        {
          std::cout << "Running Test with " << t << " threads." << std::endl;
          MT_Streaming_Kernels::Kernel_Serial_Base_Helper <
            Augmented_Rotated_Stress_Derivative_AVX_COROTATED < data_size >
            >helper (op, data_size, t);

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


      Augmented_Rotated_Stress_Derivative_NEON_COROTATED < data_size >
        op ((float *) &dPdF, (float *) &Sigma, (float *) &p, (float *) &mu,
            (float *) &alpha);

      for (int t = threads; t <= threads_max;
           t += std::max < int >(((threads_max - threads) / 30), 1))
        {
          std::cout << "Running Test with " << t << " threads." << std::endl;
          MT_Streaming_Kernels::Kernel_Serial_Base_Helper <
            Augmented_Rotated_Stress_Derivative_NEON_COROTATED < data_size >
            >helper (op, data_size, t);

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


      Augmented_Rotated_Stress_Derivative_MIC_COROTATED < data_size >
        op ((float *) &dPdF, (float *) &Sigma, (float *) &p, (float *) &mu,
            (float *) &alpha);

      for (int t = threads; t <= threads_max;
           t += std::max < int >(((threads_max - threads) / 30), 1))
        {
          std::cout << "Running Test with " << t << " threads." << std::endl;
          MT_Streaming_Kernels::Kernel_Serial_Base_Helper <
            Augmented_Rotated_Stress_Derivative_MIC_COROTATED < data_size >
            >helper (op, data_size, t);

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

    _mm_free (reinterpret_cast < void *>(dPdF));
    _mm_free (reinterpret_cast < void *>(Sigma));
    _mm_free (reinterpret_cast < void *>(p));
    _mm_free (reinterpret_cast < void *>(mu));
    _mm_free (reinterpret_cast < void *>(alpha));


  }


  return 0;

}


#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <sys/time.h>
#include "KernelCommon.h"
struct NEOHOOKEAN_TAG;
struct COROTATED_TAG;

#include "Muscle_Tension.h"

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


template < int SIZE > class Tension_SCALAR_None
{
private:
  // Generate Variables Here
  float *_local_tension;
  float *_local_stretch;
  float *_local_activation;
  float *_local_density;
  float *_local_fiber_max_stress;

public:
    explicit Tension_SCALAR_None (float *tension_in, float *stretch_in,
                                  float *activation_in, float *density_in,
                                  float
                                  *fiber_max_stress_in):_local_tension
    (tension_in), _local_stretch (stretch_in),
    _local_activation (activation_in), _local_density (density_in),
    _local_fiber_max_stress (fiber_max_stress_in)
  {
  }
  void Execute (int index)
  {
    // full array typedefs
    //typedef int (&refArray)[SIZE][3][8];
    typedef float (&fullArray1)[SIZE][16];
    typedef float (&fullArray2)[SIZE][16];
    typedef float (&fullArray3)[SIZE][16];
    typedef float (&fullArray4)[SIZE][16];
    typedef float (&fullArray5)[SIZE][16];

    // full array extractions
    //refArray _rA = reinterpret_cast < refArray >(*_A);
    fullArray1 _rtension = reinterpret_cast < fullArray1 > (*_local_tension);
    fullArray2 _rstretch = reinterpret_cast < fullArray2 > (*_local_stretch);
    fullArray3 _ractivation =
      reinterpret_cast < fullArray3 > (*_local_activation);
    fullArray4 _rdensity = reinterpret_cast < fullArray4 > (*_local_density);
    fullArray5 _rfiber_max_stress =
      reinterpret_cast < fullArray5 > (*_local_fiber_max_stress);

    const int ChunkSize = 1;
    // chunk typedef
    //typedef int (&refArray1)[3][8];
    typedef float (&refArray1)[16];
    typedef float (&refArray2)[16];
    typedef float (&refArray3)[16];
    typedef float (&refArray4)[16];
    typedef float (&refArray5)[16];

    for (int chunk_offset = 0; chunk_offset < 16; chunk_offset += ChunkSize)
      {
        //refArray1 Ak = reinterpret_cast < refArray1 > (_rA[entry][0][chunk_offset]);
        refArray1 tensionk =
          reinterpret_cast < refArray1 > (_rtension[index][chunk_offset]);
        refArray2 stretchk =
          reinterpret_cast < refArray2 > (_rstretch[index][chunk_offset]);
        refArray3 activationk =
          reinterpret_cast < refArray3 > (_ractivation[index][chunk_offset]);
        refArray4 densityk =
          reinterpret_cast < refArray4 > (_rdensity[index][chunk_offset]);
        refArray5 fiber_max_stressk =
          reinterpret_cast < refArray5 >
          (_rfiber_max_stress[index][chunk_offset]);

        Tension < float, float[16], int[16] > (tensionk, stretchk, activationk,
                                               densityk, fiber_max_stressk);
      }

  }
};

template < int SIZE > class Tension_Derivative_SCALAR_None
{
private:
  // Generate Variables Here
  float *_local_tension_derivative;
  float *_local_stretch;
  float *_local_activation;
  float *_local_density;
  float *_local_fiber_max_stress;

public:
    explicit Tension_Derivative_SCALAR_None (float *tension_derivative_in,
                                             float *stretch_in,
                                             float *activation_in,
                                             float *density_in,
                                             float
                                             *fiber_max_stress_in):_local_tension_derivative
    (tension_derivative_in), _local_stretch (stretch_in),
    _local_activation (activation_in), _local_density (density_in),
    _local_fiber_max_stress (fiber_max_stress_in)
  {
  }
  void Execute (int index)
  {
    // full array typedefs
    //typedef int (&refArray)[SIZE][3][8];
    typedef float (&fullArray1)[SIZE][16];
    typedef float (&fullArray2)[SIZE][16];
    typedef float (&fullArray3)[SIZE][16];
    typedef float (&fullArray4)[SIZE][16];
    typedef float (&fullArray5)[SIZE][16];

    // full array extractions
    //refArray _rA = reinterpret_cast < refArray >(*_A);
    fullArray1 _rtension_derivative =
      reinterpret_cast < fullArray1 > (*_local_tension_derivative);
    fullArray2 _rstretch = reinterpret_cast < fullArray2 > (*_local_stretch);
    fullArray3 _ractivation =
      reinterpret_cast < fullArray3 > (*_local_activation);
    fullArray4 _rdensity = reinterpret_cast < fullArray4 > (*_local_density);
    fullArray5 _rfiber_max_stress =
      reinterpret_cast < fullArray5 > (*_local_fiber_max_stress);

    const int ChunkSize = 1;
    // chunk typedef
    //typedef int (&refArray1)[3][8];
    typedef float (&refArray1)[16];
    typedef float (&refArray2)[16];
    typedef float (&refArray3)[16];
    typedef float (&refArray4)[16];
    typedef float (&refArray5)[16];

    for (int chunk_offset = 0; chunk_offset < 16; chunk_offset += ChunkSize)
      {
        //refArray1 Ak = reinterpret_cast < refArray1 > (_rA[entry][0][chunk_offset]);
        refArray1 tension_derivativek =
          reinterpret_cast < refArray1 >
          (_rtension_derivative[index][chunk_offset]);
        refArray2 stretchk =
          reinterpret_cast < refArray2 > (_rstretch[index][chunk_offset]);
        refArray3 activationk =
          reinterpret_cast < refArray3 > (_ractivation[index][chunk_offset]);
        refArray4 densityk =
          reinterpret_cast < refArray4 > (_rdensity[index][chunk_offset]);
        refArray5 fiber_max_stressk =
          reinterpret_cast < refArray5 >
          (_rfiber_max_stress[index][chunk_offset]);

        Tension_Derivative < float, float[16], int[16] > (tension_derivativek,
                                                          stretchk, activationk,
                                                          densityk,
                                                          fiber_max_stressk);
      }

  }
};


#ifdef ENABLE_SSE_INSTRUCTION_SET

template < int SIZE > class Tension_SSE_None
{
private:
  // Generate Variables Here
  float *_local_tension;
  float *_local_stretch;
  float *_local_activation;
  float *_local_density;
  float *_local_fiber_max_stress;

public:
    explicit Tension_SSE_None (float *tension_in, float *stretch_in,
                               float *activation_in, float *density_in,
                               float
                               *fiber_max_stress_in):_local_tension
    (tension_in), _local_stretch (stretch_in),
    _local_activation (activation_in), _local_density (density_in),
    _local_fiber_max_stress (fiber_max_stress_in)
  {
  }
  void Execute (int index)
  {
    // full array typedefs
    //typedef int (&refArray)[SIZE][3][8];
    typedef float (&fullArray1)[SIZE][16];
    typedef float (&fullArray2)[SIZE][16];
    typedef float (&fullArray3)[SIZE][16];
    typedef float (&fullArray4)[SIZE][16];
    typedef float (&fullArray5)[SIZE][16];

    // full array extractions
    //refArray _rA = reinterpret_cast < refArray >(*_A);
    fullArray1 _rtension = reinterpret_cast < fullArray1 > (*_local_tension);
    fullArray2 _rstretch = reinterpret_cast < fullArray2 > (*_local_stretch);
    fullArray3 _ractivation =
      reinterpret_cast < fullArray3 > (*_local_activation);
    fullArray4 _rdensity = reinterpret_cast < fullArray4 > (*_local_density);
    fullArray5 _rfiber_max_stress =
      reinterpret_cast < fullArray5 > (*_local_fiber_max_stress);

    const int ChunkSize = 4;
    // chunk typedef
    //typedef int (&refArray1)[3][8];
    typedef float (&refArray1)[16];
    typedef float (&refArray2)[16];
    typedef float (&refArray3)[16];
    typedef float (&refArray4)[16];
    typedef float (&refArray5)[16];

    for (int chunk_offset = 0; chunk_offset < 16; chunk_offset += ChunkSize)
      {
        //refArray1 Ak = reinterpret_cast < refArray1 > (_rA[entry][0][chunk_offset]);
        refArray1 tensionk =
          reinterpret_cast < refArray1 > (_rtension[index][chunk_offset]);
        refArray2 stretchk =
          reinterpret_cast < refArray2 > (_rstretch[index][chunk_offset]);
        refArray3 activationk =
          reinterpret_cast < refArray3 > (_ractivation[index][chunk_offset]);
        refArray4 densityk =
          reinterpret_cast < refArray4 > (_rdensity[index][chunk_offset]);
        refArray5 fiber_max_stressk =
          reinterpret_cast < refArray5 >
          (_rfiber_max_stress[index][chunk_offset]);

        Tension < __m128, float[16], int[16] > (tensionk, stretchk, activationk,
                                                densityk, fiber_max_stressk);
      }

  }
};

template < int SIZE > class Tension_Derivative_SSE_None
{
private:
  // Generate Variables Here
  float *_local_tension_derivative;
  float *_local_stretch;
  float *_local_activation;
  float *_local_density;
  float *_local_fiber_max_stress;

public:
    explicit Tension_Derivative_SSE_None (float *tension_derivative_in,
                                          float *stretch_in,
                                          float *activation_in,
                                          float *density_in,
                                          float
                                          *fiber_max_stress_in):_local_tension_derivative
    (tension_derivative_in), _local_stretch (stretch_in),
    _local_activation (activation_in), _local_density (density_in),
    _local_fiber_max_stress (fiber_max_stress_in)
  {
  }
  void Execute (int index)
  {
    // full array typedefs
    //typedef int (&refArray)[SIZE][3][8];
    typedef float (&fullArray1)[SIZE][16];
    typedef float (&fullArray2)[SIZE][16];
    typedef float (&fullArray3)[SIZE][16];
    typedef float (&fullArray4)[SIZE][16];
    typedef float (&fullArray5)[SIZE][16];

    // full array extractions
    //refArray _rA = reinterpret_cast < refArray >(*_A);
    fullArray1 _rtension_derivative =
      reinterpret_cast < fullArray1 > (*_local_tension_derivative);
    fullArray2 _rstretch = reinterpret_cast < fullArray2 > (*_local_stretch);
    fullArray3 _ractivation =
      reinterpret_cast < fullArray3 > (*_local_activation);
    fullArray4 _rdensity = reinterpret_cast < fullArray4 > (*_local_density);
    fullArray5 _rfiber_max_stress =
      reinterpret_cast < fullArray5 > (*_local_fiber_max_stress);

    const int ChunkSize = 4;
    // chunk typedef
    //typedef int (&refArray1)[3][8];
    typedef float (&refArray1)[16];
    typedef float (&refArray2)[16];
    typedef float (&refArray3)[16];
    typedef float (&refArray4)[16];
    typedef float (&refArray5)[16];

    for (int chunk_offset = 0; chunk_offset < 16; chunk_offset += ChunkSize)
      {
        //refArray1 Ak = reinterpret_cast < refArray1 > (_rA[entry][0][chunk_offset]);
        refArray1 tension_derivativek =
          reinterpret_cast < refArray1 >
          (_rtension_derivative[index][chunk_offset]);
        refArray2 stretchk =
          reinterpret_cast < refArray2 > (_rstretch[index][chunk_offset]);
        refArray3 activationk =
          reinterpret_cast < refArray3 > (_ractivation[index][chunk_offset]);
        refArray4 densityk =
          reinterpret_cast < refArray4 > (_rdensity[index][chunk_offset]);
        refArray5 fiber_max_stressk =
          reinterpret_cast < refArray5 >
          (_rfiber_max_stress[index][chunk_offset]);

        Tension_Derivative < __m128, float[16], int[16] > (tension_derivativek,
                                                           stretchk,
                                                           activationk,
                                                           densityk,
                                                           fiber_max_stressk);
      }

  }
};

#endif

#ifdef ENABLE_AVX_INSTRUCTION_SET

template < int SIZE > class Tension_AVX_None
{
private:
  // Generate Variables Here
  float *_local_tension;
  float *_local_stretch;
  float *_local_activation;
  float *_local_density;
  float *_local_fiber_max_stress;

public:
    explicit Tension_AVX_None (float *tension_in, float *stretch_in,
                               float *activation_in, float *density_in,
                               float
                               *fiber_max_stress_in):_local_tension
    (tension_in), _local_stretch (stretch_in),
    _local_activation (activation_in), _local_density (density_in),
    _local_fiber_max_stress (fiber_max_stress_in)
  {
  }
  void Execute (int index)
  {
    // full array typedefs
    //typedef int (&refArray)[SIZE][3][8];
    typedef float (&fullArray1)[SIZE][16];
    typedef float (&fullArray2)[SIZE][16];
    typedef float (&fullArray3)[SIZE][16];
    typedef float (&fullArray4)[SIZE][16];
    typedef float (&fullArray5)[SIZE][16];

    // full array extractions
    //refArray _rA = reinterpret_cast < refArray >(*_A);
    fullArray1 _rtension = reinterpret_cast < fullArray1 > (*_local_tension);
    fullArray2 _rstretch = reinterpret_cast < fullArray2 > (*_local_stretch);
    fullArray3 _ractivation =
      reinterpret_cast < fullArray3 > (*_local_activation);
    fullArray4 _rdensity = reinterpret_cast < fullArray4 > (*_local_density);
    fullArray5 _rfiber_max_stress =
      reinterpret_cast < fullArray5 > (*_local_fiber_max_stress);

    const int ChunkSize = 8;
    // chunk typedef
    //typedef int (&refArray1)[3][8];
    typedef float (&refArray1)[16];
    typedef float (&refArray2)[16];
    typedef float (&refArray3)[16];
    typedef float (&refArray4)[16];
    typedef float (&refArray5)[16];

    for (int chunk_offset = 0; chunk_offset < 16; chunk_offset += ChunkSize)
      {
        //refArray1 Ak = reinterpret_cast < refArray1 > (_rA[entry][0][chunk_offset]);
        refArray1 tensionk =
          reinterpret_cast < refArray1 > (_rtension[index][chunk_offset]);
        refArray2 stretchk =
          reinterpret_cast < refArray2 > (_rstretch[index][chunk_offset]);
        refArray3 activationk =
          reinterpret_cast < refArray3 > (_ractivation[index][chunk_offset]);
        refArray4 densityk =
          reinterpret_cast < refArray4 > (_rdensity[index][chunk_offset]);
        refArray5 fiber_max_stressk =
          reinterpret_cast < refArray5 >
          (_rfiber_max_stress[index][chunk_offset]);

        Tension < __m256, float[16], int[16] > (tensionk, stretchk, activationk,
                                                densityk, fiber_max_stressk);
      }

  }
};

template < int SIZE > class Tension_Derivative_AVX_None
{
private:
  // Generate Variables Here
  float *_local_tension_derivative;
  float *_local_stretch;
  float *_local_activation;
  float *_local_density;
  float *_local_fiber_max_stress;

public:
    explicit Tension_Derivative_AVX_None (float *tension_derivative_in,
                                          float *stretch_in,
                                          float *activation_in,
                                          float *density_in,
                                          float
                                          *fiber_max_stress_in):_local_tension_derivative
    (tension_derivative_in), _local_stretch (stretch_in),
    _local_activation (activation_in), _local_density (density_in),
    _local_fiber_max_stress (fiber_max_stress_in)
  {
  }
  void Execute (int index)
  {
    // full array typedefs
    //typedef int (&refArray)[SIZE][3][8];
    typedef float (&fullArray1)[SIZE][16];
    typedef float (&fullArray2)[SIZE][16];
    typedef float (&fullArray3)[SIZE][16];
    typedef float (&fullArray4)[SIZE][16];
    typedef float (&fullArray5)[SIZE][16];

    // full array extractions
    //refArray _rA = reinterpret_cast < refArray >(*_A);
    fullArray1 _rtension_derivative =
      reinterpret_cast < fullArray1 > (*_local_tension_derivative);
    fullArray2 _rstretch = reinterpret_cast < fullArray2 > (*_local_stretch);
    fullArray3 _ractivation =
      reinterpret_cast < fullArray3 > (*_local_activation);
    fullArray4 _rdensity = reinterpret_cast < fullArray4 > (*_local_density);
    fullArray5 _rfiber_max_stress =
      reinterpret_cast < fullArray5 > (*_local_fiber_max_stress);

    const int ChunkSize = 8;
    // chunk typedef
    //typedef int (&refArray1)[3][8];
    typedef float (&refArray1)[16];
    typedef float (&refArray2)[16];
    typedef float (&refArray3)[16];
    typedef float (&refArray4)[16];
    typedef float (&refArray5)[16];

    for (int chunk_offset = 0; chunk_offset < 16; chunk_offset += ChunkSize)
      {
        //refArray1 Ak = reinterpret_cast < refArray1 > (_rA[entry][0][chunk_offset]);
        refArray1 tension_derivativek =
          reinterpret_cast < refArray1 >
          (_rtension_derivative[index][chunk_offset]);
        refArray2 stretchk =
          reinterpret_cast < refArray2 > (_rstretch[index][chunk_offset]);
        refArray3 activationk =
          reinterpret_cast < refArray3 > (_ractivation[index][chunk_offset]);
        refArray4 densityk =
          reinterpret_cast < refArray4 > (_rdensity[index][chunk_offset]);
        refArray5 fiber_max_stressk =
          reinterpret_cast < refArray5 >
          (_rfiber_max_stress[index][chunk_offset]);

        Tension_Derivative < __m256, float[16], int[16] > (tension_derivativek,
                                                           stretchk,
                                                           activationk,
                                                           densityk,
                                                           fiber_max_stressk);
      }

  }
};

#endif

#ifdef ENABLE_NEON_INSTRUCTION_SET

template < int SIZE > class Tension_NEON_None
{
private:
  // Generate Variables Here
  float *_local_tension;
  float *_local_stretch;
  float *_local_activation;
  float *_local_density;
  float *_local_fiber_max_stress;

public:
    explicit Tension_NEON_None (float *tension_in, float *stretch_in,
                                float *activation_in, float *density_in,
                                float
                                *fiber_max_stress_in):_local_tension
    (tension_in), _local_stretch (stretch_in),
    _local_activation (activation_in), _local_density (density_in),
    _local_fiber_max_stress (fiber_max_stress_in)
  {
  }
  void Execute (int index)
  {
    // full array typedefs
    //typedef int (&refArray)[SIZE][3][8];
    typedef float (&fullArray1)[SIZE][16];
    typedef float (&fullArray2)[SIZE][16];
    typedef float (&fullArray3)[SIZE][16];
    typedef float (&fullArray4)[SIZE][16];
    typedef float (&fullArray5)[SIZE][16];

    // full array extractions
    //refArray _rA = reinterpret_cast < refArray >(*_A);
    fullArray1 _rtension = reinterpret_cast < fullArray1 > (*_local_tension);
    fullArray2 _rstretch = reinterpret_cast < fullArray2 > (*_local_stretch);
    fullArray3 _ractivation =
      reinterpret_cast < fullArray3 > (*_local_activation);
    fullArray4 _rdensity = reinterpret_cast < fullArray4 > (*_local_density);
    fullArray5 _rfiber_max_stress =
      reinterpret_cast < fullArray5 > (*_local_fiber_max_stress);

    const int ChunkSize = 4;
    // chunk typedef
    //typedef int (&refArray1)[3][8];
    typedef float (&refArray1)[16];
    typedef float (&refArray2)[16];
    typedef float (&refArray3)[16];
    typedef float (&refArray4)[16];
    typedef float (&refArray5)[16];

    for (int chunk_offset = 0; chunk_offset < 16; chunk_offset += ChunkSize)
      {
        //refArray1 Ak = reinterpret_cast < refArray1 > (_rA[entry][0][chunk_offset]);
        refArray1 tensionk =
          reinterpret_cast < refArray1 > (_rtension[index][chunk_offset]);
        refArray2 stretchk =
          reinterpret_cast < refArray2 > (_rstretch[index][chunk_offset]);
        refArray3 activationk =
          reinterpret_cast < refArray3 > (_ractivation[index][chunk_offset]);
        refArray4 densityk =
          reinterpret_cast < refArray4 > (_rdensity[index][chunk_offset]);
        refArray5 fiber_max_stressk =
          reinterpret_cast < refArray5 >
          (_rfiber_max_stress[index][chunk_offset]);

        Tension < float32x4_t, float[16], int[16] > (tensionk, stretchk,
                                                     activationk, densityk,
                                                     fiber_max_stressk);
      }

  }
};

template < int SIZE > class Tension_Derivative_NEON_None
{
private:
  // Generate Variables Here
  float *_local_tension_derivative;
  float *_local_stretch;
  float *_local_activation;
  float *_local_density;
  float *_local_fiber_max_stress;

public:
    explicit Tension_Derivative_NEON_None (float *tension_derivative_in,
                                           float *stretch_in,
                                           float *activation_in,
                                           float *density_in,
                                           float
                                           *fiber_max_stress_in):_local_tension_derivative
    (tension_derivative_in), _local_stretch (stretch_in),
    _local_activation (activation_in), _local_density (density_in),
    _local_fiber_max_stress (fiber_max_stress_in)
  {
  }
  void Execute (int index)
  {
    // full array typedefs
    //typedef int (&refArray)[SIZE][3][8];
    typedef float (&fullArray1)[SIZE][16];
    typedef float (&fullArray2)[SIZE][16];
    typedef float (&fullArray3)[SIZE][16];
    typedef float (&fullArray4)[SIZE][16];
    typedef float (&fullArray5)[SIZE][16];

    // full array extractions
    //refArray _rA = reinterpret_cast < refArray >(*_A);
    fullArray1 _rtension_derivative =
      reinterpret_cast < fullArray1 > (*_local_tension_derivative);
    fullArray2 _rstretch = reinterpret_cast < fullArray2 > (*_local_stretch);
    fullArray3 _ractivation =
      reinterpret_cast < fullArray3 > (*_local_activation);
    fullArray4 _rdensity = reinterpret_cast < fullArray4 > (*_local_density);
    fullArray5 _rfiber_max_stress =
      reinterpret_cast < fullArray5 > (*_local_fiber_max_stress);

    const int ChunkSize = 4;
    // chunk typedef
    //typedef int (&refArray1)[3][8];
    typedef float (&refArray1)[16];
    typedef float (&refArray2)[16];
    typedef float (&refArray3)[16];
    typedef float (&refArray4)[16];
    typedef float (&refArray5)[16];

    for (int chunk_offset = 0; chunk_offset < 16; chunk_offset += ChunkSize)
      {
        //refArray1 Ak = reinterpret_cast < refArray1 > (_rA[entry][0][chunk_offset]);
        refArray1 tension_derivativek =
          reinterpret_cast < refArray1 >
          (_rtension_derivative[index][chunk_offset]);
        refArray2 stretchk =
          reinterpret_cast < refArray2 > (_rstretch[index][chunk_offset]);
        refArray3 activationk =
          reinterpret_cast < refArray3 > (_ractivation[index][chunk_offset]);
        refArray4 densityk =
          reinterpret_cast < refArray4 > (_rdensity[index][chunk_offset]);
        refArray5 fiber_max_stressk =
          reinterpret_cast < refArray5 >
          (_rfiber_max_stress[index][chunk_offset]);

        Tension_Derivative < float32x4_t, float[16],
          int[16] > (tension_derivativek, stretchk, activationk, densityk,
                     fiber_max_stressk);
      }

  }
};

#endif

#ifdef ENABLE_MIC_INSTRUCTION_SET

template < int SIZE > class Tension_MIC_None
{
private:
  // Generate Variables Here
  float *_local_tension;
  float *_local_stretch;
  float *_local_activation;
  float *_local_density;
  float *_local_fiber_max_stress;

public:
    explicit Tension_MIC_None (float *tension_in, float *stretch_in,
                               float *activation_in, float *density_in,
                               float
                               *fiber_max_stress_in):_local_tension
    (tension_in), _local_stretch (stretch_in),
    _local_activation (activation_in), _local_density (density_in),
    _local_fiber_max_stress (fiber_max_stress_in)
  {
  }
  void Execute (int index)
  {
    // full array typedefs
    //typedef int (&refArray)[SIZE][3][8];
    typedef float (&fullArray1)[SIZE][16];
    typedef float (&fullArray2)[SIZE][16];
    typedef float (&fullArray3)[SIZE][16];
    typedef float (&fullArray4)[SIZE][16];
    typedef float (&fullArray5)[SIZE][16];

    // full array extractions
    //refArray _rA = reinterpret_cast < refArray >(*_A);
    fullArray1 _rtension = reinterpret_cast < fullArray1 > (*_local_tension);
    fullArray2 _rstretch = reinterpret_cast < fullArray2 > (*_local_stretch);
    fullArray3 _ractivation =
      reinterpret_cast < fullArray3 > (*_local_activation);
    fullArray4 _rdensity = reinterpret_cast < fullArray4 > (*_local_density);
    fullArray5 _rfiber_max_stress =
      reinterpret_cast < fullArray5 > (*_local_fiber_max_stress);

    const int ChunkSize = 16;
    // chunk typedef
    //typedef int (&refArray1)[3][8];
    typedef float (&refArray1)[16];
    typedef float (&refArray2)[16];
    typedef float (&refArray3)[16];
    typedef float (&refArray4)[16];
    typedef float (&refArray5)[16];

    for (int chunk_offset = 0; chunk_offset < 16; chunk_offset += ChunkSize)
      {
        //refArray1 Ak = reinterpret_cast < refArray1 > (_rA[entry][0][chunk_offset]);
        refArray1 tensionk =
          reinterpret_cast < refArray1 > (_rtension[index][chunk_offset]);
        refArray2 stretchk =
          reinterpret_cast < refArray2 > (_rstretch[index][chunk_offset]);
        refArray3 activationk =
          reinterpret_cast < refArray3 > (_ractivation[index][chunk_offset]);
        refArray4 densityk =
          reinterpret_cast < refArray4 > (_rdensity[index][chunk_offset]);
        refArray5 fiber_max_stressk =
          reinterpret_cast < refArray5 >
          (_rfiber_max_stress[index][chunk_offset]);

        Tension < __m512, float[16], int[16] > (tensionk, stretchk, activationk,
                                                densityk, fiber_max_stressk);
      }

  }
};

template < int SIZE > class Tension_Derivative_MIC_None
{
private:
  // Generate Variables Here
  float *_local_tension_derivative;
  float *_local_stretch;
  float *_local_activation;
  float *_local_density;
  float *_local_fiber_max_stress;

public:
    explicit Tension_Derivative_MIC_None (float *tension_derivative_in,
                                          float *stretch_in,
                                          float *activation_in,
                                          float *density_in,
                                          float
                                          *fiber_max_stress_in):_local_tension_derivative
    (tension_derivative_in), _local_stretch (stretch_in),
    _local_activation (activation_in), _local_density (density_in),
    _local_fiber_max_stress (fiber_max_stress_in)
  {
  }
  void Execute (int index)
  {
    // full array typedefs
    //typedef int (&refArray)[SIZE][3][8];
    typedef float (&fullArray1)[SIZE][16];
    typedef float (&fullArray2)[SIZE][16];
    typedef float (&fullArray3)[SIZE][16];
    typedef float (&fullArray4)[SIZE][16];
    typedef float (&fullArray5)[SIZE][16];

    // full array extractions
    //refArray _rA = reinterpret_cast < refArray >(*_A);
    fullArray1 _rtension_derivative =
      reinterpret_cast < fullArray1 > (*_local_tension_derivative);
    fullArray2 _rstretch = reinterpret_cast < fullArray2 > (*_local_stretch);
    fullArray3 _ractivation =
      reinterpret_cast < fullArray3 > (*_local_activation);
    fullArray4 _rdensity = reinterpret_cast < fullArray4 > (*_local_density);
    fullArray5 _rfiber_max_stress =
      reinterpret_cast < fullArray5 > (*_local_fiber_max_stress);

    const int ChunkSize = 16;
    // chunk typedef
    //typedef int (&refArray1)[3][8];
    typedef float (&refArray1)[16];
    typedef float (&refArray2)[16];
    typedef float (&refArray3)[16];
    typedef float (&refArray4)[16];
    typedef float (&refArray5)[16];

    for (int chunk_offset = 0; chunk_offset < 16; chunk_offset += ChunkSize)
      {
        //refArray1 Ak = reinterpret_cast < refArray1 > (_rA[entry][0][chunk_offset]);
        refArray1 tension_derivativek =
          reinterpret_cast < refArray1 >
          (_rtension_derivative[index][chunk_offset]);
        refArray2 stretchk =
          reinterpret_cast < refArray2 > (_rstretch[index][chunk_offset]);
        refArray3 activationk =
          reinterpret_cast < refArray3 > (_ractivation[index][chunk_offset]);
        refArray4 densityk =
          reinterpret_cast < refArray4 > (_rdensity[index][chunk_offset]);
        refArray5 fiber_max_stressk =
          reinterpret_cast < refArray5 >
          (_rfiber_max_stress[index][chunk_offset]);

        Tension_Derivative < __m512, float[16], int[16] > (tension_derivativek,
                                                           stretchk,
                                                           activationk,
                                                           densityk,
                                                           fiber_max_stressk);
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
    std::cout << "Running Thread Test for Tension " << std::endl;

//=======================================================
//
//        DEFINE ALL VARIABLES USED BY KERNEL
//
//=======================================================
    std::cout << "\nAllocating all data: ";
    std::cout.flush ();

    start_timer ();
    typedef T (&tension_type)[data_size][16];
    tension_type tension =
      reinterpret_cast < tension_type >
      (*((T *) (_mm_malloc (data_size * 16 * sizeof (T), 64))));
    typedef T (&stretch_type)[data_size][16];
    stretch_type stretch =
      reinterpret_cast < stretch_type >
      (*((T *) (_mm_malloc (data_size * 16 * sizeof (T), 64))));
    typedef T (&activation_type)[data_size][16];
    activation_type activation =
      reinterpret_cast < activation_type >
      (*((T *) (_mm_malloc (data_size * 16 * sizeof (T), 64))));
    typedef T (&density_type)[data_size][16];
    density_type density =
      reinterpret_cast < density_type >
      (*((T *) (_mm_malloc (data_size * 16 * sizeof (T), 64))));
    typedef T (&fiber_max_stress_type)[data_size][16];
    fiber_max_stress_type fiber_max_stress =
      reinterpret_cast < fiber_max_stress_type >
      (*((T *) (_mm_malloc (data_size * 16 * sizeof (T), 64))));


    for (int __a = 0; __a < data_size; __a++)
      for (int __b = 0; __b < 16; __b++)
        {
          tension[__a][__b] = Get_Random < float >();
    } for (int __a = 0; __a < data_size; __a++)
      for (int __b = 0; __b < 16; __b++)
        {
          stretch[__a][__b] = Get_Random < float >();
    } for (int __a = 0; __a < data_size; __a++)
      for (int __b = 0; __b < 16; __b++)
        {
          activation[__a][__b] = Get_Random < float >();
    } for (int __a = 0; __a < data_size; __a++)
      for (int __b = 0; __b < 16; __b++)
        {
          density[__a][__b] = Get_Random < float >();
    } for (int __a = 0; __a < data_size; __a++)
      for (int __b = 0; __b < 16; __b++)
        {
          fiber_max_stress[__a][__b] = Get_Random < float >();
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


      Tension_SCALAR_None < data_size > op ((float *) &tension,
                                            (float *) &stretch,
                                            (float *) &activation,
                                            (float *) &density,
                                            (float *) &fiber_max_stress);

      for (int t = threads; t <= threads_max;
           t += std::max < int >(((threads_max - threads) / 30), 1))
        {
          std::cout << "Running Test with " << t << " threads." << std::endl;
          MT_Streaming_Kernels::Kernel_Serial_Base_Helper <
            Tension_SCALAR_None < data_size > >helper (op, data_size, t);

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


      Tension_SSE_None < data_size > op ((float *) &tension, (float *) &stretch,
                                         (float *) &activation,
                                         (float *) &density,
                                         (float *) &fiber_max_stress);

      for (int t = threads; t <= threads_max;
           t += std::max < int >(((threads_max - threads) / 30), 1))
        {
          std::cout << "Running Test with " << t << " threads." << std::endl;
          MT_Streaming_Kernels::Kernel_Serial_Base_Helper < Tension_SSE_None <
            data_size > >helper (op, data_size, t);

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


      Tension_AVX_None < data_size > op ((float *) &tension, (float *) &stretch,
                                         (float *) &activation,
                                         (float *) &density,
                                         (float *) &fiber_max_stress);

      for (int t = threads; t <= threads_max;
           t += std::max < int >(((threads_max - threads) / 30), 1))
        {
          std::cout << "Running Test with " << t << " threads." << std::endl;
          MT_Streaming_Kernels::Kernel_Serial_Base_Helper < Tension_AVX_None <
            data_size > >helper (op, data_size, t);

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


      Tension_NEON_None < data_size > op ((float *) &tension,
                                          (float *) &stretch,
                                          (float *) &activation,
                                          (float *) &density,
                                          (float *) &fiber_max_stress);

      for (int t = threads; t <= threads_max;
           t += std::max < int >(((threads_max - threads) / 30), 1))
        {
          std::cout << "Running Test with " << t << " threads." << std::endl;
          MT_Streaming_Kernels::Kernel_Serial_Base_Helper < Tension_NEON_None <
            data_size > >helper (op, data_size, t);

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


      Tension_MIC_None < data_size > op ((float *) &tension, (float *) &stretch,
                                         (float *) &activation,
                                         (float *) &density,
                                         (float *) &fiber_max_stress);

      for (int t = threads; t <= threads_max;
           t += std::max < int >(((threads_max - threads) / 30), 1))
        {
          std::cout << "Running Test with " << t << " threads." << std::endl;
          MT_Streaming_Kernels::Kernel_Serial_Base_Helper < Tension_MIC_None <
            data_size > >helper (op, data_size, t);

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

    _mm_free (reinterpret_cast < void *>(tension));
    _mm_free (reinterpret_cast < void *>(stretch));
    _mm_free (reinterpret_cast < void *>(activation));
    _mm_free (reinterpret_cast < void *>(density));
    _mm_free (reinterpret_cast < void *>(fiber_max_stress));


  }


  {
    std::cout << "Running Thread Test for Tension_Derivative " << std::endl;

//=======================================================
//
//        DEFINE ALL VARIABLES USED BY KERNEL
//
//=======================================================
    std::cout << "\nAllocating all data: ";
    std::cout.flush ();

    start_timer ();
    typedef T (&tension_derivative_type)[data_size][16];
    tension_derivative_type tension_derivative =
      reinterpret_cast < tension_derivative_type >
      (*((T *) (_mm_malloc (data_size * 16 * sizeof (T), 64))));
    typedef T (&stretch_type)[data_size][16];
    stretch_type stretch =
      reinterpret_cast < stretch_type >
      (*((T *) (_mm_malloc (data_size * 16 * sizeof (T), 64))));
    typedef T (&activation_type)[data_size][16];
    activation_type activation =
      reinterpret_cast < activation_type >
      (*((T *) (_mm_malloc (data_size * 16 * sizeof (T), 64))));
    typedef T (&density_type)[data_size][16];
    density_type density =
      reinterpret_cast < density_type >
      (*((T *) (_mm_malloc (data_size * 16 * sizeof (T), 64))));
    typedef T (&fiber_max_stress_type)[data_size][16];
    fiber_max_stress_type fiber_max_stress =
      reinterpret_cast < fiber_max_stress_type >
      (*((T *) (_mm_malloc (data_size * 16 * sizeof (T), 64))));


    for (int __a = 0; __a < data_size; __a++)
      for (int __b = 0; __b < 16; __b++)
        {
          tension_derivative[__a][__b] = Get_Random < float >();
    } for (int __a = 0; __a < data_size; __a++)
      for (int __b = 0; __b < 16; __b++)
        {
          stretch[__a][__b] = Get_Random < float >();
    } for (int __a = 0; __a < data_size; __a++)
      for (int __b = 0; __b < 16; __b++)
        {
          activation[__a][__b] = Get_Random < float >();
    } for (int __a = 0; __a < data_size; __a++)
      for (int __b = 0; __b < 16; __b++)
        {
          density[__a][__b] = Get_Random < float >();
    } for (int __a = 0; __a < data_size; __a++)
      for (int __b = 0; __b < 16; __b++)
        {
          fiber_max_stress[__a][__b] = Get_Random < float >();
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


      Tension_Derivative_SCALAR_None < data_size >
        op ((float *) &tension_derivative, (float *) &stretch,
            (float *) &activation, (float *) &density,
            (float *) &fiber_max_stress);

      for (int t = threads; t <= threads_max;
           t += std::max < int >(((threads_max - threads) / 30), 1))
        {
          std::cout << "Running Test with " << t << " threads." << std::endl;
          MT_Streaming_Kernels::Kernel_Serial_Base_Helper <
            Tension_Derivative_SCALAR_None < data_size > >helper (op, data_size,
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


      Tension_Derivative_SSE_None < data_size >
        op ((float *) &tension_derivative, (float *) &stretch,
            (float *) &activation, (float *) &density,
            (float *) &fiber_max_stress);

      for (int t = threads; t <= threads_max;
           t += std::max < int >(((threads_max - threads) / 30), 1))
        {
          std::cout << "Running Test with " << t << " threads." << std::endl;
          MT_Streaming_Kernels::Kernel_Serial_Base_Helper <
            Tension_Derivative_SSE_None < data_size > >helper (op, data_size,
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


      Tension_Derivative_AVX_None < data_size >
        op ((float *) &tension_derivative, (float *) &stretch,
            (float *) &activation, (float *) &density,
            (float *) &fiber_max_stress);

      for (int t = threads; t <= threads_max;
           t += std::max < int >(((threads_max - threads) / 30), 1))
        {
          std::cout << "Running Test with " << t << " threads." << std::endl;
          MT_Streaming_Kernels::Kernel_Serial_Base_Helper <
            Tension_Derivative_AVX_None < data_size > >helper (op, data_size,
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


      Tension_Derivative_NEON_None < data_size >
        op ((float *) &tension_derivative, (float *) &stretch,
            (float *) &activation, (float *) &density,
            (float *) &fiber_max_stress);

      for (int t = threads; t <= threads_max;
           t += std::max < int >(((threads_max - threads) / 30), 1))
        {
          std::cout << "Running Test with " << t << " threads." << std::endl;
          MT_Streaming_Kernels::Kernel_Serial_Base_Helper <
            Tension_Derivative_NEON_None < data_size > >helper (op, data_size,
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


      Tension_Derivative_MIC_None < data_size >
        op ((float *) &tension_derivative, (float *) &stretch,
            (float *) &activation, (float *) &density,
            (float *) &fiber_max_stress);

      for (int t = threads; t <= threads_max;
           t += std::max < int >(((threads_max - threads) / 30), 1))
        {
          std::cout << "Running Test with " << t << " threads." << std::endl;
          MT_Streaming_Kernels::Kernel_Serial_Base_Helper <
            Tension_Derivative_MIC_None < data_size > >helper (op, data_size,
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

    _mm_free (reinterpret_cast < void *>(tension_derivative));
    _mm_free (reinterpret_cast < void *>(stretch));
    _mm_free (reinterpret_cast < void *>(activation));
    _mm_free (reinterpret_cast < void *>(density));
    _mm_free (reinterpret_cast < void *>(fiber_max_stress));


  }


  return 0;

}

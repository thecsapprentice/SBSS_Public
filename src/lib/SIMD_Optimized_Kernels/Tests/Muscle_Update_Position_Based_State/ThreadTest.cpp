
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <sys/time.h>
#include "KernelCommon.h"
struct NEOHOOKEAN_TAG;
struct COROTATED_TAG;

#include "Muscle_Update_Position_Based_State.h"

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


template < int SIZE > class Muscle_Update_Position_Based_State_SCALAR_None
{
private:
  // Generate Variables Here
  float *_local_u;
  int *_local_muscle_id;
  float *_local_fiber;
  float *_local_density;
  float *_local_one_over_h;
  float *_local_c1;
  float *_local_c2;
  float *_local_F_fiber;
  float **_local_activations;
  float **_local_fiber_max_stresses;

public:
    explicit Muscle_Update_Position_Based_State_SCALAR_None (float *u_in,
                                                             int *muscle_id_in,
                                                             float *fiber_in,
                                                             float *density_in,
                                                             float
                                                             *one_over_h_in,
                                                             float *c1_in,
                                                             float *c2_in,
                                                             float *F_fiber_in,
                                                             float
                                                             **activations_in,
                                                             float
                                                             **fiber_max_stresses_in):_local_u
    (u_in), _local_muscle_id (muscle_id_in), _local_fiber (fiber_in),
    _local_density (density_in), _local_one_over_h (one_over_h_in),
    _local_c1 (c1_in), _local_c2 (c2_in), _local_F_fiber (F_fiber_in),
    _local_activations (activations_in),
    _local_fiber_max_stresses (fiber_max_stresses_in)
  {
  }
  void Execute (int index)
  {
    // full array typedefs
    //typedef int (&refArray)[SIZE][3][8];
    typedef float (&fullArray1)[SIZE][3][8][16];
    typedef int (&fullArray2)[SIZE][16];
    typedef float (&fullArray3)[SIZE][3][16];
    typedef float (&fullArray4)[SIZE][16];
    typedef float (&fullArray5)[SIZE][16];
    typedef float (&fullArray6)[SIZE][16];
    typedef float (&fullArray7)[SIZE][16];
    typedef float (&fullArray8)[SIZE][3][16];
    typedef float *(&fullArray9)[SIZE];
    typedef float *(&fullArray10)[SIZE];

    // full array extractions
    //refArray _rA = reinterpret_cast < refArray >(*_A);
    fullArray1 _ru = reinterpret_cast < fullArray1 > (*_local_u);
    fullArray2 _rmuscle_id =
      reinterpret_cast < fullArray2 > (*_local_muscle_id);
    fullArray3 _rfiber = reinterpret_cast < fullArray3 > (*_local_fiber);
    fullArray4 _rdensity = reinterpret_cast < fullArray4 > (*_local_density);
    fullArray5 _rone_over_h =
      reinterpret_cast < fullArray5 > (*_local_one_over_h);
    fullArray6 _rc1 = reinterpret_cast < fullArray6 > (*_local_c1);
    fullArray7 _rc2 = reinterpret_cast < fullArray7 > (*_local_c2);
    fullArray8 _rF_fiber = reinterpret_cast < fullArray8 > (*_local_F_fiber);
    fullArray9 _ractivations =
      reinterpret_cast < fullArray9 > (*_local_activations);
    fullArray10 _rfiber_max_stresses =
      reinterpret_cast < fullArray10 > (*_local_fiber_max_stresses);

    const int ChunkSize = 1;
    // chunk typedef
    //typedef int (&refArray1)[3][8];
    typedef float (&refArray1)[3][8][16];
    typedef int (&refArray2)[16];
    typedef float (&refArray3)[3][16];
    typedef float (&refArray4)[16];
    typedef float (&refArray5)[16];
    typedef float (&refArray6)[16];
    typedef float (&refArray7)[16];
    typedef float (&refArray8)[3][16];
    typedef float *(&refArray9);
    typedef float *(&refArray10);

    for (int chunk_offset = 0; chunk_offset < 16; chunk_offset += ChunkSize)
      {
        //refArray1 Ak = reinterpret_cast < refArray1 > (_rA[entry][0][chunk_offset]);
        refArray1 uk =
          reinterpret_cast < refArray1 > (_ru[index][0][0][chunk_offset]);
        refArray2 muscle_idk =
          reinterpret_cast < refArray2 > (_rmuscle_id[index][chunk_offset]);
        refArray3 fiberk =
          reinterpret_cast < refArray3 > (_rfiber[index][0][chunk_offset]);
        refArray4 densityk =
          reinterpret_cast < refArray4 > (_rdensity[index][chunk_offset]);
        refArray5 one_over_hk =
          reinterpret_cast < refArray5 > (_rone_over_h[index][chunk_offset]);
        refArray6 c1k =
          reinterpret_cast < refArray6 > (_rc1[index][chunk_offset]);
        refArray7 c2k =
          reinterpret_cast < refArray7 > (_rc2[index][chunk_offset]);
        refArray8 F_fiberk =
          reinterpret_cast < refArray8 > (_rF_fiber[index][0][chunk_offset]);
        refArray9 activationsk =
          reinterpret_cast < refArray9 > (_ractivations[index]);
        refArray10 fiber_max_stressesk =
          reinterpret_cast < refArray10 > (_rfiber_max_stresses[index]);

        Muscle_Update_Position_Based_State < float, float[16], int[16] > (uk,
                                                                          muscle_idk,
                                                                          fiberk,
                                                                          densityk,
                                                                          one_over_hk,
                                                                          c1k,
                                                                          c2k,
                                                                          F_fiberk,
                                                                          activationsk,
                                                                          fiber_max_stressesk);
      }

  }
};


#ifdef ENABLE_SSE_INSTRUCTION_SET

template < int SIZE > class Muscle_Update_Position_Based_State_SSE_None
{
private:
  // Generate Variables Here
  float *_local_u;
  int *_local_muscle_id;
  float *_local_fiber;
  float *_local_density;
  float *_local_one_over_h;
  float *_local_c1;
  float *_local_c2;
  float *_local_F_fiber;
  float **_local_activations;
  float **_local_fiber_max_stresses;

public:
    explicit Muscle_Update_Position_Based_State_SSE_None (float *u_in,
                                                          int *muscle_id_in,
                                                          float *fiber_in,
                                                          float *density_in,
                                                          float *one_over_h_in,
                                                          float *c1_in,
                                                          float *c2_in,
                                                          float *F_fiber_in,
                                                          float
                                                          **activations_in,
                                                          float
                                                          **fiber_max_stresses_in):_local_u
    (u_in), _local_muscle_id (muscle_id_in), _local_fiber (fiber_in),
    _local_density (density_in), _local_one_over_h (one_over_h_in),
    _local_c1 (c1_in), _local_c2 (c2_in), _local_F_fiber (F_fiber_in),
    _local_activations (activations_in),
    _local_fiber_max_stresses (fiber_max_stresses_in)
  {
  }
  void Execute (int index)
  {
    // full array typedefs
    //typedef int (&refArray)[SIZE][3][8];
    typedef float (&fullArray1)[SIZE][3][8][16];
    typedef int (&fullArray2)[SIZE][16];
    typedef float (&fullArray3)[SIZE][3][16];
    typedef float (&fullArray4)[SIZE][16];
    typedef float (&fullArray5)[SIZE][16];
    typedef float (&fullArray6)[SIZE][16];
    typedef float (&fullArray7)[SIZE][16];
    typedef float (&fullArray8)[SIZE][3][16];
    typedef float *(&fullArray9)[SIZE];
    typedef float *(&fullArray10)[SIZE];

    // full array extractions
    //refArray _rA = reinterpret_cast < refArray >(*_A);
    fullArray1 _ru = reinterpret_cast < fullArray1 > (*_local_u);
    fullArray2 _rmuscle_id =
      reinterpret_cast < fullArray2 > (*_local_muscle_id);
    fullArray3 _rfiber = reinterpret_cast < fullArray3 > (*_local_fiber);
    fullArray4 _rdensity = reinterpret_cast < fullArray4 > (*_local_density);
    fullArray5 _rone_over_h =
      reinterpret_cast < fullArray5 > (*_local_one_over_h);
    fullArray6 _rc1 = reinterpret_cast < fullArray6 > (*_local_c1);
    fullArray7 _rc2 = reinterpret_cast < fullArray7 > (*_local_c2);
    fullArray8 _rF_fiber = reinterpret_cast < fullArray8 > (*_local_F_fiber);
    fullArray9 _ractivations =
      reinterpret_cast < fullArray9 > (*_local_activations);
    fullArray10 _rfiber_max_stresses =
      reinterpret_cast < fullArray10 > (*_local_fiber_max_stresses);

    const int ChunkSize = 4;
    // chunk typedef
    //typedef int (&refArray1)[3][8];
    typedef float (&refArray1)[3][8][16];
    typedef int (&refArray2)[16];
    typedef float (&refArray3)[3][16];
    typedef float (&refArray4)[16];
    typedef float (&refArray5)[16];
    typedef float (&refArray6)[16];
    typedef float (&refArray7)[16];
    typedef float (&refArray8)[3][16];
    typedef float *(&refArray9);
    typedef float *(&refArray10);

    for (int chunk_offset = 0; chunk_offset < 16; chunk_offset += ChunkSize)
      {
        //refArray1 Ak = reinterpret_cast < refArray1 > (_rA[entry][0][chunk_offset]);
        refArray1 uk =
          reinterpret_cast < refArray1 > (_ru[index][0][0][chunk_offset]);
        refArray2 muscle_idk =
          reinterpret_cast < refArray2 > (_rmuscle_id[index][chunk_offset]);
        refArray3 fiberk =
          reinterpret_cast < refArray3 > (_rfiber[index][0][chunk_offset]);
        refArray4 densityk =
          reinterpret_cast < refArray4 > (_rdensity[index][chunk_offset]);
        refArray5 one_over_hk =
          reinterpret_cast < refArray5 > (_rone_over_h[index][chunk_offset]);
        refArray6 c1k =
          reinterpret_cast < refArray6 > (_rc1[index][chunk_offset]);
        refArray7 c2k =
          reinterpret_cast < refArray7 > (_rc2[index][chunk_offset]);
        refArray8 F_fiberk =
          reinterpret_cast < refArray8 > (_rF_fiber[index][0][chunk_offset]);
        refArray9 activationsk =
          reinterpret_cast < refArray9 > (_ractivations[index]);
        refArray10 fiber_max_stressesk =
          reinterpret_cast < refArray10 > (_rfiber_max_stresses[index]);

        Muscle_Update_Position_Based_State < __m128, float[16], int[16] > (uk,
                                                                           muscle_idk,
                                                                           fiberk,
                                                                           densityk,
                                                                           one_over_hk,
                                                                           c1k,
                                                                           c2k,
                                                                           F_fiberk,
                                                                           activationsk,
                                                                           fiber_max_stressesk);
      }

  }
};

#endif

#ifdef ENABLE_AVX_INSTRUCTION_SET

template < int SIZE > class Muscle_Update_Position_Based_State_AVX_None
{
private:
  // Generate Variables Here
  float *_local_u;
  int *_local_muscle_id;
  float *_local_fiber;
  float *_local_density;
  float *_local_one_over_h;
  float *_local_c1;
  float *_local_c2;
  float *_local_F_fiber;
  float **_local_activations;
  float **_local_fiber_max_stresses;

public:
    explicit Muscle_Update_Position_Based_State_AVX_None (float *u_in,
                                                          int *muscle_id_in,
                                                          float *fiber_in,
                                                          float *density_in,
                                                          float *one_over_h_in,
                                                          float *c1_in,
                                                          float *c2_in,
                                                          float *F_fiber_in,
                                                          float
                                                          **activations_in,
                                                          float
                                                          **fiber_max_stresses_in):_local_u
    (u_in), _local_muscle_id (muscle_id_in), _local_fiber (fiber_in),
    _local_density (density_in), _local_one_over_h (one_over_h_in),
    _local_c1 (c1_in), _local_c2 (c2_in), _local_F_fiber (F_fiber_in),
    _local_activations (activations_in),
    _local_fiber_max_stresses (fiber_max_stresses_in)
  {
  }
  void Execute (int index)
  {
    // full array typedefs
    //typedef int (&refArray)[SIZE][3][8];
    typedef float (&fullArray1)[SIZE][3][8][16];
    typedef int (&fullArray2)[SIZE][16];
    typedef float (&fullArray3)[SIZE][3][16];
    typedef float (&fullArray4)[SIZE][16];
    typedef float (&fullArray5)[SIZE][16];
    typedef float (&fullArray6)[SIZE][16];
    typedef float (&fullArray7)[SIZE][16];
    typedef float (&fullArray8)[SIZE][3][16];
    typedef float *(&fullArray9)[SIZE];
    typedef float *(&fullArray10)[SIZE];

    // full array extractions
    //refArray _rA = reinterpret_cast < refArray >(*_A);
    fullArray1 _ru = reinterpret_cast < fullArray1 > (*_local_u);
    fullArray2 _rmuscle_id =
      reinterpret_cast < fullArray2 > (*_local_muscle_id);
    fullArray3 _rfiber = reinterpret_cast < fullArray3 > (*_local_fiber);
    fullArray4 _rdensity = reinterpret_cast < fullArray4 > (*_local_density);
    fullArray5 _rone_over_h =
      reinterpret_cast < fullArray5 > (*_local_one_over_h);
    fullArray6 _rc1 = reinterpret_cast < fullArray6 > (*_local_c1);
    fullArray7 _rc2 = reinterpret_cast < fullArray7 > (*_local_c2);
    fullArray8 _rF_fiber = reinterpret_cast < fullArray8 > (*_local_F_fiber);
    fullArray9 _ractivations =
      reinterpret_cast < fullArray9 > (*_local_activations);
    fullArray10 _rfiber_max_stresses =
      reinterpret_cast < fullArray10 > (*_local_fiber_max_stresses);

    const int ChunkSize = 8;
    // chunk typedef
    //typedef int (&refArray1)[3][8];
    typedef float (&refArray1)[3][8][16];
    typedef int (&refArray2)[16];
    typedef float (&refArray3)[3][16];
    typedef float (&refArray4)[16];
    typedef float (&refArray5)[16];
    typedef float (&refArray6)[16];
    typedef float (&refArray7)[16];
    typedef float (&refArray8)[3][16];
    typedef float *(&refArray9);
    typedef float *(&refArray10);

    for (int chunk_offset = 0; chunk_offset < 16; chunk_offset += ChunkSize)
      {
        //refArray1 Ak = reinterpret_cast < refArray1 > (_rA[entry][0][chunk_offset]);
        refArray1 uk =
          reinterpret_cast < refArray1 > (_ru[index][0][0][chunk_offset]);
        refArray2 muscle_idk =
          reinterpret_cast < refArray2 > (_rmuscle_id[index][chunk_offset]);
        refArray3 fiberk =
          reinterpret_cast < refArray3 > (_rfiber[index][0][chunk_offset]);
        refArray4 densityk =
          reinterpret_cast < refArray4 > (_rdensity[index][chunk_offset]);
        refArray5 one_over_hk =
          reinterpret_cast < refArray5 > (_rone_over_h[index][chunk_offset]);
        refArray6 c1k =
          reinterpret_cast < refArray6 > (_rc1[index][chunk_offset]);
        refArray7 c2k =
          reinterpret_cast < refArray7 > (_rc2[index][chunk_offset]);
        refArray8 F_fiberk =
          reinterpret_cast < refArray8 > (_rF_fiber[index][0][chunk_offset]);
        refArray9 activationsk =
          reinterpret_cast < refArray9 > (_ractivations[index]);
        refArray10 fiber_max_stressesk =
          reinterpret_cast < refArray10 > (_rfiber_max_stresses[index]);

        Muscle_Update_Position_Based_State < __m256, float[16], int[16] > (uk,
                                                                           muscle_idk,
                                                                           fiberk,
                                                                           densityk,
                                                                           one_over_hk,
                                                                           c1k,
                                                                           c2k,
                                                                           F_fiberk,
                                                                           activationsk,
                                                                           fiber_max_stressesk);
      }

  }
};

#endif

#ifdef ENABLE_NEON_INSTRUCTION_SET

template < int SIZE > class Muscle_Update_Position_Based_State_NEON_None
{
private:
  // Generate Variables Here
  float *_local_u;
  int *_local_muscle_id;
  float *_local_fiber;
  float *_local_density;
  float *_local_one_over_h;
  float *_local_c1;
  float *_local_c2;
  float *_local_F_fiber;
  float **_local_activations;
  float **_local_fiber_max_stresses;

public:
    explicit Muscle_Update_Position_Based_State_NEON_None (float *u_in,
                                                           int *muscle_id_in,
                                                           float *fiber_in,
                                                           float *density_in,
                                                           float *one_over_h_in,
                                                           float *c1_in,
                                                           float *c2_in,
                                                           float *F_fiber_in,
                                                           float
                                                           **activations_in,
                                                           float
                                                           **fiber_max_stresses_in):_local_u
    (u_in), _local_muscle_id (muscle_id_in), _local_fiber (fiber_in),
    _local_density (density_in), _local_one_over_h (one_over_h_in),
    _local_c1 (c1_in), _local_c2 (c2_in), _local_F_fiber (F_fiber_in),
    _local_activations (activations_in),
    _local_fiber_max_stresses (fiber_max_stresses_in)
  {
  }
  void Execute (int index)
  {
    // full array typedefs
    //typedef int (&refArray)[SIZE][3][8];
    typedef float (&fullArray1)[SIZE][3][8][16];
    typedef int (&fullArray2)[SIZE][16];
    typedef float (&fullArray3)[SIZE][3][16];
    typedef float (&fullArray4)[SIZE][16];
    typedef float (&fullArray5)[SIZE][16];
    typedef float (&fullArray6)[SIZE][16];
    typedef float (&fullArray7)[SIZE][16];
    typedef float (&fullArray8)[SIZE][3][16];
    typedef float *(&fullArray9)[SIZE];
    typedef float *(&fullArray10)[SIZE];

    // full array extractions
    //refArray _rA = reinterpret_cast < refArray >(*_A);
    fullArray1 _ru = reinterpret_cast < fullArray1 > (*_local_u);
    fullArray2 _rmuscle_id =
      reinterpret_cast < fullArray2 > (*_local_muscle_id);
    fullArray3 _rfiber = reinterpret_cast < fullArray3 > (*_local_fiber);
    fullArray4 _rdensity = reinterpret_cast < fullArray4 > (*_local_density);
    fullArray5 _rone_over_h =
      reinterpret_cast < fullArray5 > (*_local_one_over_h);
    fullArray6 _rc1 = reinterpret_cast < fullArray6 > (*_local_c1);
    fullArray7 _rc2 = reinterpret_cast < fullArray7 > (*_local_c2);
    fullArray8 _rF_fiber = reinterpret_cast < fullArray8 > (*_local_F_fiber);
    fullArray9 _ractivations =
      reinterpret_cast < fullArray9 > (*_local_activations);
    fullArray10 _rfiber_max_stresses =
      reinterpret_cast < fullArray10 > (*_local_fiber_max_stresses);

    const int ChunkSize = 4;
    // chunk typedef
    //typedef int (&refArray1)[3][8];
    typedef float (&refArray1)[3][8][16];
    typedef int (&refArray2)[16];
    typedef float (&refArray3)[3][16];
    typedef float (&refArray4)[16];
    typedef float (&refArray5)[16];
    typedef float (&refArray6)[16];
    typedef float (&refArray7)[16];
    typedef float (&refArray8)[3][16];
    typedef float *(&refArray9);
    typedef float *(&refArray10);

    for (int chunk_offset = 0; chunk_offset < 16; chunk_offset += ChunkSize)
      {
        //refArray1 Ak = reinterpret_cast < refArray1 > (_rA[entry][0][chunk_offset]);
        refArray1 uk =
          reinterpret_cast < refArray1 > (_ru[index][0][0][chunk_offset]);
        refArray2 muscle_idk =
          reinterpret_cast < refArray2 > (_rmuscle_id[index][chunk_offset]);
        refArray3 fiberk =
          reinterpret_cast < refArray3 > (_rfiber[index][0][chunk_offset]);
        refArray4 densityk =
          reinterpret_cast < refArray4 > (_rdensity[index][chunk_offset]);
        refArray5 one_over_hk =
          reinterpret_cast < refArray5 > (_rone_over_h[index][chunk_offset]);
        refArray6 c1k =
          reinterpret_cast < refArray6 > (_rc1[index][chunk_offset]);
        refArray7 c2k =
          reinterpret_cast < refArray7 > (_rc2[index][chunk_offset]);
        refArray8 F_fiberk =
          reinterpret_cast < refArray8 > (_rF_fiber[index][0][chunk_offset]);
        refArray9 activationsk =
          reinterpret_cast < refArray9 > (_ractivations[index]);
        refArray10 fiber_max_stressesk =
          reinterpret_cast < refArray10 > (_rfiber_max_stresses[index]);

        Muscle_Update_Position_Based_State < float32x4_t, float[16],
          int[16] > (uk, muscle_idk, fiberk, densityk, one_over_hk, c1k, c2k,
                     F_fiberk, activationsk, fiber_max_stressesk);
      }

  }
};

#endif

#ifdef ENABLE_MIC_INSTRUCTION_SET

template < int SIZE > class Muscle_Update_Position_Based_State_MIC_None
{
private:
  // Generate Variables Here
  float *_local_u;
  int *_local_muscle_id;
  float *_local_fiber;
  float *_local_density;
  float *_local_one_over_h;
  float *_local_c1;
  float *_local_c2;
  float *_local_F_fiber;
  float **_local_activations;
  float **_local_fiber_max_stresses;

public:
    explicit Muscle_Update_Position_Based_State_MIC_None (float *u_in,
                                                          int *muscle_id_in,
                                                          float *fiber_in,
                                                          float *density_in,
                                                          float *one_over_h_in,
                                                          float *c1_in,
                                                          float *c2_in,
                                                          float *F_fiber_in,
                                                          float
                                                          **activations_in,
                                                          float
                                                          **fiber_max_stresses_in):_local_u
    (u_in), _local_muscle_id (muscle_id_in), _local_fiber (fiber_in),
    _local_density (density_in), _local_one_over_h (one_over_h_in),
    _local_c1 (c1_in), _local_c2 (c2_in), _local_F_fiber (F_fiber_in),
    _local_activations (activations_in),
    _local_fiber_max_stresses (fiber_max_stresses_in)
  {
  }
  void Execute (int index)
  {
    // full array typedefs
    //typedef int (&refArray)[SIZE][3][8];
    typedef float (&fullArray1)[SIZE][3][8][16];
    typedef int (&fullArray2)[SIZE][16];
    typedef float (&fullArray3)[SIZE][3][16];
    typedef float (&fullArray4)[SIZE][16];
    typedef float (&fullArray5)[SIZE][16];
    typedef float (&fullArray6)[SIZE][16];
    typedef float (&fullArray7)[SIZE][16];
    typedef float (&fullArray8)[SIZE][3][16];
    typedef float *(&fullArray9)[SIZE];
    typedef float *(&fullArray10)[SIZE];

    // full array extractions
    //refArray _rA = reinterpret_cast < refArray >(*_A);
    fullArray1 _ru = reinterpret_cast < fullArray1 > (*_local_u);
    fullArray2 _rmuscle_id =
      reinterpret_cast < fullArray2 > (*_local_muscle_id);
    fullArray3 _rfiber = reinterpret_cast < fullArray3 > (*_local_fiber);
    fullArray4 _rdensity = reinterpret_cast < fullArray4 > (*_local_density);
    fullArray5 _rone_over_h =
      reinterpret_cast < fullArray5 > (*_local_one_over_h);
    fullArray6 _rc1 = reinterpret_cast < fullArray6 > (*_local_c1);
    fullArray7 _rc2 = reinterpret_cast < fullArray7 > (*_local_c2);
    fullArray8 _rF_fiber = reinterpret_cast < fullArray8 > (*_local_F_fiber);
    fullArray9 _ractivations =
      reinterpret_cast < fullArray9 > (*_local_activations);
    fullArray10 _rfiber_max_stresses =
      reinterpret_cast < fullArray10 > (*_local_fiber_max_stresses);

    const int ChunkSize = 16;
    // chunk typedef
    //typedef int (&refArray1)[3][8];
    typedef float (&refArray1)[3][8][16];
    typedef int (&refArray2)[16];
    typedef float (&refArray3)[3][16];
    typedef float (&refArray4)[16];
    typedef float (&refArray5)[16];
    typedef float (&refArray6)[16];
    typedef float (&refArray7)[16];
    typedef float (&refArray8)[3][16];
    typedef float *(&refArray9);
    typedef float *(&refArray10);

    for (int chunk_offset = 0; chunk_offset < 16; chunk_offset += ChunkSize)
      {
        //refArray1 Ak = reinterpret_cast < refArray1 > (_rA[entry][0][chunk_offset]);
        refArray1 uk =
          reinterpret_cast < refArray1 > (_ru[index][0][0][chunk_offset]);
        refArray2 muscle_idk =
          reinterpret_cast < refArray2 > (_rmuscle_id[index][chunk_offset]);
        refArray3 fiberk =
          reinterpret_cast < refArray3 > (_rfiber[index][0][chunk_offset]);
        refArray4 densityk =
          reinterpret_cast < refArray4 > (_rdensity[index][chunk_offset]);
        refArray5 one_over_hk =
          reinterpret_cast < refArray5 > (_rone_over_h[index][chunk_offset]);
        refArray6 c1k =
          reinterpret_cast < refArray6 > (_rc1[index][chunk_offset]);
        refArray7 c2k =
          reinterpret_cast < refArray7 > (_rc2[index][chunk_offset]);
        refArray8 F_fiberk =
          reinterpret_cast < refArray8 > (_rF_fiber[index][0][chunk_offset]);
        refArray9 activationsk =
          reinterpret_cast < refArray9 > (_ractivations[index]);
        refArray10 fiber_max_stressesk =
          reinterpret_cast < refArray10 > (_rfiber_max_stresses[index]);

        Muscle_Update_Position_Based_State < __m512, float[16], int[16] > (uk,
                                                                           muscle_idk,
                                                                           fiberk,
                                                                           densityk,
                                                                           one_over_hk,
                                                                           c1k,
                                                                           c2k,
                                                                           F_fiberk,
                                                                           activationsk,
                                                                           fiber_max_stressesk);
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
      cout << "Running Thread Test for Muscle_Update_Position_Based_State " <<
      std::endl;

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
    typedef int (&muscle_id_type)[data_size][16];
    muscle_id_type muscle_id =
      reinterpret_cast < muscle_id_type >
      (*((T *) (_mm_malloc (data_size * 16 * sizeof (int), 64))));
    typedef T (&fiber_type)[data_size][3][16];
    fiber_type fiber =
      reinterpret_cast < fiber_type >
      (*((T *) (_mm_malloc (data_size * 3 * 16 * sizeof (T), 64))));
    typedef T (&density_type)[data_size][16];
    density_type density =
      reinterpret_cast < density_type >
      (*((T *) (_mm_malloc (data_size * 16 * sizeof (T), 64))));
    typedef T (&one_over_h_type)[data_size][16];
    one_over_h_type one_over_h =
      reinterpret_cast < one_over_h_type >
      (*((T *) (_mm_malloc (data_size * 16 * sizeof (T), 64))));
    typedef T (&c1_type)[data_size][16];
    c1_type c1 =
      reinterpret_cast < c1_type >
      (*((T *) (_mm_malloc (data_size * 16 * sizeof (T), 64))));
    typedef T (&c2_type)[data_size][16];
    c2_type c2 =
      reinterpret_cast < c2_type >
      (*((T *) (_mm_malloc (data_size * 16 * sizeof (T), 64))));
    typedef T (&F_fiber_type)[data_size][3][16];
    F_fiber_type F_fiber =
      reinterpret_cast < F_fiber_type >
      (*((T *) (_mm_malloc (data_size * 3 * 16 * sizeof (T), 64))));
    typedef float *(&activations_type)[data_size];
    activations_type activations =
      reinterpret_cast < activations_type >
      (*((T *) (_mm_malloc (data_size * sizeof (float *), 64))));
    typedef float *(&fiber_max_stresses_type)[data_size];
    fiber_max_stresses_type fiber_max_stresses =
      reinterpret_cast < fiber_max_stresses_type >
      (*((T *) (_mm_malloc (data_size * sizeof (float *), 64))));


    for (int __a = 0; __a < data_size; __a++)
      for (int __b = 0; __b < 3; __b++)
        for (int __c = 0; __c < 8; __c++)
          for (int __d = 0; __d < 16; __d++)
            {
              u[__a][__b][__c][__d] = Get_Random < float >();
    } for (int __a = 0; __a < data_size; __a++)
      for (int __b = 0; __b < 16; __b++)
        {
          muscle_id[__a][__b] = Get_Random < int >(1, 99);
    } for (int __a = 0; __a < data_size; __a++)
      for (int __b = 0; __b < 3; __b++)
        for (int __c = 0; __c < 16; __c++)
          {
            fiber[__a][__b][__c] = Get_Random < float >();
    } for (int __a = 0; __a < data_size; __a++)
      for (int __b = 0; __b < 16; __b++)
        {
          density[__a][__b] = Get_Random < float >();
    } for (int __a = 0; __a < data_size; __a++)
      for (int __b = 0; __b < 16; __b++)
        {
          one_over_h[__a][__b] = Get_Random < float >();
    } for (int __a = 0; __a < data_size; __a++)
      for (int __b = 0; __b < 16; __b++)
        {
          c1[__a][__b] = Get_Random < float >();
    } for (int __a = 0; __a < data_size; __a++)
      for (int __b = 0; __b < 16; __b++)
        {
          c2[__a][__b] = Get_Random < float >();
    } for (int __a = 0; __a < data_size; __a++)
      for (int __b = 0; __b < 3; __b++)
        for (int __c = 0; __c < 16; __c++)
          {
            F_fiber[__a][__b][__c] = Get_Random < float >();
    } for (int __a = 0; __a < data_size; __a++)
      {
        activations[__a] = new float[100];
        for (int __x__ = 0; __x__ < 100; __x__++)
          activations[__a][__x__] = Get_Random < float >();;
    } for (int __a = 0; __a < data_size; __a++)
      {
        fiber_max_stresses[__a] = new float[100];
        for (int __x__ = 0; __x__ < 100; __x__++)
          fiber_max_stresses[__a][__x__] = Get_Random < float >();;
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


      Muscle_Update_Position_Based_State_SCALAR_None < data_size >
        op ((float *) &u, (int *) &muscle_id, (float *) &fiber,
            (float *) &density, (float *) &one_over_h, (float *) &c1,
            (float *) &c2, (float *) &F_fiber, (float **) &activations,
            (float **) &fiber_max_stresses);

      for (int t = threads; t <= threads_max;
           t += std::max < int >(((threads_max - threads) / 30), 1))
        {
          std::cout << "Running Test with " << t << " threads." << std::endl;
          MT_Streaming_Kernels::Kernel_Serial_Base_Helper <
            Muscle_Update_Position_Based_State_SCALAR_None < data_size >
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


      Muscle_Update_Position_Based_State_SSE_None < data_size >
        op ((float *) &u, (int *) &muscle_id, (float *) &fiber,
            (float *) &density, (float *) &one_over_h, (float *) &c1,
            (float *) &c2, (float *) &F_fiber, (float **) &activations,
            (float **) &fiber_max_stresses);

      for (int t = threads; t <= threads_max;
           t += std::max < int >(((threads_max - threads) / 30), 1))
        {
          std::cout << "Running Test with " << t << " threads." << std::endl;
          MT_Streaming_Kernels::Kernel_Serial_Base_Helper <
            Muscle_Update_Position_Based_State_SSE_None < data_size >
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


      Muscle_Update_Position_Based_State_AVX_None < data_size >
        op ((float *) &u, (int *) &muscle_id, (float *) &fiber,
            (float *) &density, (float *) &one_over_h, (float *) &c1,
            (float *) &c2, (float *) &F_fiber, (float **) &activations,
            (float **) &fiber_max_stresses);

      for (int t = threads; t <= threads_max;
           t += std::max < int >(((threads_max - threads) / 30), 1))
        {
          std::cout << "Running Test with " << t << " threads." << std::endl;
          MT_Streaming_Kernels::Kernel_Serial_Base_Helper <
            Muscle_Update_Position_Based_State_AVX_None < data_size >
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


      Muscle_Update_Position_Based_State_NEON_None < data_size >
        op ((float *) &u, (int *) &muscle_id, (float *) &fiber,
            (float *) &density, (float *) &one_over_h, (float *) &c1,
            (float *) &c2, (float *) &F_fiber, (float **) &activations,
            (float **) &fiber_max_stresses);

      for (int t = threads; t <= threads_max;
           t += std::max < int >(((threads_max - threads) / 30), 1))
        {
          std::cout << "Running Test with " << t << " threads." << std::endl;
          MT_Streaming_Kernels::Kernel_Serial_Base_Helper <
            Muscle_Update_Position_Based_State_NEON_None < data_size >
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


      Muscle_Update_Position_Based_State_MIC_None < data_size >
        op ((float *) &u, (int *) &muscle_id, (float *) &fiber,
            (float *) &density, (float *) &one_over_h, (float *) &c1,
            (float *) &c2, (float *) &F_fiber, (float **) &activations,
            (float **) &fiber_max_stresses);

      for (int t = threads; t <= threads_max;
           t += std::max < int >(((threads_max - threads) / 30), 1))
        {
          std::cout << "Running Test with " << t << " threads." << std::endl;
          MT_Streaming_Kernels::Kernel_Serial_Base_Helper <
            Muscle_Update_Position_Based_State_MIC_None < data_size >
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

    _mm_free (reinterpret_cast < void *>(u));
    _mm_free (reinterpret_cast < void *>(muscle_id));
    _mm_free (reinterpret_cast < void *>(fiber));
    _mm_free (reinterpret_cast < void *>(density));
    _mm_free (reinterpret_cast < void *>(one_over_h));
    _mm_free (reinterpret_cast < void *>(c1));
    _mm_free (reinterpret_cast < void *>(c2));
    _mm_free (reinterpret_cast < void *>(F_fiber));
    _mm_free (reinterpret_cast < void *>(activations));
    _mm_free (reinterpret_cast < void *>(fiber_max_stresses));


  }


  return 0;

}

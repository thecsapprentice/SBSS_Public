
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <sys/time.h>
#include "KernelCommon.h"
struct NEOHOOKEAN_TAG;
struct COROTATED_TAG;

#include "Muscle_Update_Position_Based_State.h"

#define NUM_TRIALS 1000000

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

int
main (int argc, char *argv[])
{
  typedef float T;

  std::cout << "Preparing to Run " << NUM_TRIALS << " of all kernels." << std::
    endl;

  int seed = 1;
  if (argc == 2)
    seed = atoi (argv[1]);
  srand (seed);



  {
    std::
      cout << "Running Stream Test for Muscle_Update_Position_Based_State " <<
      std::endl;

//=======================================================
//
//        DEFINE ALL VARIABLES USED BY KERNEL
//
//=======================================================

    T u[3][8][16] __attribute__ ((aligned (64)));
    int muscle_id[16] __attribute__ ((aligned (64)));
    T fiber[3][16] __attribute__ ((aligned (64)));
    T density[16] __attribute__ ((aligned (64)));
    T one_over_h[16] __attribute__ ((aligned (64)));
    T c1[16] __attribute__ ((aligned (64)));
    T c1_reference[16] __attribute__ ((aligned (64)));
    T c1_original[16] __attribute__ ((aligned (64)));
    T c2[16] __attribute__ ((aligned (64)));
    T c2_reference[16] __attribute__ ((aligned (64)));
    T c2_original[16] __attribute__ ((aligned (64)));
    T F_fiber[3][16] __attribute__ ((aligned (64)));
    T F_fiber_reference[3][16] __attribute__ ((aligned (64)));
    T F_fiber_original[3][16] __attribute__ ((aligned (64)));
    float *activations __attribute__ ((aligned (64)));
    float *fiber_max_stresses __attribute__ ((aligned (64)));


    for (int __a = 0; __a < 3; __a++)
      for (int __b = 0; __b < 8; __b++)
        for (int __c = 0; __c < 16; __c++)
          u[__a][__b][__c] = Get_Random < float >();
    for (int __a = 0; __a < 16; __a++)
      muscle_id[__a] = Get_Random < int >(1, 99);
    for (int __a = 0; __a < 3; __a++)
      for (int __b = 0; __b < 16; __b++)
        fiber[__a][__b] = Get_Random < float >();
    for (int __a = 0; __a < 16; __a++)
      density[__a] = Get_Random < float >();
    for (int __a = 0; __a < 16; __a++)
      one_over_h[__a] = Get_Random < float >();
    for (int __a = 0; __a < 16; __a++)
      {
        c1_original[__a] = Get_Random < float >();
        c1[__a] = c1_original[__a];
        c1_reference[__a] = c1_original[__a];
      }
    for (int __a = 0; __a < 16; __a++)
      {
        c2_original[__a] = Get_Random < float >();
        c2[__a] = c2_original[__a];
        c2_reference[__a] = c2_original[__a];
      }
    for (int __a = 0; __a < 3; __a++)
      for (int __b = 0; __b < 16; __b++)
        {
          F_fiber_original[__a][__b] = Get_Random < float >();
          F_fiber[__a][__b] = F_fiber_original[__a][__b];
          F_fiber_reference[__a][__b] = F_fiber_original[__a][__b];
        }

    activations = new float[100];
    for (int __x__ = 0; __x__ < 100; __x__++)
      activations[__x__] = Get_Random < float >();;
    fiber_max_stresses = new float[100];
    for (int __x__ = 0; __x__ < 100; __x__++)
      fiber_max_stresses[__x__] = Get_Random < float >();;

//=======================================================
//
//             COMPUTE SCALAR RESULTS
//
//=======================================================

    {
      std::cout << "	Running " << NUM_TRIALS << " of SCALAR :  ";
      start_timer ();
      for (int n = 0; n < NUM_TRIALS; n++)
        {
          typedef T (&refArray1)[3][8][16];
          typedef int (&refArray2)[16];
          typedef T (&refArray3)[3][16];
          typedef T (&refArray4)[16];
          typedef T (&refArray5)[16];
          typedef T (&refArray6)[16];
          typedef T (&refArray7)[16];
          typedef T (&refArray8)[3][16];
          typedef float *(&refArray9);
          typedef float *(&refArray10);
          for (int i = 0; i < 16; i += 1)
            {
              refArray1 uk = reinterpret_cast < refArray1 > (u[0][0][i]);
              refArray2 muscle_idk =
                reinterpret_cast < refArray2 > (muscle_id[i]);
              refArray3 fiberk = reinterpret_cast < refArray3 > (fiber[0][i]);
              refArray4 densityk = reinterpret_cast < refArray4 > (density[i]);
              refArray5 one_over_hk =
                reinterpret_cast < refArray5 > (one_over_h[i]);
              refArray6 c1k = reinterpret_cast < refArray6 > (c1[i]);
              refArray7 c2k = reinterpret_cast < refArray7 > (c2[i]);
              refArray8 F_fiberk =
                reinterpret_cast < refArray8 > (F_fiber[0][i]);
              refArray9 activationsk =
                reinterpret_cast < refArray9 > (activations);
              refArray10 fiber_max_stressesk =
                reinterpret_cast < refArray10 > (fiber_max_stresses);
              Muscle_Update_Position_Based_State < float, float[16],
                int[16] > (uk, muscle_idk, fiberk, densityk, one_over_hk, c1k,
                           c2k, F_fiberk, activationsk, fiber_max_stressesk);
            }

        }
      stop_timer ();
      std::cout << get_time () << "s" << std::endl;
    }

//=======================================================
//
//             COMPUTE SSE RESULTS
//
//=======================================================

#ifdef ENABLE_SSE_INSTRUCTION_SET
    {
      std::cout << "	Running " << NUM_TRIALS << " of SSE :  ";
      start_timer ();
      for (int n = 0; n < NUM_TRIALS; n++)
        {
          typedef T (&refArray1)[3][8][16];
          typedef int (&refArray2)[16];
          typedef T (&refArray3)[3][16];
          typedef T (&refArray4)[16];
          typedef T (&refArray5)[16];
          typedef T (&refArray6)[16];
          typedef T (&refArray7)[16];
          typedef T (&refArray8)[3][16];
          typedef float *(&refArray9);
          typedef float *(&refArray10);
          for (int i = 0; i < 16; i += 4)
            {
              refArray1 uk = reinterpret_cast < refArray1 > (u[0][0][i]);
              refArray2 muscle_idk =
                reinterpret_cast < refArray2 > (muscle_id[i]);
              refArray3 fiberk = reinterpret_cast < refArray3 > (fiber[0][i]);
              refArray4 densityk = reinterpret_cast < refArray4 > (density[i]);
              refArray5 one_over_hk =
                reinterpret_cast < refArray5 > (one_over_h[i]);
              refArray6 c1k = reinterpret_cast < refArray6 > (c1[i]);
              refArray7 c2k = reinterpret_cast < refArray7 > (c2[i]);
              refArray8 F_fiberk =
                reinterpret_cast < refArray8 > (F_fiber[0][i]);
              refArray9 activationsk =
                reinterpret_cast < refArray9 > (activations);
              refArray10 fiber_max_stressesk =
                reinterpret_cast < refArray10 > (fiber_max_stresses);
              Muscle_Update_Position_Based_State < __m128, float[16],
                int[16] > (uk, muscle_idk, fiberk, densityk, one_over_hk, c1k,
                           c2k, F_fiberk, activationsk, fiber_max_stressesk);
            }

        }
      stop_timer ();
      std::cout << get_time () << "s" << std::endl;
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
      start_timer ();
      for (int n = 0; n < NUM_TRIALS; n++)
        {
          typedef T (&refArray1)[3][8][16];
          typedef int (&refArray2)[16];
          typedef T (&refArray3)[3][16];
          typedef T (&refArray4)[16];
          typedef T (&refArray5)[16];
          typedef T (&refArray6)[16];
          typedef T (&refArray7)[16];
          typedef T (&refArray8)[3][16];
          typedef float *(&refArray9);
          typedef float *(&refArray10);
          for (int i = 0; i < 16; i += 8)
            {
              refArray1 uk = reinterpret_cast < refArray1 > (u[0][0][i]);
              refArray2 muscle_idk =
                reinterpret_cast < refArray2 > (muscle_id[i]);
              refArray3 fiberk = reinterpret_cast < refArray3 > (fiber[0][i]);
              refArray4 densityk = reinterpret_cast < refArray4 > (density[i]);
              refArray5 one_over_hk =
                reinterpret_cast < refArray5 > (one_over_h[i]);
              refArray6 c1k = reinterpret_cast < refArray6 > (c1[i]);
              refArray7 c2k = reinterpret_cast < refArray7 > (c2[i]);
              refArray8 F_fiberk =
                reinterpret_cast < refArray8 > (F_fiber[0][i]);
              refArray9 activationsk =
                reinterpret_cast < refArray9 > (activations);
              refArray10 fiber_max_stressesk =
                reinterpret_cast < refArray10 > (fiber_max_stresses);
              Muscle_Update_Position_Based_State < __m256, float[16],
                int[16] > (uk, muscle_idk, fiberk, densityk, one_over_hk, c1k,
                           c2k, F_fiberk, activationsk, fiber_max_stressesk);
            }

        }
      stop_timer ();
      std::cout << get_time () << "s" << std::endl;
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
      start_timer ();
      for (int n = 0; n < NUM_TRIALS; n++)
        {
          typedef T (&refArray1)[3][8][16];
          typedef int (&refArray2)[16];
          typedef T (&refArray3)[3][16];
          typedef T (&refArray4)[16];
          typedef T (&refArray5)[16];
          typedef T (&refArray6)[16];
          typedef T (&refArray7)[16];
          typedef T (&refArray8)[3][16];
          typedef float *(&refArray9);
          typedef float *(&refArray10);
          for (int i = 0; i < 16; i += 4)
            {
              refArray1 uk = reinterpret_cast < refArray1 > (u[0][0][i]);
              refArray2 muscle_idk =
                reinterpret_cast < refArray2 > (muscle_id[i]);
              refArray3 fiberk = reinterpret_cast < refArray3 > (fiber[0][i]);
              refArray4 densityk = reinterpret_cast < refArray4 > (density[i]);
              refArray5 one_over_hk =
                reinterpret_cast < refArray5 > (one_over_h[i]);
              refArray6 c1k = reinterpret_cast < refArray6 > (c1[i]);
              refArray7 c2k = reinterpret_cast < refArray7 > (c2[i]);
              refArray8 F_fiberk =
                reinterpret_cast < refArray8 > (F_fiber[0][i]);
              refArray9 activationsk =
                reinterpret_cast < refArray9 > (activations);
              refArray10 fiber_max_stressesk =
                reinterpret_cast < refArray10 > (fiber_max_stresses);
              Muscle_Update_Position_Based_State < float32x4_t, float[16],
                int[16] > (uk, muscle_idk, fiberk, densityk, one_over_hk, c1k,
                           c2k, F_fiberk, activationsk, fiber_max_stressesk);
            }

        }
      stop_timer ();
      std::cout << get_time () << "s" << std::endl;
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
      start_timer ();
      for (int n = 0; n < NUM_TRIALS; n++)
        {
          typedef T (&refArray1)[3][8][16];
          typedef int (&refArray2)[16];
          typedef T (&refArray3)[3][16];
          typedef T (&refArray4)[16];
          typedef T (&refArray5)[16];
          typedef T (&refArray6)[16];
          typedef T (&refArray7)[16];
          typedef T (&refArray8)[3][16];
          typedef float *(&refArray9);
          typedef float *(&refArray10);
          for (int i = 0; i < 16; i += 16)
            {
              refArray1 uk = reinterpret_cast < refArray1 > (u[0][0][i]);
              refArray2 muscle_idk =
                reinterpret_cast < refArray2 > (muscle_id[i]);
              refArray3 fiberk = reinterpret_cast < refArray3 > (fiber[0][i]);
              refArray4 densityk = reinterpret_cast < refArray4 > (density[i]);
              refArray5 one_over_hk =
                reinterpret_cast < refArray5 > (one_over_h[i]);
              refArray6 c1k = reinterpret_cast < refArray6 > (c1[i]);
              refArray7 c2k = reinterpret_cast < refArray7 > (c2[i]);
              refArray8 F_fiberk =
                reinterpret_cast < refArray8 > (F_fiber[0][i]);
              refArray9 activationsk =
                reinterpret_cast < refArray9 > (activations);
              refArray10 fiber_max_stressesk =
                reinterpret_cast < refArray10 > (fiber_max_stresses);
              Muscle_Update_Position_Based_State < __m512, float[16],
                int[16] > (uk, muscle_idk, fiberk, densityk, one_over_hk, c1k,
                           c2k, F_fiberk, activationsk, fiber_max_stressesk);
            }

        }
      stop_timer ();
      std::cout << get_time () << "s" << std::endl;
    }
#endif

  }



  return 0;

}

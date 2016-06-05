
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <sys/time.h>
#include "KernelCommon.h"
struct NEOHOOKEAN_TAG;
struct COROTATED_TAG;

#include "Compute_Diagonal_Contribution.h"

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
      cout << "Running Stream Test for Compute_Diagonal_Contribution " << std::
      endl;

//=======================================================
//
//        DEFINE ALL VARIABLES USED BY KERNEL
//
//=======================================================

    T one_over_h[16] __attribute__ ((aligned (64)));
    T mu_stab[16] __attribute__ ((aligned (64)));
    T cell_volume[16] __attribute__ ((aligned (64)));
    T U[9][16] __attribute__ ((aligned (64)));
    T V[9][16] __attribute__ ((aligned (64)));
    T dPdF[12][16] __attribute__ ((aligned (64)));
    T d[3][8][16] __attribute__ ((aligned (64)));
    T d_reference[3][8][16] __attribute__ ((aligned (64)));
    T d_original[3][8][16] __attribute__ ((aligned (64)));


    for (int __a = 0; __a < 16; __a++)
      one_over_h[__a] = Get_Random < float >();
    for (int __a = 0; __a < 16; __a++)
      mu_stab[__a] = Get_Random < float >();
    for (int __a = 0; __a < 16; __a++)
      cell_volume[__a] = Get_Random < float >();
    for (int __a = 0; __a < 9; __a++)
      for (int __b = 0; __b < 16; __b++)
        U[__a][__b] = Get_Random < float >();
    for (int __a = 0; __a < 9; __a++)
      for (int __b = 0; __b < 16; __b++)
        V[__a][__b] = Get_Random < float >();
    for (int __a = 0; __a < 12; __a++)
      for (int __b = 0; __b < 16; __b++)
        dPdF[__a][__b] = Get_Random < float >();
    for (int __a = 0; __a < 3; __a++)
      for (int __b = 0; __b < 8; __b++)
        for (int __c = 0; __c < 16; __c++)
          {
            d_original[__a][__b][__c] = Get_Random < float >();
            d[__a][__b][__c] = d_original[__a][__b][__c];
            d_reference[__a][__b][__c] = d_original[__a][__b][__c];
          }


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
          typedef T (&refArray1)[16];
          typedef T (&refArray2)[16];
          typedef T (&refArray3)[16];
          typedef T (&refArray4)[9][16];
          typedef T (&refArray5)[9][16];
          typedef T (&refArray6)[12][16];
          typedef T (&refArray7)[3][8][16];
          for (int i = 0; i < 16; i += 1)
            {
              refArray1 one_over_hk =
                reinterpret_cast < refArray1 > (one_over_h[i]);
              refArray2 mu_stabk = reinterpret_cast < refArray2 > (mu_stab[i]);
              refArray3 cell_volumek =
                reinterpret_cast < refArray3 > (cell_volume[i]);
              refArray4 Uk = reinterpret_cast < refArray4 > (U[0][i]);
              refArray5 Vk = reinterpret_cast < refArray5 > (V[0][i]);
              refArray6 dPdFk = reinterpret_cast < refArray6 > (dPdF[0][i]);
              refArray7 dk = reinterpret_cast < refArray7 > (d[0][0][i]);
              Compute_Diagonal_Contribution < float, float[16],
                int[16] > (one_over_hk, mu_stabk, cell_volumek, Uk, Vk, dPdFk,
                           dk);
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
          typedef T (&refArray1)[16];
          typedef T (&refArray2)[16];
          typedef T (&refArray3)[16];
          typedef T (&refArray4)[9][16];
          typedef T (&refArray5)[9][16];
          typedef T (&refArray6)[12][16];
          typedef T (&refArray7)[3][8][16];
          for (int i = 0; i < 16; i += 4)
            {
              refArray1 one_over_hk =
                reinterpret_cast < refArray1 > (one_over_h[i]);
              refArray2 mu_stabk = reinterpret_cast < refArray2 > (mu_stab[i]);
              refArray3 cell_volumek =
                reinterpret_cast < refArray3 > (cell_volume[i]);
              refArray4 Uk = reinterpret_cast < refArray4 > (U[0][i]);
              refArray5 Vk = reinterpret_cast < refArray5 > (V[0][i]);
              refArray6 dPdFk = reinterpret_cast < refArray6 > (dPdF[0][i]);
              refArray7 dk = reinterpret_cast < refArray7 > (d[0][0][i]);
              Compute_Diagonal_Contribution < __m128, float[16],
                int[16] > (one_over_hk, mu_stabk, cell_volumek, Uk, Vk, dPdFk,
                           dk);
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
          typedef T (&refArray1)[16];
          typedef T (&refArray2)[16];
          typedef T (&refArray3)[16];
          typedef T (&refArray4)[9][16];
          typedef T (&refArray5)[9][16];
          typedef T (&refArray6)[12][16];
          typedef T (&refArray7)[3][8][16];
          for (int i = 0; i < 16; i += 8)
            {
              refArray1 one_over_hk =
                reinterpret_cast < refArray1 > (one_over_h[i]);
              refArray2 mu_stabk = reinterpret_cast < refArray2 > (mu_stab[i]);
              refArray3 cell_volumek =
                reinterpret_cast < refArray3 > (cell_volume[i]);
              refArray4 Uk = reinterpret_cast < refArray4 > (U[0][i]);
              refArray5 Vk = reinterpret_cast < refArray5 > (V[0][i]);
              refArray6 dPdFk = reinterpret_cast < refArray6 > (dPdF[0][i]);
              refArray7 dk = reinterpret_cast < refArray7 > (d[0][0][i]);
              Compute_Diagonal_Contribution < __m256, float[16],
                int[16] > (one_over_hk, mu_stabk, cell_volumek, Uk, Vk, dPdFk,
                           dk);
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
          typedef T (&refArray1)[16];
          typedef T (&refArray2)[16];
          typedef T (&refArray3)[16];
          typedef T (&refArray4)[9][16];
          typedef T (&refArray5)[9][16];
          typedef T (&refArray6)[12][16];
          typedef T (&refArray7)[3][8][16];
          for (int i = 0; i < 16; i += 4)
            {
              refArray1 one_over_hk =
                reinterpret_cast < refArray1 > (one_over_h[i]);
              refArray2 mu_stabk = reinterpret_cast < refArray2 > (mu_stab[i]);
              refArray3 cell_volumek =
                reinterpret_cast < refArray3 > (cell_volume[i]);
              refArray4 Uk = reinterpret_cast < refArray4 > (U[0][i]);
              refArray5 Vk = reinterpret_cast < refArray5 > (V[0][i]);
              refArray6 dPdFk = reinterpret_cast < refArray6 > (dPdF[0][i]);
              refArray7 dk = reinterpret_cast < refArray7 > (d[0][0][i]);
              Compute_Diagonal_Contribution < float32x4_t, float[16],
                int[16] > (one_over_hk, mu_stabk, cell_volumek, Uk, Vk, dPdFk,
                           dk);
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
          typedef T (&refArray1)[16];
          typedef T (&refArray2)[16];
          typedef T (&refArray3)[16];
          typedef T (&refArray4)[9][16];
          typedef T (&refArray5)[9][16];
          typedef T (&refArray6)[12][16];
          typedef T (&refArray7)[3][8][16];
          for (int i = 0; i < 16; i += 16)
            {
              refArray1 one_over_hk =
                reinterpret_cast < refArray1 > (one_over_h[i]);
              refArray2 mu_stabk = reinterpret_cast < refArray2 > (mu_stab[i]);
              refArray3 cell_volumek =
                reinterpret_cast < refArray3 > (cell_volume[i]);
              refArray4 Uk = reinterpret_cast < refArray4 > (U[0][i]);
              refArray5 Vk = reinterpret_cast < refArray5 > (V[0][i]);
              refArray6 dPdFk = reinterpret_cast < refArray6 > (dPdF[0][i]);
              refArray7 dk = reinterpret_cast < refArray7 > (d[0][0][i]);
              Compute_Diagonal_Contribution < __m512, float[16],
                int[16] > (one_over_hk, mu_stabk, cell_volumek, Uk, Vk, dPdFk,
                           dk);
            }

        }
      stop_timer ();
      std::cout << get_time () << "s" << std::endl;
    }
#endif

  }



  return 0;

}

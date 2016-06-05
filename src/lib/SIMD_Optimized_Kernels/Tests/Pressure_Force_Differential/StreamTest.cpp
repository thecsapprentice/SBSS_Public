
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <sys/time.h>
#include "KernelCommon.h"
struct NEOHOOKEAN_TAG;
struct COROTATED_TAG;

#include "Pressure_Force_Differential.h"

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
    std::cout << "Running Stream Test for Pressure_Force_Differential " << std::
      endl;

//=======================================================
//
//        DEFINE ALL VARIABLES USED BY KERNEL
//
//=======================================================

    T dq[16] __attribute__ ((aligned (64)));
    T dq_reference[16] __attribute__ ((aligned (64)));
    T dq_original[16] __attribute__ ((aligned (64)));
    T Q_hat[3][16] __attribute__ ((aligned (64)));
    T dF_hat[9][16] __attribute__ ((aligned (64)));
    T dp[16] __attribute__ ((aligned (64)));
    T alpha[16] __attribute__ ((aligned (64)));
    T alpha_squared_over_kappa[16] __attribute__ ((aligned (64)));


    for (int __a = 0; __a < 16; __a++)
      {
        dq_original[__a] = Get_Random < float >();
        dq[__a] = dq_original[__a];
        dq_reference[__a] = dq_original[__a];
      }
    for (int __a = 0; __a < 3; __a++)
      for (int __b = 0; __b < 16; __b++)
        Q_hat[__a][__b] = Get_Random < float >();
    for (int __a = 0; __a < 9; __a++)
      for (int __b = 0; __b < 16; __b++)
        dF_hat[__a][__b] = Get_Random < float >();
    for (int __a = 0; __a < 16; __a++)
      dp[__a] = Get_Random < float >();
    for (int __a = 0; __a < 16; __a++)
      alpha[__a] = Get_Random < float >();
    for (int __a = 0; __a < 16; __a++)
      alpha_squared_over_kappa[__a] = Get_Random < float >();

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
          typedef T (&refArray2)[3][16];
          typedef T (&refArray3)[9][16];
          typedef T (&refArray4)[16];
          typedef T (&refArray5)[16];
          typedef T (&refArray6)[16];
          for (int i = 0; i < 16; i += 1)
            {
              refArray1 dqk = reinterpret_cast < refArray1 > (dq[i]);
              refArray2 Q_hatk = reinterpret_cast < refArray2 > (Q_hat[0][i]);
              refArray3 dF_hatk = reinterpret_cast < refArray3 > (dF_hat[0][i]);
              refArray4 dpk = reinterpret_cast < refArray4 > (dp[i]);
              refArray5 alphak = reinterpret_cast < refArray5 > (alpha[i]);
              refArray6 alpha_squared_over_kappak =
                reinterpret_cast < refArray6 > (alpha_squared_over_kappa[i]);
              Pressure_Force_Differential < float, float[16], int[16] > (dqk,
                                                                         Q_hatk,
                                                                         dF_hatk,
                                                                         dpk,
                                                                         alphak,
                                                                         alpha_squared_over_kappak);
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
          typedef T (&refArray2)[3][16];
          typedef T (&refArray3)[9][16];
          typedef T (&refArray4)[16];
          typedef T (&refArray5)[16];
          typedef T (&refArray6)[16];
          for (int i = 0; i < 16; i += 4)
            {
              refArray1 dqk = reinterpret_cast < refArray1 > (dq[i]);
              refArray2 Q_hatk = reinterpret_cast < refArray2 > (Q_hat[0][i]);
              refArray3 dF_hatk = reinterpret_cast < refArray3 > (dF_hat[0][i]);
              refArray4 dpk = reinterpret_cast < refArray4 > (dp[i]);
              refArray5 alphak = reinterpret_cast < refArray5 > (alpha[i]);
              refArray6 alpha_squared_over_kappak =
                reinterpret_cast < refArray6 > (alpha_squared_over_kappa[i]);
              Pressure_Force_Differential < __m128, float[16], int[16] > (dqk,
                                                                          Q_hatk,
                                                                          dF_hatk,
                                                                          dpk,
                                                                          alphak,
                                                                          alpha_squared_over_kappak);
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
          typedef T (&refArray2)[3][16];
          typedef T (&refArray3)[9][16];
          typedef T (&refArray4)[16];
          typedef T (&refArray5)[16];
          typedef T (&refArray6)[16];
          for (int i = 0; i < 16; i += 8)
            {
              refArray1 dqk = reinterpret_cast < refArray1 > (dq[i]);
              refArray2 Q_hatk = reinterpret_cast < refArray2 > (Q_hat[0][i]);
              refArray3 dF_hatk = reinterpret_cast < refArray3 > (dF_hat[0][i]);
              refArray4 dpk = reinterpret_cast < refArray4 > (dp[i]);
              refArray5 alphak = reinterpret_cast < refArray5 > (alpha[i]);
              refArray6 alpha_squared_over_kappak =
                reinterpret_cast < refArray6 > (alpha_squared_over_kappa[i]);
              Pressure_Force_Differential < __m256, float[16], int[16] > (dqk,
                                                                          Q_hatk,
                                                                          dF_hatk,
                                                                          dpk,
                                                                          alphak,
                                                                          alpha_squared_over_kappak);
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
          typedef T (&refArray2)[3][16];
          typedef T (&refArray3)[9][16];
          typedef T (&refArray4)[16];
          typedef T (&refArray5)[16];
          typedef T (&refArray6)[16];
          for (int i = 0; i < 16; i += 4)
            {
              refArray1 dqk = reinterpret_cast < refArray1 > (dq[i]);
              refArray2 Q_hatk = reinterpret_cast < refArray2 > (Q_hat[0][i]);
              refArray3 dF_hatk = reinterpret_cast < refArray3 > (dF_hat[0][i]);
              refArray4 dpk = reinterpret_cast < refArray4 > (dp[i]);
              refArray5 alphak = reinterpret_cast < refArray5 > (alpha[i]);
              refArray6 alpha_squared_over_kappak =
                reinterpret_cast < refArray6 > (alpha_squared_over_kappa[i]);
              Pressure_Force_Differential < float32x4_t, float[16],
                int[16] > (dqk, Q_hatk, dF_hatk, dpk, alphak,
                           alpha_squared_over_kappak);
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
          typedef T (&refArray2)[3][16];
          typedef T (&refArray3)[9][16];
          typedef T (&refArray4)[16];
          typedef T (&refArray5)[16];
          typedef T (&refArray6)[16];
          for (int i = 0; i < 16; i += 16)
            {
              refArray1 dqk = reinterpret_cast < refArray1 > (dq[i]);
              refArray2 Q_hatk = reinterpret_cast < refArray2 > (Q_hat[0][i]);
              refArray3 dF_hatk = reinterpret_cast < refArray3 > (dF_hat[0][i]);
              refArray4 dpk = reinterpret_cast < refArray4 > (dp[i]);
              refArray5 alphak = reinterpret_cast < refArray5 > (alpha[i]);
              refArray6 alpha_squared_over_kappak =
                reinterpret_cast < refArray6 > (alpha_squared_over_kappa[i]);
              Pressure_Force_Differential < __m512, float[16], int[16] > (dqk,
                                                                          Q_hatk,
                                                                          dF_hatk,
                                                                          dpk,
                                                                          alphak,
                                                                          alpha_squared_over_kappak);
            }

        }
      stop_timer ();
      std::cout << get_time () << "s" << std::endl;
    }
#endif

  }



  return 0;

}

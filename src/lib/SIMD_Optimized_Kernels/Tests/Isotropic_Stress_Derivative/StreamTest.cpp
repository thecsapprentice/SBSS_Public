
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <sys/time.h>
#include "KernelCommon.h"
struct NEOHOOKEAN_TAG;
struct COROTATED_TAG;

#include "Isotropic_Stress_Derivative.h"

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
    std::cout << "Running Stream Test for Isotropic_Stress_Derivative " << std::
      endl;

//=======================================================
//
//        DEFINE ALL VARIABLES USED BY KERNEL
//
//=======================================================

    T dPdF[12][16] __attribute__ ((aligned (64)));
    T dPdF_reference[12][16] __attribute__ ((aligned (64)));
    T dPdF_original[12][16] __attribute__ ((aligned (64)));
    T Sigma[3][16] __attribute__ ((aligned (64)));
    T p[16] __attribute__ ((aligned (64)));
    T mu[16] __attribute__ ((aligned (64)));
    T kappa[16] __attribute__ ((aligned (64)));
    T alpha[16] __attribute__ ((aligned (64)));
    bool apply_definiteness_fix __attribute__ ((aligned (64)));


    for (int __a = 0; __a < 12; __a++)
      for (int __b = 0; __b < 16; __b++)
        {
          dPdF_original[__a][__b] = Get_Random < float >();
          dPdF[__a][__b] = dPdF_original[__a][__b];
          dPdF_reference[__a][__b] = dPdF_original[__a][__b];
        }
    for (int __a = 0; __a < 3; __a++)
      for (int __b = 0; __b < 16; __b++)
        Sigma[__a][__b] = Get_Random < float >();
    for (int __a = 0; __a < 16; __a++)
      p[__a] = Get_Random < float >();
    for (int __a = 0; __a < 16; __a++)
      mu[__a] = Get_Random < float >();
    for (int __a = 0; __a < 16; __a++)
      kappa[__a] = Get_Random < float >();
    for (int __a = 0; __a < 16; __a++)
      alpha[__a] = Get_Random < float >();
    apply_definiteness_fix = true;

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
          typedef T (&refArray1)[12][16];
          typedef T (&refArray2)[3][16];
          typedef T (&refArray3)[16];
          typedef T (&refArray4)[16];
          typedef T (&refArray5)[16];
          typedef T (&refArray6)[16];
          typedef bool (&refArray7);
          for (int i = 0; i < 16; i += 1)
            {
              refArray1 dPdFk = reinterpret_cast < refArray1 > (dPdF[0][i]);
              refArray2 Sigmak = reinterpret_cast < refArray2 > (Sigma[0][i]);
              refArray3 pk = reinterpret_cast < refArray3 > (p[i]);
              refArray4 muk = reinterpret_cast < refArray4 > (mu[i]);
              refArray5 kappak = reinterpret_cast < refArray5 > (kappa[i]);
              refArray6 alphak = reinterpret_cast < refArray6 > (alpha[i]);
              refArray7 apply_definiteness_fixk =
                reinterpret_cast < refArray7 > (apply_definiteness_fix);
              Isotropic_Stress_Derivative < NEOHOOKEAN_TAG, float, float[16],
                int[16] >::Run (dPdFk, Sigmak, pk, muk, kappak, alphak,
                                apply_definiteness_fixk);
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
          typedef T (&refArray1)[12][16];
          typedef T (&refArray2)[3][16];
          typedef T (&refArray3)[16];
          typedef T (&refArray4)[16];
          typedef T (&refArray5)[16];
          typedef T (&refArray6)[16];
          typedef bool (&refArray7);
          for (int i = 0; i < 16; i += 4)
            {
              refArray1 dPdFk = reinterpret_cast < refArray1 > (dPdF[0][i]);
              refArray2 Sigmak = reinterpret_cast < refArray2 > (Sigma[0][i]);
              refArray3 pk = reinterpret_cast < refArray3 > (p[i]);
              refArray4 muk = reinterpret_cast < refArray4 > (mu[i]);
              refArray5 kappak = reinterpret_cast < refArray5 > (kappa[i]);
              refArray6 alphak = reinterpret_cast < refArray6 > (alpha[i]);
              refArray7 apply_definiteness_fixk =
                reinterpret_cast < refArray7 > (apply_definiteness_fix);
              Isotropic_Stress_Derivative < NEOHOOKEAN_TAG, __m128, float[16],
                int[16] >::Run (dPdFk, Sigmak, pk, muk, kappak, alphak,
                                apply_definiteness_fixk);
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
          typedef T (&refArray1)[12][16];
          typedef T (&refArray2)[3][16];
          typedef T (&refArray3)[16];
          typedef T (&refArray4)[16];
          typedef T (&refArray5)[16];
          typedef T (&refArray6)[16];
          typedef bool (&refArray7);
          for (int i = 0; i < 16; i += 8)
            {
              refArray1 dPdFk = reinterpret_cast < refArray1 > (dPdF[0][i]);
              refArray2 Sigmak = reinterpret_cast < refArray2 > (Sigma[0][i]);
              refArray3 pk = reinterpret_cast < refArray3 > (p[i]);
              refArray4 muk = reinterpret_cast < refArray4 > (mu[i]);
              refArray5 kappak = reinterpret_cast < refArray5 > (kappa[i]);
              refArray6 alphak = reinterpret_cast < refArray6 > (alpha[i]);
              refArray7 apply_definiteness_fixk =
                reinterpret_cast < refArray7 > (apply_definiteness_fix);
              Isotropic_Stress_Derivative < NEOHOOKEAN_TAG, __m256, float[16],
                int[16] >::Run (dPdFk, Sigmak, pk, muk, kappak, alphak,
                                apply_definiteness_fixk);
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
          typedef T (&refArray1)[12][16];
          typedef T (&refArray2)[3][16];
          typedef T (&refArray3)[16];
          typedef T (&refArray4)[16];
          typedef T (&refArray5)[16];
          typedef T (&refArray6)[16];
          typedef bool (&refArray7);
          for (int i = 0; i < 16; i += 4)
            {
              refArray1 dPdFk = reinterpret_cast < refArray1 > (dPdF[0][i]);
              refArray2 Sigmak = reinterpret_cast < refArray2 > (Sigma[0][i]);
              refArray3 pk = reinterpret_cast < refArray3 > (p[i]);
              refArray4 muk = reinterpret_cast < refArray4 > (mu[i]);
              refArray5 kappak = reinterpret_cast < refArray5 > (kappa[i]);
              refArray6 alphak = reinterpret_cast < refArray6 > (alpha[i]);
              refArray7 apply_definiteness_fixk =
                reinterpret_cast < refArray7 > (apply_definiteness_fix);
              Isotropic_Stress_Derivative < NEOHOOKEAN_TAG, float32x4_t,
                float[16], int[16] >::Run (dPdFk, Sigmak, pk, muk, kappak,
                                           alphak, apply_definiteness_fixk);
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
          typedef T (&refArray1)[12][16];
          typedef T (&refArray2)[3][16];
          typedef T (&refArray3)[16];
          typedef T (&refArray4)[16];
          typedef T (&refArray5)[16];
          typedef T (&refArray6)[16];
          typedef bool (&refArray7);
          for (int i = 0; i < 16; i += 16)
            {
              refArray1 dPdFk = reinterpret_cast < refArray1 > (dPdF[0][i]);
              refArray2 Sigmak = reinterpret_cast < refArray2 > (Sigma[0][i]);
              refArray3 pk = reinterpret_cast < refArray3 > (p[i]);
              refArray4 muk = reinterpret_cast < refArray4 > (mu[i]);
              refArray5 kappak = reinterpret_cast < refArray5 > (kappa[i]);
              refArray6 alphak = reinterpret_cast < refArray6 > (alpha[i]);
              refArray7 apply_definiteness_fixk =
                reinterpret_cast < refArray7 > (apply_definiteness_fix);
              Isotropic_Stress_Derivative < NEOHOOKEAN_TAG, __m512, float[16],
                int[16] >::Run (dPdFk, Sigmak, pk, muk, kappak, alphak,
                                apply_definiteness_fixk);
            }

        }
      stop_timer ();
      std::cout << get_time () << "s" << std::endl;
    }
#endif

  }



  {
    std::cout << "Running Stream Test for Isotropic_Stress_Derivative " << std::
      endl;

//=======================================================
//
//        DEFINE ALL VARIABLES USED BY KERNEL
//
//=======================================================

    T dPdF[12][16] __attribute__ ((aligned (64)));
    T dPdF_reference[12][16] __attribute__ ((aligned (64)));
    T dPdF_original[12][16] __attribute__ ((aligned (64)));
    T Sigma[3][16] __attribute__ ((aligned (64)));
    T p[16] __attribute__ ((aligned (64)));
    T mu[16] __attribute__ ((aligned (64)));
    T kappa[16] __attribute__ ((aligned (64)));
    T alpha[16] __attribute__ ((aligned (64)));
    bool apply_definiteness_fix __attribute__ ((aligned (64)));


    for (int __a = 0; __a < 12; __a++)
      for (int __b = 0; __b < 16; __b++)
        {
          dPdF_original[__a][__b] = Get_Random < float >();
          dPdF[__a][__b] = dPdF_original[__a][__b];
          dPdF_reference[__a][__b] = dPdF_original[__a][__b];
        }
    for (int __a = 0; __a < 3; __a++)
      for (int __b = 0; __b < 16; __b++)
        Sigma[__a][__b] = Get_Random < float >();
    for (int __a = 0; __a < 16; __a++)
      p[__a] = Get_Random < float >();
    for (int __a = 0; __a < 16; __a++)
      mu[__a] = Get_Random < float >();
    for (int __a = 0; __a < 16; __a++)
      kappa[__a] = Get_Random < float >();
    for (int __a = 0; __a < 16; __a++)
      alpha[__a] = Get_Random < float >();
    apply_definiteness_fix = true;

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
          typedef T (&refArray1)[12][16];
          typedef T (&refArray2)[3][16];
          typedef T (&refArray3)[16];
          typedef T (&refArray4)[16];
          typedef T (&refArray5)[16];
          typedef T (&refArray6)[16];
          typedef bool (&refArray7);
          for (int i = 0; i < 16; i += 1)
            {
              refArray1 dPdFk = reinterpret_cast < refArray1 > (dPdF[0][i]);
              refArray2 Sigmak = reinterpret_cast < refArray2 > (Sigma[0][i]);
              refArray3 pk = reinterpret_cast < refArray3 > (p[i]);
              refArray4 muk = reinterpret_cast < refArray4 > (mu[i]);
              refArray5 kappak = reinterpret_cast < refArray5 > (kappa[i]);
              refArray6 alphak = reinterpret_cast < refArray6 > (alpha[i]);
              refArray7 apply_definiteness_fixk =
                reinterpret_cast < refArray7 > (apply_definiteness_fix);
              Isotropic_Stress_Derivative < COROTATED_TAG, float, float[16],
                int[16] >::Run (dPdFk, Sigmak, pk, muk, kappak, alphak,
                                apply_definiteness_fixk);
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
          typedef T (&refArray1)[12][16];
          typedef T (&refArray2)[3][16];
          typedef T (&refArray3)[16];
          typedef T (&refArray4)[16];
          typedef T (&refArray5)[16];
          typedef T (&refArray6)[16];
          typedef bool (&refArray7);
          for (int i = 0; i < 16; i += 4)
            {
              refArray1 dPdFk = reinterpret_cast < refArray1 > (dPdF[0][i]);
              refArray2 Sigmak = reinterpret_cast < refArray2 > (Sigma[0][i]);
              refArray3 pk = reinterpret_cast < refArray3 > (p[i]);
              refArray4 muk = reinterpret_cast < refArray4 > (mu[i]);
              refArray5 kappak = reinterpret_cast < refArray5 > (kappa[i]);
              refArray6 alphak = reinterpret_cast < refArray6 > (alpha[i]);
              refArray7 apply_definiteness_fixk =
                reinterpret_cast < refArray7 > (apply_definiteness_fix);
              Isotropic_Stress_Derivative < COROTATED_TAG, __m128, float[16],
                int[16] >::Run (dPdFk, Sigmak, pk, muk, kappak, alphak,
                                apply_definiteness_fixk);
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
          typedef T (&refArray1)[12][16];
          typedef T (&refArray2)[3][16];
          typedef T (&refArray3)[16];
          typedef T (&refArray4)[16];
          typedef T (&refArray5)[16];
          typedef T (&refArray6)[16];
          typedef bool (&refArray7);
          for (int i = 0; i < 16; i += 8)
            {
              refArray1 dPdFk = reinterpret_cast < refArray1 > (dPdF[0][i]);
              refArray2 Sigmak = reinterpret_cast < refArray2 > (Sigma[0][i]);
              refArray3 pk = reinterpret_cast < refArray3 > (p[i]);
              refArray4 muk = reinterpret_cast < refArray4 > (mu[i]);
              refArray5 kappak = reinterpret_cast < refArray5 > (kappa[i]);
              refArray6 alphak = reinterpret_cast < refArray6 > (alpha[i]);
              refArray7 apply_definiteness_fixk =
                reinterpret_cast < refArray7 > (apply_definiteness_fix);
              Isotropic_Stress_Derivative < COROTATED_TAG, __m256, float[16],
                int[16] >::Run (dPdFk, Sigmak, pk, muk, kappak, alphak,
                                apply_definiteness_fixk);
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
          typedef T (&refArray1)[12][16];
          typedef T (&refArray2)[3][16];
          typedef T (&refArray3)[16];
          typedef T (&refArray4)[16];
          typedef T (&refArray5)[16];
          typedef T (&refArray6)[16];
          typedef bool (&refArray7);
          for (int i = 0; i < 16; i += 4)
            {
              refArray1 dPdFk = reinterpret_cast < refArray1 > (dPdF[0][i]);
              refArray2 Sigmak = reinterpret_cast < refArray2 > (Sigma[0][i]);
              refArray3 pk = reinterpret_cast < refArray3 > (p[i]);
              refArray4 muk = reinterpret_cast < refArray4 > (mu[i]);
              refArray5 kappak = reinterpret_cast < refArray5 > (kappa[i]);
              refArray6 alphak = reinterpret_cast < refArray6 > (alpha[i]);
              refArray7 apply_definiteness_fixk =
                reinterpret_cast < refArray7 > (apply_definiteness_fix);
              Isotropic_Stress_Derivative < COROTATED_TAG, float32x4_t,
                float[16], int[16] >::Run (dPdFk, Sigmak, pk, muk, kappak,
                                           alphak, apply_definiteness_fixk);
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
          typedef T (&refArray1)[12][16];
          typedef T (&refArray2)[3][16];
          typedef T (&refArray3)[16];
          typedef T (&refArray4)[16];
          typedef T (&refArray5)[16];
          typedef T (&refArray6)[16];
          typedef bool (&refArray7);
          for (int i = 0; i < 16; i += 16)
            {
              refArray1 dPdFk = reinterpret_cast < refArray1 > (dPdF[0][i]);
              refArray2 Sigmak = reinterpret_cast < refArray2 > (Sigma[0][i]);
              refArray3 pk = reinterpret_cast < refArray3 > (p[i]);
              refArray4 muk = reinterpret_cast < refArray4 > (mu[i]);
              refArray5 kappak = reinterpret_cast < refArray5 > (kappa[i]);
              refArray6 alphak = reinterpret_cast < refArray6 > (alpha[i]);
              refArray7 apply_definiteness_fixk =
                reinterpret_cast < refArray7 > (apply_definiteness_fix);
              Isotropic_Stress_Derivative < COROTATED_TAG, __m512, float[16],
                int[16] >::Run (dPdFk, Sigmak, pk, muk, kappak, alphak,
                                apply_definiteness_fixk);
            }

        }
      stop_timer ();
      std::cout << get_time () << "s" << std::endl;
    }
#endif

  }



  return 0;

}

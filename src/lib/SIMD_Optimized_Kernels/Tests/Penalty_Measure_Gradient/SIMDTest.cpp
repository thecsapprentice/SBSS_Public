
#include <cstdlib>
#include <iostream>
#include "KernelCommon.h"

struct NEOHOOKEAN_TAG;
struct COROTATED_TAG;

#include "Penalty_Measure_Gradient.h"
#include "Penalty_Measure_Gradient_Reference.h"

template < class T > T Get_Random (const T a = (T) - 1., const T b = (T) 1.)
{
  return ((b - a) * (T) rand ()) / (T) RAND_MAX + a;
}

int
main (int argc, char *argv[])
{
  typedef float   T;

  int             seed = 1;
  if (argc == 2)
    seed = atoi (argv[1]);
  srand (seed);

  std::cout.precision (10);
  std::cout.setf (std::ios::fixed, std::ios::floatfield);



  {
    std::cout << "Running SIMD Test for Penalty_Measure_Gradient " << std::endl;


//=======================================================
//
//        DEFINE ALL VARIABLES USED BY KERNEL
//
//=======================================================

    T               Sigma[3][16] __attribute__ ((aligned (64)));
    T               Q_hat[3][16] __attribute__ ((aligned (64)));
    T               Q_hat_reference[3][16] __attribute__ ((aligned (64)));
    T               Q_hat_original[3][16] __attribute__ ((aligned (64)));


    for (int __a = 0; __a < 3; __a++)
      for (int __b = 0; __b < 16; __b++)
        Sigma[__a][__b] = Get_Random < float >();
    for (int __a = 0; __a < 3; __a++)
      for (int __b = 0; __b < 16; __b++)
      {
        Q_hat_original[__a][__b] = Get_Random < float >();
        Q_hat[__a][__b] = Q_hat_original[__a][__b];
        Q_hat_reference[__a][__b] = Q_hat_original[__a][__b];
      }


//=======================================================
//
//             COMPUTE REFERENCE RESULTS
//
//=======================================================

    T               __mSigma[3] __attribute__ ((aligned (4)));
    T               __mQ_hat[3] __attribute__ ((aligned (4)));
    T               __mQ_hat_reference[3] __attribute__ ((aligned (4)));
    T               __mQ_hat_original[3] __attribute__ ((aligned (4)));
    for (int k = 0; k < 16; k++)
    {
      for (int __a = 0; __a < 3; __a++)
        __mSigma[__a] = Sigma[__a][k];
      for (int __a = 0; __a < 3; __a++)
        __mQ_hat_reference[__a] = Q_hat_reference[__a][k];
      Penalty_Measure_Gradient < NEOHOOKEAN_TAG, float, float,
        int >::Run (__mSigma, __mQ_hat_reference);
      for (int __a = 0; __a < 3; __a++)
        Q_hat_reference[__a][k] = __mQ_hat_reference[__a];
    }

//=======================================================
//
//               COMPUTE SCALAR RESULTS
//
//=======================================================

    {
      typedef         T (&refArray1)[3][16];
      typedef         T (&refArray2)[3][16];
      for (int __a = 0; __a < 3; __a++)
        for (int __b = 0; __b < 16; __b++)
          Q_hat[__a][__b] = Q_hat_original[__a][__b];
      for (int i = 0; i < 16; i += 1)
      {
        refArray1       Sigmak = reinterpret_cast < refArray1 > (Sigma[0][i]);
        refArray2       Q_hatk = reinterpret_cast < refArray2 > (Q_hat[0][i]);
        Penalty_Measure_Gradient < NEOHOOKEAN_TAG, float, float[16],
          int[16] >::Run (Sigmak, Q_hatk);
      }
      for (int __a = 0; __a < 3; __a++)
        for (int __b = 0; __b < 16; __b++)
          if (std::
              abs ((Q_hat[__a][__b] -
                    Q_hat_reference[__a][__b]) / (Q_hat_reference[__a][__b])) >
              1)
          {
            std::cerr << "Mismatch detected in SCALAR implementation" << std::
              endl;
            std::cerr << "Variable Q_hat:" << std::endl;
            std::
              cerr << "seed=" << seed << ", __a=" << __a << ", __b=" << __b <<
              std::endl;
            std::cerr << "Q_hat SCALAR=  " << Q_hat[__a][__b] << std::endl;
            std::
              cerr << "Q_hat Reference=  " << Q_hat_reference[__a][__b] << std::
              endl;
            std::cerr << "Q_hat Rel Difference=  " << std::
              abs ((Q_hat[__a][__b] -
                    Q_hat_reference[__a][__b]) /
                   (Q_hat_reference[__a][__b])) << std::endl;
            std::cerr << "Q_hat Abs Difference=  " << std::
              abs (Q_hat[__a][__b] - Q_hat_reference[__a][__b]) << std::endl;
            return 1;
          }

    }

//=======================================================
//
//               COMPUTE SSE RESULTS
//
//=======================================================

#ifdef ENABLE_SSE_INSTRUCTION_SET
    {
      typedef         T (&refArray1)[3][16];
      typedef         T (&refArray2)[3][16];
      for (int __a = 0; __a < 3; __a++)
        for (int __b = 0; __b < 16; __b++)
          Q_hat[__a][__b] = Q_hat_original[__a][__b];
      for (int i = 0; i < 16; i += 4)
      {
        refArray1       Sigmak = reinterpret_cast < refArray1 > (Sigma[0][i]);
        refArray2       Q_hatk = reinterpret_cast < refArray2 > (Q_hat[0][i]);
        Penalty_Measure_Gradient < NEOHOOKEAN_TAG, __m128, float[16],
          int[16] >::Run (Sigmak, Q_hatk);
      }
      for (int __a = 0; __a < 3; __a++)
        for (int __b = 0; __b < 16; __b++)
          if (std::
              abs ((Q_hat[__a][__b] -
                    Q_hat_reference[__a][__b]) / (Q_hat_reference[__a][__b])) >
              1)
          {
            std::cerr << "Mismatch detected in SSE implementation" << std::endl;
            std::cerr << "Variable Q_hat:" << std::endl;
            std::
              cerr << "seed=" << seed << ", __a=" << __a << ", __b=" << __b <<
              std::endl;
            std::cerr << "Q_hat SSE=  " << Q_hat[__a][__b] << std::endl;
            std::
              cerr << "Q_hat Reference=  " << Q_hat_reference[__a][__b] << std::
              endl;
            std::cerr << "Q_hat Rel Difference=  " << std::
              abs ((Q_hat[__a][__b] -
                    Q_hat_reference[__a][__b]) /
                   (Q_hat_reference[__a][__b])) << std::endl;
            std::cerr << "Q_hat Abs Difference=  " << std::
              abs (Q_hat[__a][__b] - Q_hat_reference[__a][__b]) << std::endl;
            return 1;
          }

    }
#endif

//=======================================================
//
//               COMPUTE AVX RESULTS
//
//=======================================================

#ifdef ENABLE_AVX_INSTRUCTION_SET
    {
      typedef         T (&refArray1)[3][16];
      typedef         T (&refArray2)[3][16];
      for (int __a = 0; __a < 3; __a++)
        for (int __b = 0; __b < 16; __b++)
          Q_hat[__a][__b] = Q_hat_original[__a][__b];
      for (int i = 0; i < 16; i += 8)
      {
        refArray1       Sigmak = reinterpret_cast < refArray1 > (Sigma[0][i]);
        refArray2       Q_hatk = reinterpret_cast < refArray2 > (Q_hat[0][i]);
        Penalty_Measure_Gradient < NEOHOOKEAN_TAG, __m256, float[16],
          int[16] >::Run (Sigmak, Q_hatk);
      }
      for (int __a = 0; __a < 3; __a++)
        for (int __b = 0; __b < 16; __b++)
          if (std::
              abs ((Q_hat[__a][__b] -
                    Q_hat_reference[__a][__b]) / (Q_hat_reference[__a][__b])) >
              1)
          {
            std::cerr << "Mismatch detected in AVX implementation" << std::endl;
            std::cerr << "Variable Q_hat:" << std::endl;
            std::
              cerr << "seed=" << seed << ", __a=" << __a << ", __b=" << __b <<
              std::endl;
            std::cerr << "Q_hat AVX=  " << Q_hat[__a][__b] << std::endl;
            std::
              cerr << "Q_hat Reference=  " << Q_hat_reference[__a][__b] << std::
              endl;
            std::cerr << "Q_hat Rel Difference=  " << std::
              abs ((Q_hat[__a][__b] -
                    Q_hat_reference[__a][__b]) /
                   (Q_hat_reference[__a][__b])) << std::endl;
            std::cerr << "Q_hat Abs Difference=  " << std::
              abs (Q_hat[__a][__b] - Q_hat_reference[__a][__b]) << std::endl;
            return 1;
          }

    }
#endif

//=======================================================
//
//               COMPUTE NEON RESULTS
//
//=======================================================

#ifdef ENABLE_NEON_INSTRUCTION_SET
    {
      typedef         T (&refArray1)[3][16];
      typedef         T (&refArray2)[3][16];
      for (int __a = 0; __a < 3; __a++)
        for (int __b = 0; __b < 16; __b++)
          Q_hat[__a][__b] = Q_hat_original[__a][__b];
      for (int i = 0; i < 16; i += 4)
      {
        refArray1       Sigmak = reinterpret_cast < refArray1 > (Sigma[0][i]);
        refArray2       Q_hatk = reinterpret_cast < refArray2 > (Q_hat[0][i]);
        Penalty_Measure_Gradient < NEOHOOKEAN_TAG, float32x4_t, float[16],
          int[16] >::Run (Sigmak, Q_hatk);
      }
      for (int __a = 0; __a < 3; __a++)
        for (int __b = 0; __b < 16; __b++)
          if (std::
              abs ((Q_hat[__a][__b] -
                    Q_hat_reference[__a][__b]) / (Q_hat_reference[__a][__b])) >
              1)
          {
            std::cerr << "Mismatch detected in NEON implementation" << std::
              endl;
            std::cerr << "Variable Q_hat:" << std::endl;
            std::
              cerr << "seed=" << seed << ", __a=" << __a << ", __b=" << __b <<
              std::endl;
            std::cerr << "Q_hat NEON=  " << Q_hat[__a][__b] << std::endl;
            std::
              cerr << "Q_hat Reference=  " << Q_hat_reference[__a][__b] << std::
              endl;
            std::cerr << "Q_hat Rel Difference=  " << std::
              abs ((Q_hat[__a][__b] -
                    Q_hat_reference[__a][__b]) /
                   (Q_hat_reference[__a][__b])) << std::endl;
            std::cerr << "Q_hat Abs Difference=  " << std::
              abs (Q_hat[__a][__b] - Q_hat_reference[__a][__b]) << std::endl;
            return 1;
          }

    }
#endif

//=======================================================
//
//               COMPUTE MIC RESULTS
//
//=======================================================

#ifdef ENABLE_MIC_INSTRUCTION_SET
    {
      typedef         T (&refArray1)[3][16];
      typedef         T (&refArray2)[3][16];
      for (int __a = 0; __a < 3; __a++)
        for (int __b = 0; __b < 16; __b++)
          Q_hat[__a][__b] = Q_hat_original[__a][__b];
      for (int i = 0; i < 16; i += 16)
      {
        refArray1       Sigmak = reinterpret_cast < refArray1 > (Sigma[0][i]);
        refArray2       Q_hatk = reinterpret_cast < refArray2 > (Q_hat[0][i]);
        Penalty_Measure_Gradient < NEOHOOKEAN_TAG, __m512, float[16],
          int[16] >::Run (Sigmak, Q_hatk);
      }
      for (int __a = 0; __a < 3; __a++)
        for (int __b = 0; __b < 16; __b++)
          if (std::
              abs ((Q_hat[__a][__b] -
                    Q_hat_reference[__a][__b]) / (Q_hat_reference[__a][__b])) >
              1)
          {
            std::cerr << "Mismatch detected in MIC implementation" << std::endl;
            std::cerr << "Variable Q_hat:" << std::endl;
            std::
              cerr << "seed=" << seed << ", __a=" << __a << ", __b=" << __b <<
              std::endl;
            std::cerr << "Q_hat MIC=  " << Q_hat[__a][__b] << std::endl;
            std::
              cerr << "Q_hat Reference=  " << Q_hat_reference[__a][__b] << std::
              endl;
            std::cerr << "Q_hat Rel Difference=  " << std::
              abs ((Q_hat[__a][__b] -
                    Q_hat_reference[__a][__b]) /
                   (Q_hat_reference[__a][__b])) << std::endl;
            std::cerr << "Q_hat Abs Difference=  " << std::
              abs (Q_hat[__a][__b] - Q_hat_reference[__a][__b]) << std::endl;
            return 1;
          }

    }
#endif

  }



  {
    std::cout << "Running SIMD Test for Penalty_Measure_Gradient " << std::endl;


//=======================================================
//
//        DEFINE ALL VARIABLES USED BY KERNEL
//
//=======================================================

    T               Sigma[3][16] __attribute__ ((aligned (64)));
    T               Q_hat[3][16] __attribute__ ((aligned (64)));
    T               Q_hat_reference[3][16] __attribute__ ((aligned (64)));
    T               Q_hat_original[3][16] __attribute__ ((aligned (64)));


    for (int __a = 0; __a < 3; __a++)
      for (int __b = 0; __b < 16; __b++)
        Sigma[__a][__b] = Get_Random < float >();
    for (int __a = 0; __a < 3; __a++)
      for (int __b = 0; __b < 16; __b++)
      {
        Q_hat_original[__a][__b] = Get_Random < float >();
        Q_hat[__a][__b] = Q_hat_original[__a][__b];
        Q_hat_reference[__a][__b] = Q_hat_original[__a][__b];
      }


//=======================================================
//
//             COMPUTE REFERENCE RESULTS
//
//=======================================================

    T               __mSigma[3] __attribute__ ((aligned (4)));
    T               __mQ_hat[3] __attribute__ ((aligned (4)));
    T               __mQ_hat_reference[3] __attribute__ ((aligned (4)));
    T               __mQ_hat_original[3] __attribute__ ((aligned (4)));
    for (int k = 0; k < 16; k++)
    {
      for (int __a = 0; __a < 3; __a++)
        __mSigma[__a] = Sigma[__a][k];
      for (int __a = 0; __a < 3; __a++)
        __mQ_hat_reference[__a] = Q_hat_reference[__a][k];
      Penalty_Measure_Gradient < COROTATED_TAG, float, float,
        int >::Run (__mSigma, __mQ_hat_reference);
      for (int __a = 0; __a < 3; __a++)
        Q_hat_reference[__a][k] = __mQ_hat_reference[__a];
    }

//=======================================================
//
//               COMPUTE SCALAR RESULTS
//
//=======================================================

    {
      typedef         T (&refArray1)[3][16];
      typedef         T (&refArray2)[3][16];
      for (int __a = 0; __a < 3; __a++)
        for (int __b = 0; __b < 16; __b++)
          Q_hat[__a][__b] = Q_hat_original[__a][__b];
      for (int i = 0; i < 16; i += 1)
      {
        refArray1       Sigmak = reinterpret_cast < refArray1 > (Sigma[0][i]);
        refArray2       Q_hatk = reinterpret_cast < refArray2 > (Q_hat[0][i]);
        Penalty_Measure_Gradient < COROTATED_TAG, float, float[16],
          int[16] >::Run (Sigmak, Q_hatk);
      }
      for (int __a = 0; __a < 3; __a++)
        for (int __b = 0; __b < 16; __b++)
          if (std::
              abs ((Q_hat[__a][__b] -
                    Q_hat_reference[__a][__b]) / (Q_hat_reference[__a][__b])) >
              1)
          {
            std::cerr << "Mismatch detected in SCALAR implementation" << std::
              endl;
            std::cerr << "Variable Q_hat:" << std::endl;
            std::
              cerr << "seed=" << seed << ", __a=" << __a << ", __b=" << __b <<
              std::endl;
            std::cerr << "Q_hat SCALAR=  " << Q_hat[__a][__b] << std::endl;
            std::
              cerr << "Q_hat Reference=  " << Q_hat_reference[__a][__b] << std::
              endl;
            std::cerr << "Q_hat Rel Difference=  " << std::
              abs ((Q_hat[__a][__b] -
                    Q_hat_reference[__a][__b]) /
                   (Q_hat_reference[__a][__b])) << std::endl;
            std::cerr << "Q_hat Abs Difference=  " << std::
              abs (Q_hat[__a][__b] - Q_hat_reference[__a][__b]) << std::endl;
            return 1;
          }

    }

//=======================================================
//
//               COMPUTE SSE RESULTS
//
//=======================================================

#ifdef ENABLE_SSE_INSTRUCTION_SET
    {
      typedef         T (&refArray1)[3][16];
      typedef         T (&refArray2)[3][16];
      for (int __a = 0; __a < 3; __a++)
        for (int __b = 0; __b < 16; __b++)
          Q_hat[__a][__b] = Q_hat_original[__a][__b];
      for (int i = 0; i < 16; i += 4)
      {
        refArray1       Sigmak = reinterpret_cast < refArray1 > (Sigma[0][i]);
        refArray2       Q_hatk = reinterpret_cast < refArray2 > (Q_hat[0][i]);
        Penalty_Measure_Gradient < COROTATED_TAG, __m128, float[16],
          int[16] >::Run (Sigmak, Q_hatk);
      }
      for (int __a = 0; __a < 3; __a++)
        for (int __b = 0; __b < 16; __b++)
          if (std::
              abs ((Q_hat[__a][__b] -
                    Q_hat_reference[__a][__b]) / (Q_hat_reference[__a][__b])) >
              1)
          {
            std::cerr << "Mismatch detected in SSE implementation" << std::endl;
            std::cerr << "Variable Q_hat:" << std::endl;
            std::
              cerr << "seed=" << seed << ", __a=" << __a << ", __b=" << __b <<
              std::endl;
            std::cerr << "Q_hat SSE=  " << Q_hat[__a][__b] << std::endl;
            std::
              cerr << "Q_hat Reference=  " << Q_hat_reference[__a][__b] << std::
              endl;
            std::cerr << "Q_hat Rel Difference=  " << std::
              abs ((Q_hat[__a][__b] -
                    Q_hat_reference[__a][__b]) /
                   (Q_hat_reference[__a][__b])) << std::endl;
            std::cerr << "Q_hat Abs Difference=  " << std::
              abs (Q_hat[__a][__b] - Q_hat_reference[__a][__b]) << std::endl;
            return 1;
          }

    }
#endif

//=======================================================
//
//               COMPUTE AVX RESULTS
//
//=======================================================

#ifdef ENABLE_AVX_INSTRUCTION_SET
    {
      typedef         T (&refArray1)[3][16];
      typedef         T (&refArray2)[3][16];
      for (int __a = 0; __a < 3; __a++)
        for (int __b = 0; __b < 16; __b++)
          Q_hat[__a][__b] = Q_hat_original[__a][__b];
      for (int i = 0; i < 16; i += 8)
      {
        refArray1       Sigmak = reinterpret_cast < refArray1 > (Sigma[0][i]);
        refArray2       Q_hatk = reinterpret_cast < refArray2 > (Q_hat[0][i]);
        Penalty_Measure_Gradient < COROTATED_TAG, __m256, float[16],
          int[16] >::Run (Sigmak, Q_hatk);
      }
      for (int __a = 0; __a < 3; __a++)
        for (int __b = 0; __b < 16; __b++)
          if (std::
              abs ((Q_hat[__a][__b] -
                    Q_hat_reference[__a][__b]) / (Q_hat_reference[__a][__b])) >
              1)
          {
            std::cerr << "Mismatch detected in AVX implementation" << std::endl;
            std::cerr << "Variable Q_hat:" << std::endl;
            std::
              cerr << "seed=" << seed << ", __a=" << __a << ", __b=" << __b <<
              std::endl;
            std::cerr << "Q_hat AVX=  " << Q_hat[__a][__b] << std::endl;
            std::
              cerr << "Q_hat Reference=  " << Q_hat_reference[__a][__b] << std::
              endl;
            std::cerr << "Q_hat Rel Difference=  " << std::
              abs ((Q_hat[__a][__b] -
                    Q_hat_reference[__a][__b]) /
                   (Q_hat_reference[__a][__b])) << std::endl;
            std::cerr << "Q_hat Abs Difference=  " << std::
              abs (Q_hat[__a][__b] - Q_hat_reference[__a][__b]) << std::endl;
            return 1;
          }

    }
#endif

//=======================================================
//
//               COMPUTE NEON RESULTS
//
//=======================================================

#ifdef ENABLE_NEON_INSTRUCTION_SET
    {
      typedef         T (&refArray1)[3][16];
      typedef         T (&refArray2)[3][16];
      for (int __a = 0; __a < 3; __a++)
        for (int __b = 0; __b < 16; __b++)
          Q_hat[__a][__b] = Q_hat_original[__a][__b];
      for (int i = 0; i < 16; i += 4)
      {
        refArray1       Sigmak = reinterpret_cast < refArray1 > (Sigma[0][i]);
        refArray2       Q_hatk = reinterpret_cast < refArray2 > (Q_hat[0][i]);
        Penalty_Measure_Gradient < COROTATED_TAG, float32x4_t, float[16],
          int[16] >::Run (Sigmak, Q_hatk);
      }
      for (int __a = 0; __a < 3; __a++)
        for (int __b = 0; __b < 16; __b++)
          if (std::
              abs ((Q_hat[__a][__b] -
                    Q_hat_reference[__a][__b]) / (Q_hat_reference[__a][__b])) >
              1)
          {
            std::cerr << "Mismatch detected in NEON implementation" << std::
              endl;
            std::cerr << "Variable Q_hat:" << std::endl;
            std::
              cerr << "seed=" << seed << ", __a=" << __a << ", __b=" << __b <<
              std::endl;
            std::cerr << "Q_hat NEON=  " << Q_hat[__a][__b] << std::endl;
            std::
              cerr << "Q_hat Reference=  " << Q_hat_reference[__a][__b] << std::
              endl;
            std::cerr << "Q_hat Rel Difference=  " << std::
              abs ((Q_hat[__a][__b] -
                    Q_hat_reference[__a][__b]) /
                   (Q_hat_reference[__a][__b])) << std::endl;
            std::cerr << "Q_hat Abs Difference=  " << std::
              abs (Q_hat[__a][__b] - Q_hat_reference[__a][__b]) << std::endl;
            return 1;
          }

    }
#endif

//=======================================================
//
//               COMPUTE MIC RESULTS
//
//=======================================================

#ifdef ENABLE_MIC_INSTRUCTION_SET
    {
      typedef         T (&refArray1)[3][16];
      typedef         T (&refArray2)[3][16];
      for (int __a = 0; __a < 3; __a++)
        for (int __b = 0; __b < 16; __b++)
          Q_hat[__a][__b] = Q_hat_original[__a][__b];
      for (int i = 0; i < 16; i += 16)
      {
        refArray1       Sigmak = reinterpret_cast < refArray1 > (Sigma[0][i]);
        refArray2       Q_hatk = reinterpret_cast < refArray2 > (Q_hat[0][i]);
        Penalty_Measure_Gradient < COROTATED_TAG, __m512, float[16],
          int[16] >::Run (Sigmak, Q_hatk);
      }
      for (int __a = 0; __a < 3; __a++)
        for (int __b = 0; __b < 16; __b++)
          if (std::
              abs ((Q_hat[__a][__b] -
                    Q_hat_reference[__a][__b]) / (Q_hat_reference[__a][__b])) >
              1)
          {
            std::cerr << "Mismatch detected in MIC implementation" << std::endl;
            std::cerr << "Variable Q_hat:" << std::endl;
            std::
              cerr << "seed=" << seed << ", __a=" << __a << ", __b=" << __b <<
              std::endl;
            std::cerr << "Q_hat MIC=  " << Q_hat[__a][__b] << std::endl;
            std::
              cerr << "Q_hat Reference=  " << Q_hat_reference[__a][__b] << std::
              endl;
            std::cerr << "Q_hat Rel Difference=  " << std::
              abs ((Q_hat[__a][__b] -
                    Q_hat_reference[__a][__b]) /
                   (Q_hat_reference[__a][__b])) << std::endl;
            std::cerr << "Q_hat Abs Difference=  " << std::
              abs (Q_hat[__a][__b] - Q_hat_reference[__a][__b]) << std::endl;
            return 1;
          }

    }
#endif

  }



  std::cout << "SIMD check successful!" << std::endl;

  return 0;

}

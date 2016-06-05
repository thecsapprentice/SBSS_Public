
#include <cstdlib>
#include <iostream>
#include "KernelCommon.h"

struct NEOHOOKEAN_TAG;
struct COROTATED_TAG;

#include "Stress_Tensor_Differential.h"
#include "Stress_Tensor_Differential_Reference.h"

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
    std::cout << "Running SIMD Test for Stress_Tensor_Differential " << std::
      endl;


//=======================================================
//
//        DEFINE ALL VARIABLES USED BY KERNEL
//
//=======================================================

    T               dP_hat[9][16] __attribute__ ((aligned (64)));
    T               dP_hat_reference[9][16] __attribute__ ((aligned (64)));
    T               dP_hat_original[9][16] __attribute__ ((aligned (64)));
    T               dPdF[12][16] __attribute__ ((aligned (64)));
    T               dF_hat[9][16] __attribute__ ((aligned (64)));
    T               Q_hat[3][16] __attribute__ ((aligned (64)));
    T               dp[16] __attribute__ ((aligned (64)));
    T               alpha[16] __attribute__ ((aligned (64)));


    for (int __a = 0; __a < 9; __a++)
      for (int __b = 0; __b < 16; __b++)
      {
        dP_hat_original[__a][__b] = Get_Random < float >();
        dP_hat[__a][__b] = dP_hat_original[__a][__b];
        dP_hat_reference[__a][__b] = dP_hat_original[__a][__b];
      }
    for (int __a = 0; __a < 12; __a++)
      for (int __b = 0; __b < 16; __b++)
        dPdF[__a][__b] = Get_Random < float >();
    for (int __a = 0; __a < 9; __a++)
      for (int __b = 0; __b < 16; __b++)
        dF_hat[__a][__b] = Get_Random < float >();
    for (int __a = 0; __a < 3; __a++)
      for (int __b = 0; __b < 16; __b++)
        Q_hat[__a][__b] = Get_Random < float >();
    for (int __a = 0; __a < 16; __a++)
      dp[__a] = Get_Random < float >();
    for (int __a = 0; __a < 16; __a++)
      alpha[__a] = Get_Random < float >();

//=======================================================
//
//             COMPUTE REFERENCE RESULTS
//
//=======================================================

    T               __mdP_hat[9] __attribute__ ((aligned (4)));
    T               __mdP_hat_reference[9] __attribute__ ((aligned (4)));
    T               __mdP_hat_original[9] __attribute__ ((aligned (4)));
    T               __mdPdF[12] __attribute__ ((aligned (4)));
    T               __mdF_hat[9] __attribute__ ((aligned (4)));
    T               __mQ_hat[3] __attribute__ ((aligned (4)));
    T __mdp __attribute__ ((aligned (4)));
    T __malpha __attribute__ ((aligned (4)));
    for (int k = 0; k < 16; k++)
    {
      for (int __a = 0; __a < 9; __a++)
        __mdP_hat_reference[__a] = dP_hat_reference[__a][k];
      for (int __a = 0; __a < 12; __a++)
        __mdPdF[__a] = dPdF[__a][k];
      for (int __a = 0; __a < 9; __a++)
        __mdF_hat[__a] = dF_hat[__a][k];
      for (int __a = 0; __a < 3; __a++)
        __mQ_hat[__a] = Q_hat[__a][k];
      __mdp = dp[k];
      __malpha = alpha[k];
      Stress_Tensor_Differential < float, float, int >(__mdP_hat_reference,
                                                       __mdPdF, __mdF_hat,
                                                       __mQ_hat, __mdp,
                                                       __malpha);
      for (int __a = 0; __a < 9; __a++)
        dP_hat_reference[__a][k] = __mdP_hat_reference[__a];
    }

//=======================================================
//
//               COMPUTE SCALAR RESULTS
//
//=======================================================

    {
      typedef         T (&refArray1)[9][16];
      typedef         T (&refArray2)[12][16];
      typedef         T (&refArray3)[9][16];
      typedef         T (&refArray4)[3][16];
      typedef         T (&refArray5)[16];
      typedef         T (&refArray6)[16];
      for (int __a = 0; __a < 9; __a++)
        for (int __b = 0; __b < 16; __b++)
          dP_hat[__a][__b] = dP_hat_original[__a][__b];
      for (int i = 0; i < 16; i += 1)
      {
        refArray1       dP_hatk = reinterpret_cast < refArray1 > (dP_hat[0][i]);
        refArray2       dPdFk = reinterpret_cast < refArray2 > (dPdF[0][i]);
        refArray3       dF_hatk = reinterpret_cast < refArray3 > (dF_hat[0][i]);
        refArray4       Q_hatk = reinterpret_cast < refArray4 > (Q_hat[0][i]);
        refArray5       dpk = reinterpret_cast < refArray5 > (dp[i]);
        refArray6       alphak = reinterpret_cast < refArray6 > (alpha[i]);
        Stress_Tensor_Differential < float, float[16], int[16] > (dP_hatk,
                                                                  dPdFk,
                                                                  dF_hatk,
                                                                  Q_hatk, dpk,
                                                                  alphak);
      }
      for (int __a = 0; __a < 9; __a++)
        for (int __b = 0; __b < 16; __b++)
          if (std::
              abs ((dP_hat[__a][__b] -
                    dP_hat_reference[__a][__b]) /
                   (dP_hat_reference[__a][__b])) > 1)
          {
            std::cerr << "Mismatch detected in SCALAR implementation" << std::
              endl;
            std::cerr << "Variable dP_hat:" << std::endl;
            std::
              cerr << "seed=" << seed << ", __a=" << __a << ", __b=" << __b <<
              std::endl;
            std::cerr << "dP_hat SCALAR=  " << dP_hat[__a][__b] << std::endl;
            std::
              cerr << "dP_hat Reference=  " << dP_hat_reference[__a][__b] <<
              std::endl;
            std::cerr << "dP_hat Rel Difference=  " << std::
              abs ((dP_hat[__a][__b] -
                    dP_hat_reference[__a][__b]) /
                   (dP_hat_reference[__a][__b])) << std::endl;
            std::cerr << "dP_hat Abs Difference=  " << std::
              abs (dP_hat[__a][__b] - dP_hat_reference[__a][__b]) << std::endl;
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
      typedef         T (&refArray1)[9][16];
      typedef         T (&refArray2)[12][16];
      typedef         T (&refArray3)[9][16];
      typedef         T (&refArray4)[3][16];
      typedef         T (&refArray5)[16];
      typedef         T (&refArray6)[16];
      for (int __a = 0; __a < 9; __a++)
        for (int __b = 0; __b < 16; __b++)
          dP_hat[__a][__b] = dP_hat_original[__a][__b];
      for (int i = 0; i < 16; i += 4)
      {
        refArray1       dP_hatk = reinterpret_cast < refArray1 > (dP_hat[0][i]);
        refArray2       dPdFk = reinterpret_cast < refArray2 > (dPdF[0][i]);
        refArray3       dF_hatk = reinterpret_cast < refArray3 > (dF_hat[0][i]);
        refArray4       Q_hatk = reinterpret_cast < refArray4 > (Q_hat[0][i]);
        refArray5       dpk = reinterpret_cast < refArray5 > (dp[i]);
        refArray6       alphak = reinterpret_cast < refArray6 > (alpha[i]);
        Stress_Tensor_Differential < __m128, float[16], int[16] > (dP_hatk,
                                                                   dPdFk,
                                                                   dF_hatk,
                                                                   Q_hatk, dpk,
                                                                   alphak);
      }
      for (int __a = 0; __a < 9; __a++)
        for (int __b = 0; __b < 16; __b++)
          if (std::
              abs ((dP_hat[__a][__b] -
                    dP_hat_reference[__a][__b]) /
                   (dP_hat_reference[__a][__b])) > 1)
          {
            std::cerr << "Mismatch detected in SSE implementation" << std::endl;
            std::cerr << "Variable dP_hat:" << std::endl;
            std::
              cerr << "seed=" << seed << ", __a=" << __a << ", __b=" << __b <<
              std::endl;
            std::cerr << "dP_hat SSE=  " << dP_hat[__a][__b] << std::endl;
            std::
              cerr << "dP_hat Reference=  " << dP_hat_reference[__a][__b] <<
              std::endl;
            std::cerr << "dP_hat Rel Difference=  " << std::
              abs ((dP_hat[__a][__b] -
                    dP_hat_reference[__a][__b]) /
                   (dP_hat_reference[__a][__b])) << std::endl;
            std::cerr << "dP_hat Abs Difference=  " << std::
              abs (dP_hat[__a][__b] - dP_hat_reference[__a][__b]) << std::endl;
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
      typedef         T (&refArray1)[9][16];
      typedef         T (&refArray2)[12][16];
      typedef         T (&refArray3)[9][16];
      typedef         T (&refArray4)[3][16];
      typedef         T (&refArray5)[16];
      typedef         T (&refArray6)[16];
      for (int __a = 0; __a < 9; __a++)
        for (int __b = 0; __b < 16; __b++)
          dP_hat[__a][__b] = dP_hat_original[__a][__b];
      for (int i = 0; i < 16; i += 8)
      {
        refArray1       dP_hatk = reinterpret_cast < refArray1 > (dP_hat[0][i]);
        refArray2       dPdFk = reinterpret_cast < refArray2 > (dPdF[0][i]);
        refArray3       dF_hatk = reinterpret_cast < refArray3 > (dF_hat[0][i]);
        refArray4       Q_hatk = reinterpret_cast < refArray4 > (Q_hat[0][i]);
        refArray5       dpk = reinterpret_cast < refArray5 > (dp[i]);
        refArray6       alphak = reinterpret_cast < refArray6 > (alpha[i]);
        Stress_Tensor_Differential < __m256, float[16], int[16] > (dP_hatk,
                                                                   dPdFk,
                                                                   dF_hatk,
                                                                   Q_hatk, dpk,
                                                                   alphak);
      }
      for (int __a = 0; __a < 9; __a++)
        for (int __b = 0; __b < 16; __b++)
          if (std::
              abs ((dP_hat[__a][__b] -
                    dP_hat_reference[__a][__b]) /
                   (dP_hat_reference[__a][__b])) > 1)
          {
            std::cerr << "Mismatch detected in AVX implementation" << std::endl;
            std::cerr << "Variable dP_hat:" << std::endl;
            std::
              cerr << "seed=" << seed << ", __a=" << __a << ", __b=" << __b <<
              std::endl;
            std::cerr << "dP_hat AVX=  " << dP_hat[__a][__b] << std::endl;
            std::
              cerr << "dP_hat Reference=  " << dP_hat_reference[__a][__b] <<
              std::endl;
            std::cerr << "dP_hat Rel Difference=  " << std::
              abs ((dP_hat[__a][__b] -
                    dP_hat_reference[__a][__b]) /
                   (dP_hat_reference[__a][__b])) << std::endl;
            std::cerr << "dP_hat Abs Difference=  " << std::
              abs (dP_hat[__a][__b] - dP_hat_reference[__a][__b]) << std::endl;
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
      typedef         T (&refArray1)[9][16];
      typedef         T (&refArray2)[12][16];
      typedef         T (&refArray3)[9][16];
      typedef         T (&refArray4)[3][16];
      typedef         T (&refArray5)[16];
      typedef         T (&refArray6)[16];
      for (int __a = 0; __a < 9; __a++)
        for (int __b = 0; __b < 16; __b++)
          dP_hat[__a][__b] = dP_hat_original[__a][__b];
      for (int i = 0; i < 16; i += 4)
      {
        refArray1       dP_hatk = reinterpret_cast < refArray1 > (dP_hat[0][i]);
        refArray2       dPdFk = reinterpret_cast < refArray2 > (dPdF[0][i]);
        refArray3       dF_hatk = reinterpret_cast < refArray3 > (dF_hat[0][i]);
        refArray4       Q_hatk = reinterpret_cast < refArray4 > (Q_hat[0][i]);
        refArray5       dpk = reinterpret_cast < refArray5 > (dp[i]);
        refArray6       alphak = reinterpret_cast < refArray6 > (alpha[i]);
        Stress_Tensor_Differential < float32x4_t, float[16], int[16] > (dP_hatk,
                                                                        dPdFk,
                                                                        dF_hatk,
                                                                        Q_hatk,
                                                                        dpk,
                                                                        alphak);
      }
      for (int __a = 0; __a < 9; __a++)
        for (int __b = 0; __b < 16; __b++)
          if (std::
              abs ((dP_hat[__a][__b] -
                    dP_hat_reference[__a][__b]) /
                   (dP_hat_reference[__a][__b])) > 1)
          {
            std::cerr << "Mismatch detected in NEON implementation" << std::
              endl;
            std::cerr << "Variable dP_hat:" << std::endl;
            std::
              cerr << "seed=" << seed << ", __a=" << __a << ", __b=" << __b <<
              std::endl;
            std::cerr << "dP_hat NEON=  " << dP_hat[__a][__b] << std::endl;
            std::
              cerr << "dP_hat Reference=  " << dP_hat_reference[__a][__b] <<
              std::endl;
            std::cerr << "dP_hat Rel Difference=  " << std::
              abs ((dP_hat[__a][__b] -
                    dP_hat_reference[__a][__b]) /
                   (dP_hat_reference[__a][__b])) << std::endl;
            std::cerr << "dP_hat Abs Difference=  " << std::
              abs (dP_hat[__a][__b] - dP_hat_reference[__a][__b]) << std::endl;
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
      typedef         T (&refArray1)[9][16];
      typedef         T (&refArray2)[12][16];
      typedef         T (&refArray3)[9][16];
      typedef         T (&refArray4)[3][16];
      typedef         T (&refArray5)[16];
      typedef         T (&refArray6)[16];
      for (int __a = 0; __a < 9; __a++)
        for (int __b = 0; __b < 16; __b++)
          dP_hat[__a][__b] = dP_hat_original[__a][__b];
      for (int i = 0; i < 16; i += 16)
      {
        refArray1       dP_hatk = reinterpret_cast < refArray1 > (dP_hat[0][i]);
        refArray2       dPdFk = reinterpret_cast < refArray2 > (dPdF[0][i]);
        refArray3       dF_hatk = reinterpret_cast < refArray3 > (dF_hat[0][i]);
        refArray4       Q_hatk = reinterpret_cast < refArray4 > (Q_hat[0][i]);
        refArray5       dpk = reinterpret_cast < refArray5 > (dp[i]);
        refArray6       alphak = reinterpret_cast < refArray6 > (alpha[i]);
        Stress_Tensor_Differential < __m512, float[16], int[16] > (dP_hatk,
                                                                   dPdFk,
                                                                   dF_hatk,
                                                                   Q_hatk, dpk,
                                                                   alphak);
      }
      for (int __a = 0; __a < 9; __a++)
        for (int __b = 0; __b < 16; __b++)
          if (std::
              abs ((dP_hat[__a][__b] -
                    dP_hat_reference[__a][__b]) /
                   (dP_hat_reference[__a][__b])) > 1)
          {
            std::cerr << "Mismatch detected in MIC implementation" << std::endl;
            std::cerr << "Variable dP_hat:" << std::endl;
            std::
              cerr << "seed=" << seed << ", __a=" << __a << ", __b=" << __b <<
              std::endl;
            std::cerr << "dP_hat MIC=  " << dP_hat[__a][__b] << std::endl;
            std::
              cerr << "dP_hat Reference=  " << dP_hat_reference[__a][__b] <<
              std::endl;
            std::cerr << "dP_hat Rel Difference=  " << std::
              abs ((dP_hat[__a][__b] -
                    dP_hat_reference[__a][__b]) /
                   (dP_hat_reference[__a][__b])) << std::endl;
            std::cerr << "dP_hat Abs Difference=  " << std::
              abs (dP_hat[__a][__b] - dP_hat_reference[__a][__b]) << std::endl;
            return 1;
          }

    }
#endif

  }



  std::cout << "SIMD check successful!" << std::endl;

  return 0;

}

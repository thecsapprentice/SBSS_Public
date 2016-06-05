
#include <cstdlib>
#include <iostream>
#include "KernelCommon.h"

struct NEOHOOKEAN_TAG;
struct COROTATED_TAG;

#include "Pressure_Force_Differential.h"
#include "Pressure_Force_Differential_Reference.h"

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
    std::cout << "Running SIMD Test for Pressure_Force_Differential " << std::
      endl;


//=======================================================
//
//        DEFINE ALL VARIABLES USED BY KERNEL
//
//=======================================================

    T               dq[16] __attribute__ ((aligned (64)));
    T               dq_reference[16] __attribute__ ((aligned (64)));
    T               dq_original[16] __attribute__ ((aligned (64)));
    T               Q_hat[3][16] __attribute__ ((aligned (64)));
    T               dF_hat[9][16] __attribute__ ((aligned (64)));
    T               dp[16] __attribute__ ((aligned (64)));
    T               alpha[16] __attribute__ ((aligned (64)));
    T               alpha_squared_over_kappa[16] __attribute__ ((aligned (64)));


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
//             COMPUTE REFERENCE RESULTS
//
//=======================================================

    T __mdq __attribute__ ((aligned (4)));
    T __mdq_reference __attribute__ ((aligned (4)));
    T __mdq_original __attribute__ ((aligned (4)));
    T               __mQ_hat[3] __attribute__ ((aligned (4)));
    T               __mdF_hat[9] __attribute__ ((aligned (4)));
    T __mdp __attribute__ ((aligned (4)));
    T __malpha __attribute__ ((aligned (4)));
    T __malpha_squared_over_kappa __attribute__ ((aligned (4)));
    for (int k = 0; k < 16; k++)
    {
      __mdq_reference = dq_reference[k];
      for (int __a = 0; __a < 3; __a++)
        __mQ_hat[__a] = Q_hat[__a][k];
      for (int __a = 0; __a < 9; __a++)
        __mdF_hat[__a] = dF_hat[__a][k];
      __mdp = dp[k];
      __malpha = alpha[k];
      __malpha_squared_over_kappa = alpha_squared_over_kappa[k];
      Pressure_Force_Differential < float, float, int >(__mdq_reference,
                                                        __mQ_hat, __mdF_hat,
                                                        __mdp, __malpha,
                                                        __malpha_squared_over_kappa);
      dq_reference[k] = __mdq_reference;
    }

//=======================================================
//
//               COMPUTE SCALAR RESULTS
//
//=======================================================

    {
      typedef         T (&refArray1)[16];
      typedef         T (&refArray2)[3][16];
      typedef         T (&refArray3)[9][16];
      typedef         T (&refArray4)[16];
      typedef         T (&refArray5)[16];
      typedef         T (&refArray6)[16];
      for (int __a = 0; __a < 16; __a++)
        dq[__a] = dq_original[__a];
      for (int i = 0; i < 16; i += 1)
      {
        refArray1       dqk = reinterpret_cast < refArray1 > (dq[i]);
        refArray2       Q_hatk = reinterpret_cast < refArray2 > (Q_hat[0][i]);
        refArray3       dF_hatk = reinterpret_cast < refArray3 > (dF_hat[0][i]);
        refArray4       dpk = reinterpret_cast < refArray4 > (dp[i]);
        refArray5       alphak = reinterpret_cast < refArray5 > (alpha[i]);
        refArray6       alpha_squared_over_kappak =
          reinterpret_cast < refArray6 > (alpha_squared_over_kappa[i]);
        Pressure_Force_Differential < float, float[16], int[16] > (dqk, Q_hatk,
                                                                   dF_hatk, dpk,
                                                                   alphak,
                                                                   alpha_squared_over_kappak);
      }
      for (int __a = 0; __a < 16; __a++)
        if (std::abs ((dq[__a] - dq_reference[__a]) / (dq_reference[__a])) > 1)
        {
          std::cerr << "Mismatch detected in SCALAR implementation" << std::
            endl;
          std::cerr << "Variable dq:" << std::endl;
          std::cerr << "seed=" << seed << ", __a=" << __a << std::endl;
          std::cerr << "dq SCALAR=  " << dq[__a] << std::endl;
          std::cerr << "dq Reference=  " << dq_reference[__a] << std::endl;
          std::cerr << "dq Rel Difference=  " << std::
            abs ((dq[__a] -
                  dq_reference[__a]) / (dq_reference[__a])) << std::endl;
          std::cerr << "dq Abs Difference=  " << std::abs (dq[__a] -
                                                           dq_reference[__a]) <<
            std::endl;
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
      typedef         T (&refArray1)[16];
      typedef         T (&refArray2)[3][16];
      typedef         T (&refArray3)[9][16];
      typedef         T (&refArray4)[16];
      typedef         T (&refArray5)[16];
      typedef         T (&refArray6)[16];
      for (int __a = 0; __a < 16; __a++)
        dq[__a] = dq_original[__a];
      for (int i = 0; i < 16; i += 4)
      {
        refArray1       dqk = reinterpret_cast < refArray1 > (dq[i]);
        refArray2       Q_hatk = reinterpret_cast < refArray2 > (Q_hat[0][i]);
        refArray3       dF_hatk = reinterpret_cast < refArray3 > (dF_hat[0][i]);
        refArray4       dpk = reinterpret_cast < refArray4 > (dp[i]);
        refArray5       alphak = reinterpret_cast < refArray5 > (alpha[i]);
        refArray6       alpha_squared_over_kappak =
          reinterpret_cast < refArray6 > (alpha_squared_over_kappa[i]);
        Pressure_Force_Differential < __m128, float[16], int[16] > (dqk, Q_hatk,
                                                                    dF_hatk,
                                                                    dpk, alphak,
                                                                    alpha_squared_over_kappak);
      }
      for (int __a = 0; __a < 16; __a++)
        if (std::abs ((dq[__a] - dq_reference[__a]) / (dq_reference[__a])) > 1)
        {
          std::cerr << "Mismatch detected in SSE implementation" << std::endl;
          std::cerr << "Variable dq:" << std::endl;
          std::cerr << "seed=" << seed << ", __a=" << __a << std::endl;
          std::cerr << "dq SSE=  " << dq[__a] << std::endl;
          std::cerr << "dq Reference=  " << dq_reference[__a] << std::endl;
          std::cerr << "dq Rel Difference=  " << std::
            abs ((dq[__a] -
                  dq_reference[__a]) / (dq_reference[__a])) << std::endl;
          std::cerr << "dq Abs Difference=  " << std::abs (dq[__a] -
                                                           dq_reference[__a]) <<
            std::endl;
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
      typedef         T (&refArray1)[16];
      typedef         T (&refArray2)[3][16];
      typedef         T (&refArray3)[9][16];
      typedef         T (&refArray4)[16];
      typedef         T (&refArray5)[16];
      typedef         T (&refArray6)[16];
      for (int __a = 0; __a < 16; __a++)
        dq[__a] = dq_original[__a];
      for (int i = 0; i < 16; i += 8)
      {
        refArray1       dqk = reinterpret_cast < refArray1 > (dq[i]);
        refArray2       Q_hatk = reinterpret_cast < refArray2 > (Q_hat[0][i]);
        refArray3       dF_hatk = reinterpret_cast < refArray3 > (dF_hat[0][i]);
        refArray4       dpk = reinterpret_cast < refArray4 > (dp[i]);
        refArray5       alphak = reinterpret_cast < refArray5 > (alpha[i]);
        refArray6       alpha_squared_over_kappak =
          reinterpret_cast < refArray6 > (alpha_squared_over_kappa[i]);
        Pressure_Force_Differential < __m256, float[16], int[16] > (dqk, Q_hatk,
                                                                    dF_hatk,
                                                                    dpk, alphak,
                                                                    alpha_squared_over_kappak);
      }
      for (int __a = 0; __a < 16; __a++)
        if (std::abs ((dq[__a] - dq_reference[__a]) / (dq_reference[__a])) > 1)
        {
          std::cerr << "Mismatch detected in AVX implementation" << std::endl;
          std::cerr << "Variable dq:" << std::endl;
          std::cerr << "seed=" << seed << ", __a=" << __a << std::endl;
          std::cerr << "dq AVX=  " << dq[__a] << std::endl;
          std::cerr << "dq Reference=  " << dq_reference[__a] << std::endl;
          std::cerr << "dq Rel Difference=  " << std::
            abs ((dq[__a] -
                  dq_reference[__a]) / (dq_reference[__a])) << std::endl;
          std::cerr << "dq Abs Difference=  " << std::abs (dq[__a] -
                                                           dq_reference[__a]) <<
            std::endl;
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
      typedef         T (&refArray1)[16];
      typedef         T (&refArray2)[3][16];
      typedef         T (&refArray3)[9][16];
      typedef         T (&refArray4)[16];
      typedef         T (&refArray5)[16];
      typedef         T (&refArray6)[16];
      for (int __a = 0; __a < 16; __a++)
        dq[__a] = dq_original[__a];
      for (int i = 0; i < 16; i += 4)
      {
        refArray1       dqk = reinterpret_cast < refArray1 > (dq[i]);
        refArray2       Q_hatk = reinterpret_cast < refArray2 > (Q_hat[0][i]);
        refArray3       dF_hatk = reinterpret_cast < refArray3 > (dF_hat[0][i]);
        refArray4       dpk = reinterpret_cast < refArray4 > (dp[i]);
        refArray5       alphak = reinterpret_cast < refArray5 > (alpha[i]);
        refArray6       alpha_squared_over_kappak =
          reinterpret_cast < refArray6 > (alpha_squared_over_kappa[i]);
        Pressure_Force_Differential < float32x4_t, float[16], int[16] > (dqk,
                                                                         Q_hatk,
                                                                         dF_hatk,
                                                                         dpk,
                                                                         alphak,
                                                                         alpha_squared_over_kappak);
      }
      for (int __a = 0; __a < 16; __a++)
        if (std::abs ((dq[__a] - dq_reference[__a]) / (dq_reference[__a])) > 1)
        {
          std::cerr << "Mismatch detected in NEON implementation" << std::endl;
          std::cerr << "Variable dq:" << std::endl;
          std::cerr << "seed=" << seed << ", __a=" << __a << std::endl;
          std::cerr << "dq NEON=  " << dq[__a] << std::endl;
          std::cerr << "dq Reference=  " << dq_reference[__a] << std::endl;
          std::cerr << "dq Rel Difference=  " << std::
            abs ((dq[__a] -
                  dq_reference[__a]) / (dq_reference[__a])) << std::endl;
          std::cerr << "dq Abs Difference=  " << std::abs (dq[__a] -
                                                           dq_reference[__a]) <<
            std::endl;
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
      typedef         T (&refArray1)[16];
      typedef         T (&refArray2)[3][16];
      typedef         T (&refArray3)[9][16];
      typedef         T (&refArray4)[16];
      typedef         T (&refArray5)[16];
      typedef         T (&refArray6)[16];
      for (int __a = 0; __a < 16; __a++)
        dq[__a] = dq_original[__a];
      for (int i = 0; i < 16; i += 16)
      {
        refArray1       dqk = reinterpret_cast < refArray1 > (dq[i]);
        refArray2       Q_hatk = reinterpret_cast < refArray2 > (Q_hat[0][i]);
        refArray3       dF_hatk = reinterpret_cast < refArray3 > (dF_hat[0][i]);
        refArray4       dpk = reinterpret_cast < refArray4 > (dp[i]);
        refArray5       alphak = reinterpret_cast < refArray5 > (alpha[i]);
        refArray6       alpha_squared_over_kappak =
          reinterpret_cast < refArray6 > (alpha_squared_over_kappa[i]);
        Pressure_Force_Differential < __m512, float[16], int[16] > (dqk, Q_hatk,
                                                                    dF_hatk,
                                                                    dpk, alphak,
                                                                    alpha_squared_over_kappak);
      }
      for (int __a = 0; __a < 16; __a++)
        if (std::abs ((dq[__a] - dq_reference[__a]) / (dq_reference[__a])) > 1)
        {
          std::cerr << "Mismatch detected in MIC implementation" << std::endl;
          std::cerr << "Variable dq:" << std::endl;
          std::cerr << "seed=" << seed << ", __a=" << __a << std::endl;
          std::cerr << "dq MIC=  " << dq[__a] << std::endl;
          std::cerr << "dq Reference=  " << dq_reference[__a] << std::endl;
          std::cerr << "dq Rel Difference=  " << std::
            abs ((dq[__a] -
                  dq_reference[__a]) / (dq_reference[__a])) << std::endl;
          std::cerr << "dq Abs Difference=  " << std::abs (dq[__a] -
                                                           dq_reference[__a]) <<
            std::endl;
          return 1;
        }

    }
#endif

  }



  std::cout << "SIMD check successful!" << std::endl;

  return 0;

}

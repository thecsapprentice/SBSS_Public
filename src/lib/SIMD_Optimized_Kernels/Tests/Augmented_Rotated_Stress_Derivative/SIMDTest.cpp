
#include <cstdlib>
#include <iostream>
#include "KernelCommon.h"

struct NEOHOOKEAN_TAG;
struct COROTATED_TAG;

#include "Augmented_Rotated_Stress_Derivative.h"
#include "Augmented_Rotated_Stress_Derivative_Reference.h"

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
    std::
      cout << "Running SIMD Test for Augmented_Rotated_Stress_Derivative " <<
      std::endl;


//=======================================================
//
//        DEFINE ALL VARIABLES USED BY KERNEL
//
//=======================================================

    T               dPdF[12][16] __attribute__ ((aligned (64)));
    T               dPdF_reference[12][16] __attribute__ ((aligned (64)));
    T               dPdF_original[12][16] __attribute__ ((aligned (64)));
    T               Sigma[3][16] __attribute__ ((aligned (64)));
    T               p[16] __attribute__ ((aligned (64)));
    T               mu[16] __attribute__ ((aligned (64)));
    T               alpha[16] __attribute__ ((aligned (64)));


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
      alpha[__a] = Get_Random < float >();

//=======================================================
//
//             COMPUTE REFERENCE RESULTS
//
//=======================================================

    T               __mdPdF[12] __attribute__ ((aligned (4)));
    T               __mdPdF_reference[12] __attribute__ ((aligned (4)));
    T               __mdPdF_original[12] __attribute__ ((aligned (4)));
    T               __mSigma[3] __attribute__ ((aligned (4)));
    T __mp __attribute__ ((aligned (4)));
    T __mmu __attribute__ ((aligned (4)));
    T __malpha __attribute__ ((aligned (4)));
    for (int k = 0; k < 16; k++)
    {
      for (int __a = 0; __a < 12; __a++)
        __mdPdF_reference[__a] = dPdF_reference[__a][k];
      for (int __a = 0; __a < 3; __a++)
        __mSigma[__a] = Sigma[__a][k];
      __mp = p[k];
      __mmu = mu[k];
      __malpha = alpha[k];
      Augmented_Rotated_Stress_Derivative < NEOHOOKEAN_TAG, float, float,
        int >::Run (__mdPdF_reference, __mSigma, __mp, __mmu, __malpha);
      for (int __a = 0; __a < 12; __a++)
        dPdF_reference[__a][k] = __mdPdF_reference[__a];
    }

//=======================================================
//
//               COMPUTE SCALAR RESULTS
//
//=======================================================

    {
      typedef         T (&refArray1)[12][16];
      typedef         T (&refArray2)[3][16];
      typedef         T (&refArray3)[16];
      typedef         T (&refArray4)[16];
      typedef         T (&refArray5)[16];
      for (int __a = 0; __a < 12; __a++)
        for (int __b = 0; __b < 16; __b++)
          dPdF[__a][__b] = dPdF_original[__a][__b];
      for (int i = 0; i < 16; i += 1)
      {
        refArray1       dPdFk = reinterpret_cast < refArray1 > (dPdF[0][i]);
        refArray2       Sigmak = reinterpret_cast < refArray2 > (Sigma[0][i]);
        refArray3       pk = reinterpret_cast < refArray3 > (p[i]);
        refArray4       muk = reinterpret_cast < refArray4 > (mu[i]);
        refArray5       alphak = reinterpret_cast < refArray5 > (alpha[i]);
        Augmented_Rotated_Stress_Derivative < NEOHOOKEAN_TAG, float, float[16],
          int[16] >::Run (dPdFk, Sigmak, pk, muk, alphak);
      }
      for (int __a = 0; __a < 12; __a++)
        for (int __b = 0; __b < 16; __b++)
          if (std::
              abs ((dPdF[__a][__b] -
                    dPdF_reference[__a][__b]) / (dPdF_reference[__a][__b])) > 1)
          {
            std::cerr << "Mismatch detected in SCALAR implementation" << std::
              endl;
            std::cerr << "Variable dPdF:" << std::endl;
            std::
              cerr << "seed=" << seed << ", __a=" << __a << ", __b=" << __b <<
              std::endl;
            std::cerr << "dPdF SCALAR=  " << dPdF[__a][__b] << std::endl;
            std::
              cerr << "dPdF Reference=  " << dPdF_reference[__a][__b] << std::
              endl;
            std::cerr << "dPdF Rel Difference=  " << std::
              abs ((dPdF[__a][__b] -
                    dPdF_reference[__a][__b]) /
                   (dPdF_reference[__a][__b])) << std::endl;
            std::cerr << "dPdF Abs Difference=  " << std::abs (dPdF[__a][__b] -
                                                               dPdF_reference
                                                               [__a][__b]) <<
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
      typedef         T (&refArray1)[12][16];
      typedef         T (&refArray2)[3][16];
      typedef         T (&refArray3)[16];
      typedef         T (&refArray4)[16];
      typedef         T (&refArray5)[16];
      for (int __a = 0; __a < 12; __a++)
        for (int __b = 0; __b < 16; __b++)
          dPdF[__a][__b] = dPdF_original[__a][__b];
      for (int i = 0; i < 16; i += 4)
      {
        refArray1       dPdFk = reinterpret_cast < refArray1 > (dPdF[0][i]);
        refArray2       Sigmak = reinterpret_cast < refArray2 > (Sigma[0][i]);
        refArray3       pk = reinterpret_cast < refArray3 > (p[i]);
        refArray4       muk = reinterpret_cast < refArray4 > (mu[i]);
        refArray5       alphak = reinterpret_cast < refArray5 > (alpha[i]);
        Augmented_Rotated_Stress_Derivative < NEOHOOKEAN_TAG, __m128, float[16],
          int[16] >::Run (dPdFk, Sigmak, pk, muk, alphak);
      }
      for (int __a = 0; __a < 12; __a++)
        for (int __b = 0; __b < 16; __b++)
          if (std::
              abs ((dPdF[__a][__b] -
                    dPdF_reference[__a][__b]) / (dPdF_reference[__a][__b])) > 1)
          {
            std::cerr << "Mismatch detected in SSE implementation" << std::endl;
            std::cerr << "Variable dPdF:" << std::endl;
            std::
              cerr << "seed=" << seed << ", __a=" << __a << ", __b=" << __b <<
              std::endl;
            std::cerr << "dPdF SSE=  " << dPdF[__a][__b] << std::endl;
            std::
              cerr << "dPdF Reference=  " << dPdF_reference[__a][__b] << std::
              endl;
            std::cerr << "dPdF Rel Difference=  " << std::
              abs ((dPdF[__a][__b] -
                    dPdF_reference[__a][__b]) /
                   (dPdF_reference[__a][__b])) << std::endl;
            std::cerr << "dPdF Abs Difference=  " << std::abs (dPdF[__a][__b] -
                                                               dPdF_reference
                                                               [__a][__b]) <<
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
      typedef         T (&refArray1)[12][16];
      typedef         T (&refArray2)[3][16];
      typedef         T (&refArray3)[16];
      typedef         T (&refArray4)[16];
      typedef         T (&refArray5)[16];
      for (int __a = 0; __a < 12; __a++)
        for (int __b = 0; __b < 16; __b++)
          dPdF[__a][__b] = dPdF_original[__a][__b];
      for (int i = 0; i < 16; i += 8)
      {
        refArray1       dPdFk = reinterpret_cast < refArray1 > (dPdF[0][i]);
        refArray2       Sigmak = reinterpret_cast < refArray2 > (Sigma[0][i]);
        refArray3       pk = reinterpret_cast < refArray3 > (p[i]);
        refArray4       muk = reinterpret_cast < refArray4 > (mu[i]);
        refArray5       alphak = reinterpret_cast < refArray5 > (alpha[i]);
        Augmented_Rotated_Stress_Derivative < NEOHOOKEAN_TAG, __m256, float[16],
          int[16] >::Run (dPdFk, Sigmak, pk, muk, alphak);
      }
      for (int __a = 0; __a < 12; __a++)
        for (int __b = 0; __b < 16; __b++)
          if (std::
              abs ((dPdF[__a][__b] -
                    dPdF_reference[__a][__b]) / (dPdF_reference[__a][__b])) > 1)
          {
            std::cerr << "Mismatch detected in AVX implementation" << std::endl;
            std::cerr << "Variable dPdF:" << std::endl;
            std::
              cerr << "seed=" << seed << ", __a=" << __a << ", __b=" << __b <<
              std::endl;
            std::cerr << "dPdF AVX=  " << dPdF[__a][__b] << std::endl;
            std::
              cerr << "dPdF Reference=  " << dPdF_reference[__a][__b] << std::
              endl;
            std::cerr << "dPdF Rel Difference=  " << std::
              abs ((dPdF[__a][__b] -
                    dPdF_reference[__a][__b]) /
                   (dPdF_reference[__a][__b])) << std::endl;
            std::cerr << "dPdF Abs Difference=  " << std::abs (dPdF[__a][__b] -
                                                               dPdF_reference
                                                               [__a][__b]) <<
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
      typedef         T (&refArray1)[12][16];
      typedef         T (&refArray2)[3][16];
      typedef         T (&refArray3)[16];
      typedef         T (&refArray4)[16];
      typedef         T (&refArray5)[16];
      for (int __a = 0; __a < 12; __a++)
        for (int __b = 0; __b < 16; __b++)
          dPdF[__a][__b] = dPdF_original[__a][__b];
      for (int i = 0; i < 16; i += 4)
      {
        refArray1       dPdFk = reinterpret_cast < refArray1 > (dPdF[0][i]);
        refArray2       Sigmak = reinterpret_cast < refArray2 > (Sigma[0][i]);
        refArray3       pk = reinterpret_cast < refArray3 > (p[i]);
        refArray4       muk = reinterpret_cast < refArray4 > (mu[i]);
        refArray5       alphak = reinterpret_cast < refArray5 > (alpha[i]);
        Augmented_Rotated_Stress_Derivative < NEOHOOKEAN_TAG, float32x4_t,
          float[16], int[16] >::Run (dPdFk, Sigmak, pk, muk, alphak);
      }
      for (int __a = 0; __a < 12; __a++)
        for (int __b = 0; __b < 16; __b++)
          if (std::
              abs ((dPdF[__a][__b] -
                    dPdF_reference[__a][__b]) / (dPdF_reference[__a][__b])) > 1)
          {
            std::cerr << "Mismatch detected in NEON implementation" << std::
              endl;
            std::cerr << "Variable dPdF:" << std::endl;
            std::
              cerr << "seed=" << seed << ", __a=" << __a << ", __b=" << __b <<
              std::endl;
            std::cerr << "dPdF NEON=  " << dPdF[__a][__b] << std::endl;
            std::
              cerr << "dPdF Reference=  " << dPdF_reference[__a][__b] << std::
              endl;
            std::cerr << "dPdF Rel Difference=  " << std::
              abs ((dPdF[__a][__b] -
                    dPdF_reference[__a][__b]) /
                   (dPdF_reference[__a][__b])) << std::endl;
            std::cerr << "dPdF Abs Difference=  " << std::abs (dPdF[__a][__b] -
                                                               dPdF_reference
                                                               [__a][__b]) <<
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
      typedef         T (&refArray1)[12][16];
      typedef         T (&refArray2)[3][16];
      typedef         T (&refArray3)[16];
      typedef         T (&refArray4)[16];
      typedef         T (&refArray5)[16];
      for (int __a = 0; __a < 12; __a++)
        for (int __b = 0; __b < 16; __b++)
          dPdF[__a][__b] = dPdF_original[__a][__b];
      for (int i = 0; i < 16; i += 16)
      {
        refArray1       dPdFk = reinterpret_cast < refArray1 > (dPdF[0][i]);
        refArray2       Sigmak = reinterpret_cast < refArray2 > (Sigma[0][i]);
        refArray3       pk = reinterpret_cast < refArray3 > (p[i]);
        refArray4       muk = reinterpret_cast < refArray4 > (mu[i]);
        refArray5       alphak = reinterpret_cast < refArray5 > (alpha[i]);
        Augmented_Rotated_Stress_Derivative < NEOHOOKEAN_TAG, __m512, float[16],
          int[16] >::Run (dPdFk, Sigmak, pk, muk, alphak);
      }
      for (int __a = 0; __a < 12; __a++)
        for (int __b = 0; __b < 16; __b++)
          if (std::
              abs ((dPdF[__a][__b] -
                    dPdF_reference[__a][__b]) / (dPdF_reference[__a][__b])) > 1)
          {
            std::cerr << "Mismatch detected in MIC implementation" << std::endl;
            std::cerr << "Variable dPdF:" << std::endl;
            std::
              cerr << "seed=" << seed << ", __a=" << __a << ", __b=" << __b <<
              std::endl;
            std::cerr << "dPdF MIC=  " << dPdF[__a][__b] << std::endl;
            std::
              cerr << "dPdF Reference=  " << dPdF_reference[__a][__b] << std::
              endl;
            std::cerr << "dPdF Rel Difference=  " << std::
              abs ((dPdF[__a][__b] -
                    dPdF_reference[__a][__b]) /
                   (dPdF_reference[__a][__b])) << std::endl;
            std::cerr << "dPdF Abs Difference=  " << std::abs (dPdF[__a][__b] -
                                                               dPdF_reference
                                                               [__a][__b]) <<
              std::endl;
            return 1;
          }

    }
#endif

  }



  {
    std::
      cout << "Running SIMD Test for Augmented_Rotated_Stress_Derivative " <<
      std::endl;


//=======================================================
//
//        DEFINE ALL VARIABLES USED BY KERNEL
//
//=======================================================

    T               dPdF[12][16] __attribute__ ((aligned (64)));
    T               dPdF_reference[12][16] __attribute__ ((aligned (64)));
    T               dPdF_original[12][16] __attribute__ ((aligned (64)));
    T               Sigma[3][16] __attribute__ ((aligned (64)));
    T               p[16] __attribute__ ((aligned (64)));
    T               mu[16] __attribute__ ((aligned (64)));
    T               alpha[16] __attribute__ ((aligned (64)));


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
      alpha[__a] = Get_Random < float >();

//=======================================================
//
//             COMPUTE REFERENCE RESULTS
//
//=======================================================

    T               __mdPdF[12] __attribute__ ((aligned (4)));
    T               __mdPdF_reference[12] __attribute__ ((aligned (4)));
    T               __mdPdF_original[12] __attribute__ ((aligned (4)));
    T               __mSigma[3] __attribute__ ((aligned (4)));
    T __mp __attribute__ ((aligned (4)));
    T __mmu __attribute__ ((aligned (4)));
    T __malpha __attribute__ ((aligned (4)));
    for (int k = 0; k < 16; k++)
    {
      for (int __a = 0; __a < 12; __a++)
        __mdPdF_reference[__a] = dPdF_reference[__a][k];
      for (int __a = 0; __a < 3; __a++)
        __mSigma[__a] = Sigma[__a][k];
      __mp = p[k];
      __mmu = mu[k];
      __malpha = alpha[k];
      Augmented_Rotated_Stress_Derivative < COROTATED_TAG, float, float,
        int >::Run (__mdPdF_reference, __mSigma, __mp, __mmu, __malpha);
      for (int __a = 0; __a < 12; __a++)
        dPdF_reference[__a][k] = __mdPdF_reference[__a];
    }

//=======================================================
//
//               COMPUTE SCALAR RESULTS
//
//=======================================================

    {
      typedef         T (&refArray1)[12][16];
      typedef         T (&refArray2)[3][16];
      typedef         T (&refArray3)[16];
      typedef         T (&refArray4)[16];
      typedef         T (&refArray5)[16];
      for (int __a = 0; __a < 12; __a++)
        for (int __b = 0; __b < 16; __b++)
          dPdF[__a][__b] = dPdF_original[__a][__b];
      for (int i = 0; i < 16; i += 1)
      {
        refArray1       dPdFk = reinterpret_cast < refArray1 > (dPdF[0][i]);
        refArray2       Sigmak = reinterpret_cast < refArray2 > (Sigma[0][i]);
        refArray3       pk = reinterpret_cast < refArray3 > (p[i]);
        refArray4       muk = reinterpret_cast < refArray4 > (mu[i]);
        refArray5       alphak = reinterpret_cast < refArray5 > (alpha[i]);
        Augmented_Rotated_Stress_Derivative < COROTATED_TAG, float, float[16],
          int[16] >::Run (dPdFk, Sigmak, pk, muk, alphak);
      }
      for (int __a = 0; __a < 12; __a++)
        for (int __b = 0; __b < 16; __b++)
          if (std::
              abs ((dPdF[__a][__b] -
                    dPdF_reference[__a][__b]) / (dPdF_reference[__a][__b])) > 1)
          {
            std::cerr << "Mismatch detected in SCALAR implementation" << std::
              endl;
            std::cerr << "Variable dPdF:" << std::endl;
            std::
              cerr << "seed=" << seed << ", __a=" << __a << ", __b=" << __b <<
              std::endl;
            std::cerr << "dPdF SCALAR=  " << dPdF[__a][__b] << std::endl;
            std::
              cerr << "dPdF Reference=  " << dPdF_reference[__a][__b] << std::
              endl;
            std::cerr << "dPdF Rel Difference=  " << std::
              abs ((dPdF[__a][__b] -
                    dPdF_reference[__a][__b]) /
                   (dPdF_reference[__a][__b])) << std::endl;
            std::cerr << "dPdF Abs Difference=  " << std::abs (dPdF[__a][__b] -
                                                               dPdF_reference
                                                               [__a][__b]) <<
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
      typedef         T (&refArray1)[12][16];
      typedef         T (&refArray2)[3][16];
      typedef         T (&refArray3)[16];
      typedef         T (&refArray4)[16];
      typedef         T (&refArray5)[16];
      for (int __a = 0; __a < 12; __a++)
        for (int __b = 0; __b < 16; __b++)
          dPdF[__a][__b] = dPdF_original[__a][__b];
      for (int i = 0; i < 16; i += 4)
      {
        refArray1       dPdFk = reinterpret_cast < refArray1 > (dPdF[0][i]);
        refArray2       Sigmak = reinterpret_cast < refArray2 > (Sigma[0][i]);
        refArray3       pk = reinterpret_cast < refArray3 > (p[i]);
        refArray4       muk = reinterpret_cast < refArray4 > (mu[i]);
        refArray5       alphak = reinterpret_cast < refArray5 > (alpha[i]);
        Augmented_Rotated_Stress_Derivative < COROTATED_TAG, __m128, float[16],
          int[16] >::Run (dPdFk, Sigmak, pk, muk, alphak);
      }
      for (int __a = 0; __a < 12; __a++)
        for (int __b = 0; __b < 16; __b++)
          if (std::
              abs ((dPdF[__a][__b] -
                    dPdF_reference[__a][__b]) / (dPdF_reference[__a][__b])) > 1)
          {
            std::cerr << "Mismatch detected in SSE implementation" << std::endl;
            std::cerr << "Variable dPdF:" << std::endl;
            std::
              cerr << "seed=" << seed << ", __a=" << __a << ", __b=" << __b <<
              std::endl;
            std::cerr << "dPdF SSE=  " << dPdF[__a][__b] << std::endl;
            std::
              cerr << "dPdF Reference=  " << dPdF_reference[__a][__b] << std::
              endl;
            std::cerr << "dPdF Rel Difference=  " << std::
              abs ((dPdF[__a][__b] -
                    dPdF_reference[__a][__b]) /
                   (dPdF_reference[__a][__b])) << std::endl;
            std::cerr << "dPdF Abs Difference=  " << std::abs (dPdF[__a][__b] -
                                                               dPdF_reference
                                                               [__a][__b]) <<
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
      typedef         T (&refArray1)[12][16];
      typedef         T (&refArray2)[3][16];
      typedef         T (&refArray3)[16];
      typedef         T (&refArray4)[16];
      typedef         T (&refArray5)[16];
      for (int __a = 0; __a < 12; __a++)
        for (int __b = 0; __b < 16; __b++)
          dPdF[__a][__b] = dPdF_original[__a][__b];
      for (int i = 0; i < 16; i += 8)
      {
        refArray1       dPdFk = reinterpret_cast < refArray1 > (dPdF[0][i]);
        refArray2       Sigmak = reinterpret_cast < refArray2 > (Sigma[0][i]);
        refArray3       pk = reinterpret_cast < refArray3 > (p[i]);
        refArray4       muk = reinterpret_cast < refArray4 > (mu[i]);
        refArray5       alphak = reinterpret_cast < refArray5 > (alpha[i]);
        Augmented_Rotated_Stress_Derivative < COROTATED_TAG, __m256, float[16],
          int[16] >::Run (dPdFk, Sigmak, pk, muk, alphak);
      }
      for (int __a = 0; __a < 12; __a++)
        for (int __b = 0; __b < 16; __b++)
          if (std::
              abs ((dPdF[__a][__b] -
                    dPdF_reference[__a][__b]) / (dPdF_reference[__a][__b])) > 1)
          {
            std::cerr << "Mismatch detected in AVX implementation" << std::endl;
            std::cerr << "Variable dPdF:" << std::endl;
            std::
              cerr << "seed=" << seed << ", __a=" << __a << ", __b=" << __b <<
              std::endl;
            std::cerr << "dPdF AVX=  " << dPdF[__a][__b] << std::endl;
            std::
              cerr << "dPdF Reference=  " << dPdF_reference[__a][__b] << std::
              endl;
            std::cerr << "dPdF Rel Difference=  " << std::
              abs ((dPdF[__a][__b] -
                    dPdF_reference[__a][__b]) /
                   (dPdF_reference[__a][__b])) << std::endl;
            std::cerr << "dPdF Abs Difference=  " << std::abs (dPdF[__a][__b] -
                                                               dPdF_reference
                                                               [__a][__b]) <<
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
      typedef         T (&refArray1)[12][16];
      typedef         T (&refArray2)[3][16];
      typedef         T (&refArray3)[16];
      typedef         T (&refArray4)[16];
      typedef         T (&refArray5)[16];
      for (int __a = 0; __a < 12; __a++)
        for (int __b = 0; __b < 16; __b++)
          dPdF[__a][__b] = dPdF_original[__a][__b];
      for (int i = 0; i < 16; i += 4)
      {
        refArray1       dPdFk = reinterpret_cast < refArray1 > (dPdF[0][i]);
        refArray2       Sigmak = reinterpret_cast < refArray2 > (Sigma[0][i]);
        refArray3       pk = reinterpret_cast < refArray3 > (p[i]);
        refArray4       muk = reinterpret_cast < refArray4 > (mu[i]);
        refArray5       alphak = reinterpret_cast < refArray5 > (alpha[i]);
        Augmented_Rotated_Stress_Derivative < COROTATED_TAG, float32x4_t,
          float[16], int[16] >::Run (dPdFk, Sigmak, pk, muk, alphak);
      }
      for (int __a = 0; __a < 12; __a++)
        for (int __b = 0; __b < 16; __b++)
          if (std::
              abs ((dPdF[__a][__b] -
                    dPdF_reference[__a][__b]) / (dPdF_reference[__a][__b])) > 1)
          {
            std::cerr << "Mismatch detected in NEON implementation" << std::
              endl;
            std::cerr << "Variable dPdF:" << std::endl;
            std::
              cerr << "seed=" << seed << ", __a=" << __a << ", __b=" << __b <<
              std::endl;
            std::cerr << "dPdF NEON=  " << dPdF[__a][__b] << std::endl;
            std::
              cerr << "dPdF Reference=  " << dPdF_reference[__a][__b] << std::
              endl;
            std::cerr << "dPdF Rel Difference=  " << std::
              abs ((dPdF[__a][__b] -
                    dPdF_reference[__a][__b]) /
                   (dPdF_reference[__a][__b])) << std::endl;
            std::cerr << "dPdF Abs Difference=  " << std::abs (dPdF[__a][__b] -
                                                               dPdF_reference
                                                               [__a][__b]) <<
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
      typedef         T (&refArray1)[12][16];
      typedef         T (&refArray2)[3][16];
      typedef         T (&refArray3)[16];
      typedef         T (&refArray4)[16];
      typedef         T (&refArray5)[16];
      for (int __a = 0; __a < 12; __a++)
        for (int __b = 0; __b < 16; __b++)
          dPdF[__a][__b] = dPdF_original[__a][__b];
      for (int i = 0; i < 16; i += 16)
      {
        refArray1       dPdFk = reinterpret_cast < refArray1 > (dPdF[0][i]);
        refArray2       Sigmak = reinterpret_cast < refArray2 > (Sigma[0][i]);
        refArray3       pk = reinterpret_cast < refArray3 > (p[i]);
        refArray4       muk = reinterpret_cast < refArray4 > (mu[i]);
        refArray5       alphak = reinterpret_cast < refArray5 > (alpha[i]);
        Augmented_Rotated_Stress_Derivative < COROTATED_TAG, __m512, float[16],
          int[16] >::Run (dPdFk, Sigmak, pk, muk, alphak);
      }
      for (int __a = 0; __a < 12; __a++)
        for (int __b = 0; __b < 16; __b++)
          if (std::
              abs ((dPdF[__a][__b] -
                    dPdF_reference[__a][__b]) / (dPdF_reference[__a][__b])) > 1)
          {
            std::cerr << "Mismatch detected in MIC implementation" << std::endl;
            std::cerr << "Variable dPdF:" << std::endl;
            std::
              cerr << "seed=" << seed << ", __a=" << __a << ", __b=" << __b <<
              std::endl;
            std::cerr << "dPdF MIC=  " << dPdF[__a][__b] << std::endl;
            std::
              cerr << "dPdF Reference=  " << dPdF_reference[__a][__b] << std::
              endl;
            std::cerr << "dPdF Rel Difference=  " << std::
              abs ((dPdF[__a][__b] -
                    dPdF_reference[__a][__b]) /
                   (dPdF_reference[__a][__b])) << std::endl;
            std::cerr << "dPdF Abs Difference=  " << std::abs (dPdF[__a][__b] -
                                                               dPdF_reference
                                                               [__a][__b]) <<
              std::endl;
            return 1;
          }

    }
#endif

  }



  std::cout << "SIMD check successful!" << std::endl;

  return 0;

}

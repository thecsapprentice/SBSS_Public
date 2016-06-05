
#include <cstdlib>
#include <iostream>
#include "KernelCommon.h"

struct NEOHOOKEAN_TAG;
struct COROTATED_TAG;

#include "Rotated_Stress_Derivative_Contraction.h"
#include "Rotated_Stress_Derivative_Contraction_Reference.h"

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
      cout << "Running SIMD Test for Rotated_Stress_Derivative_Contraction " <<
      std::endl;


//=======================================================
//
//        DEFINE ALL VARIABLES USED BY KERNEL
//
//=======================================================

    T               dPdF[12][16] __attribute__ ((aligned (64)));
    T               dF_Hat[9][16] __attribute__ ((aligned (64)));
    T               dP_Hat[9][16] __attribute__ ((aligned (64)));
    T               dP_Hat_reference[9][16] __attribute__ ((aligned (64)));
    T               dP_Hat_original[9][16] __attribute__ ((aligned (64)));


    for (int __a = 0; __a < 12; __a++)
      for (int __b = 0; __b < 16; __b++)
        dPdF[__a][__b] = Get_Random < float >();
    for (int __a = 0; __a < 9; __a++)
      for (int __b = 0; __b < 16; __b++)
        dF_Hat[__a][__b] = Get_Random < float >();
    for (int __a = 0; __a < 9; __a++)
      for (int __b = 0; __b < 16; __b++)
      {
        dP_Hat_original[__a][__b] = Get_Random < float >();
        dP_Hat[__a][__b] = dP_Hat_original[__a][__b];
        dP_Hat_reference[__a][__b] = dP_Hat_original[__a][__b];
      }


//=======================================================
//
//             COMPUTE REFERENCE RESULTS
//
//=======================================================

    T               __mdPdF[12] __attribute__ ((aligned (4)));
    T               __mdF_Hat[9] __attribute__ ((aligned (4)));
    T               __mdP_Hat[9] __attribute__ ((aligned (4)));
    T               __mdP_Hat_reference[9] __attribute__ ((aligned (4)));
    T               __mdP_Hat_original[9] __attribute__ ((aligned (4)));
    for (int k = 0; k < 16; k++)
    {
      for (int __a = 0; __a < 12; __a++)
        __mdPdF[__a] = dPdF[__a][k];
      for (int __a = 0; __a < 9; __a++)
        __mdF_Hat[__a] = dF_Hat[__a][k];
      for (int __a = 0; __a < 9; __a++)
        __mdP_Hat_reference[__a] = dP_Hat_reference[__a][k];
      Rotated_Stress_Derivative_Contraction < float, float, int >(__mdPdF,
                                                                  __mdF_Hat,
                                                                  __mdP_Hat_reference);
      for (int __a = 0; __a < 9; __a++)
        dP_Hat_reference[__a][k] = __mdP_Hat_reference[__a];
    }

//=======================================================
//
//               COMPUTE SCALAR RESULTS
//
//=======================================================

    {
      typedef         T (&refArray1)[12][16];
      typedef         T (&refArray2)[9][16];
      typedef         T (&refArray3)[9][16];
      for (int __a = 0; __a < 9; __a++)
        for (int __b = 0; __b < 16; __b++)
          dP_Hat[__a][__b] = dP_Hat_original[__a][__b];
      for (int i = 0; i < 16; i += 1)
      {
        refArray1       dPdFk = reinterpret_cast < refArray1 > (dPdF[0][i]);
        refArray2       dF_Hatk = reinterpret_cast < refArray2 > (dF_Hat[0][i]);
        refArray3       dP_Hatk = reinterpret_cast < refArray3 > (dP_Hat[0][i]);
        Rotated_Stress_Derivative_Contraction < float, float[16],
          int[16] > (dPdFk, dF_Hatk, dP_Hatk);
      }
      for (int __a = 0; __a < 9; __a++)
        for (int __b = 0; __b < 16; __b++)
          if (std::
              abs ((dP_Hat[__a][__b] -
                    dP_Hat_reference[__a][__b]) /
                   (dP_Hat_reference[__a][__b])) > 1)
          {
            std::cerr << "Mismatch detected in SCALAR implementation" << std::
              endl;
            std::cerr << "Variable dP_Hat:" << std::endl;
            std::
              cerr << "seed=" << seed << ", __a=" << __a << ", __b=" << __b <<
              std::endl;
            std::cerr << "dP_Hat SCALAR=  " << dP_Hat[__a][__b] << std::endl;
            std::
              cerr << "dP_Hat Reference=  " << dP_Hat_reference[__a][__b] <<
              std::endl;
            std::cerr << "dP_Hat Rel Difference=  " << std::
              abs ((dP_Hat[__a][__b] -
                    dP_Hat_reference[__a][__b]) /
                   (dP_Hat_reference[__a][__b])) << std::endl;
            std::cerr << "dP_Hat Abs Difference=  " << std::
              abs (dP_Hat[__a][__b] - dP_Hat_reference[__a][__b]) << std::endl;
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
      typedef         T (&refArray2)[9][16];
      typedef         T (&refArray3)[9][16];
      for (int __a = 0; __a < 9; __a++)
        for (int __b = 0; __b < 16; __b++)
          dP_Hat[__a][__b] = dP_Hat_original[__a][__b];
      for (int i = 0; i < 16; i += 4)
      {
        refArray1       dPdFk = reinterpret_cast < refArray1 > (dPdF[0][i]);
        refArray2       dF_Hatk = reinterpret_cast < refArray2 > (dF_Hat[0][i]);
        refArray3       dP_Hatk = reinterpret_cast < refArray3 > (dP_Hat[0][i]);
        Rotated_Stress_Derivative_Contraction < __m128, float[16],
          int[16] > (dPdFk, dF_Hatk, dP_Hatk);
      }
      for (int __a = 0; __a < 9; __a++)
        for (int __b = 0; __b < 16; __b++)
          if (std::
              abs ((dP_Hat[__a][__b] -
                    dP_Hat_reference[__a][__b]) /
                   (dP_Hat_reference[__a][__b])) > 1)
          {
            std::cerr << "Mismatch detected in SSE implementation" << std::endl;
            std::cerr << "Variable dP_Hat:" << std::endl;
            std::
              cerr << "seed=" << seed << ", __a=" << __a << ", __b=" << __b <<
              std::endl;
            std::cerr << "dP_Hat SSE=  " << dP_Hat[__a][__b] << std::endl;
            std::
              cerr << "dP_Hat Reference=  " << dP_Hat_reference[__a][__b] <<
              std::endl;
            std::cerr << "dP_Hat Rel Difference=  " << std::
              abs ((dP_Hat[__a][__b] -
                    dP_Hat_reference[__a][__b]) /
                   (dP_Hat_reference[__a][__b])) << std::endl;
            std::cerr << "dP_Hat Abs Difference=  " << std::
              abs (dP_Hat[__a][__b] - dP_Hat_reference[__a][__b]) << std::endl;
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
      typedef         T (&refArray2)[9][16];
      typedef         T (&refArray3)[9][16];
      for (int __a = 0; __a < 9; __a++)
        for (int __b = 0; __b < 16; __b++)
          dP_Hat[__a][__b] = dP_Hat_original[__a][__b];
      for (int i = 0; i < 16; i += 8)
      {
        refArray1       dPdFk = reinterpret_cast < refArray1 > (dPdF[0][i]);
        refArray2       dF_Hatk = reinterpret_cast < refArray2 > (dF_Hat[0][i]);
        refArray3       dP_Hatk = reinterpret_cast < refArray3 > (dP_Hat[0][i]);
        Rotated_Stress_Derivative_Contraction < __m256, float[16],
          int[16] > (dPdFk, dF_Hatk, dP_Hatk);
      }
      for (int __a = 0; __a < 9; __a++)
        for (int __b = 0; __b < 16; __b++)
          if (std::
              abs ((dP_Hat[__a][__b] -
                    dP_Hat_reference[__a][__b]) /
                   (dP_Hat_reference[__a][__b])) > 1)
          {
            std::cerr << "Mismatch detected in AVX implementation" << std::endl;
            std::cerr << "Variable dP_Hat:" << std::endl;
            std::
              cerr << "seed=" << seed << ", __a=" << __a << ", __b=" << __b <<
              std::endl;
            std::cerr << "dP_Hat AVX=  " << dP_Hat[__a][__b] << std::endl;
            std::
              cerr << "dP_Hat Reference=  " << dP_Hat_reference[__a][__b] <<
              std::endl;
            std::cerr << "dP_Hat Rel Difference=  " << std::
              abs ((dP_Hat[__a][__b] -
                    dP_Hat_reference[__a][__b]) /
                   (dP_Hat_reference[__a][__b])) << std::endl;
            std::cerr << "dP_Hat Abs Difference=  " << std::
              abs (dP_Hat[__a][__b] - dP_Hat_reference[__a][__b]) << std::endl;
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
      typedef         T (&refArray2)[9][16];
      typedef         T (&refArray3)[9][16];
      for (int __a = 0; __a < 9; __a++)
        for (int __b = 0; __b < 16; __b++)
          dP_Hat[__a][__b] = dP_Hat_original[__a][__b];
      for (int i = 0; i < 16; i += 4)
      {
        refArray1       dPdFk = reinterpret_cast < refArray1 > (dPdF[0][i]);
        refArray2       dF_Hatk = reinterpret_cast < refArray2 > (dF_Hat[0][i]);
        refArray3       dP_Hatk = reinterpret_cast < refArray3 > (dP_Hat[0][i]);
        Rotated_Stress_Derivative_Contraction < float32x4_t, float[16],
          int[16] > (dPdFk, dF_Hatk, dP_Hatk);
      }
      for (int __a = 0; __a < 9; __a++)
        for (int __b = 0; __b < 16; __b++)
          if (std::
              abs ((dP_Hat[__a][__b] -
                    dP_Hat_reference[__a][__b]) /
                   (dP_Hat_reference[__a][__b])) > 1)
          {
            std::cerr << "Mismatch detected in NEON implementation" << std::
              endl;
            std::cerr << "Variable dP_Hat:" << std::endl;
            std::
              cerr << "seed=" << seed << ", __a=" << __a << ", __b=" << __b <<
              std::endl;
            std::cerr << "dP_Hat NEON=  " << dP_Hat[__a][__b] << std::endl;
            std::
              cerr << "dP_Hat Reference=  " << dP_Hat_reference[__a][__b] <<
              std::endl;
            std::cerr << "dP_Hat Rel Difference=  " << std::
              abs ((dP_Hat[__a][__b] -
                    dP_Hat_reference[__a][__b]) /
                   (dP_Hat_reference[__a][__b])) << std::endl;
            std::cerr << "dP_Hat Abs Difference=  " << std::
              abs (dP_Hat[__a][__b] - dP_Hat_reference[__a][__b]) << std::endl;
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
      typedef         T (&refArray2)[9][16];
      typedef         T (&refArray3)[9][16];
      for (int __a = 0; __a < 9; __a++)
        for (int __b = 0; __b < 16; __b++)
          dP_Hat[__a][__b] = dP_Hat_original[__a][__b];
      for (int i = 0; i < 16; i += 16)
      {
        refArray1       dPdFk = reinterpret_cast < refArray1 > (dPdF[0][i]);
        refArray2       dF_Hatk = reinterpret_cast < refArray2 > (dF_Hat[0][i]);
        refArray3       dP_Hatk = reinterpret_cast < refArray3 > (dP_Hat[0][i]);
        Rotated_Stress_Derivative_Contraction < __m512, float[16],
          int[16] > (dPdFk, dF_Hatk, dP_Hatk);
      }
      for (int __a = 0; __a < 9; __a++)
        for (int __b = 0; __b < 16; __b++)
          if (std::
              abs ((dP_Hat[__a][__b] -
                    dP_Hat_reference[__a][__b]) /
                   (dP_Hat_reference[__a][__b])) > 1)
          {
            std::cerr << "Mismatch detected in MIC implementation" << std::endl;
            std::cerr << "Variable dP_Hat:" << std::endl;
            std::
              cerr << "seed=" << seed << ", __a=" << __a << ", __b=" << __b <<
              std::endl;
            std::cerr << "dP_Hat MIC=  " << dP_Hat[__a][__b] << std::endl;
            std::
              cerr << "dP_Hat Reference=  " << dP_Hat_reference[__a][__b] <<
              std::endl;
            std::cerr << "dP_Hat Rel Difference=  " << std::
              abs ((dP_Hat[__a][__b] -
                    dP_Hat_reference[__a][__b]) /
                   (dP_Hat_reference[__a][__b])) << std::endl;
            std::cerr << "dP_Hat Abs Difference=  " << std::
              abs (dP_Hat[__a][__b] - dP_Hat_reference[__a][__b]) << std::endl;
            return 1;
          }

    }
#endif

  }



  std::cout << "SIMD check successful!" << std::endl;

  return 0;

}

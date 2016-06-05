
#include <cstdlib>
#include <iostream>
#include "KernelCommon.h"

struct NEOHOOKEAN_TAG;
struct COROTATED_TAG;

#include "Unweighted_Gradient.h"
#include "Unweighted_Gradient_Reference.h"

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
    std::cout << "Running SIMD Test for Unweighted_Gradient " << std::endl;


//=======================================================
//
//        DEFINE ALL VARIABLES USED BY KERNEL
//
//=======================================================

    T               u[3][8][16] __attribute__ ((aligned (64)));
    T               F[9][16] __attribute__ ((aligned (64)));
    T               F_reference[9][16] __attribute__ ((aligned (64)));
    T               F_original[9][16] __attribute__ ((aligned (64)));
    T               one_over_h[16] __attribute__ ((aligned (64)));


    for (int __a = 0; __a < 3; __a++)
      for (int __b = 0; __b < 8; __b++)
        for (int __c = 0; __c < 16; __c++)
          u[__a][__b][__c] = Get_Random < float >();
    for (int __a = 0; __a < 9; __a++)
      for (int __b = 0; __b < 16; __b++)
      {
        F_original[__a][__b] = Get_Random < float >();
        F[__a][__b] = F_original[__a][__b];
        F_reference[__a][__b] = F_original[__a][__b];
      }
    for (int __a = 0; __a < 16; __a++)
      one_over_h[__a] = Get_Random < float >();

//=======================================================
//
//             COMPUTE REFERENCE RESULTS
//
//=======================================================

    T               __mu[3][8] __attribute__ ((aligned (4)));
    T               __mF[9] __attribute__ ((aligned (4)));
    T               __mF_reference[9] __attribute__ ((aligned (4)));
    T               __mF_original[9] __attribute__ ((aligned (4)));
    T __mone_over_h __attribute__ ((aligned (4)));
    for (int k = 0; k < 16; k++)
    {
      for (int __a = 0; __a < 3; __a++)
        for (int __b = 0; __b < 8; __b++)
          __mu[__a][__b] = u[__a][__b][k];
      for (int __a = 0; __a < 9; __a++)
        __mF_reference[__a] = F_reference[__a][k];

      __mone_over_h = one_over_h[k];
      Unweighted_Gradient < float, float, int >(__mu, __mF_reference,
                                                __mone_over_h);
      for (int __a = 0; __a < 9; __a++)
        F_reference[__a][k] = __mF_reference[__a];
    }

//=======================================================
//
//               COMPUTE SCALAR RESULTS
//
//=======================================================

    {
      typedef         T (&refArray1)[3][8][16];
      typedef         T (&refArray2)[9][16];
      typedef         T (&refArray3)[16];
      for (int __a = 0; __a < 9; __a++)
        for (int __b = 0; __b < 16; __b++)
          F[__a][__b] = F_original[__a][__b];
      for (int i = 0; i < 16; i += 1)
      {
        refArray1       uk = reinterpret_cast < refArray1 > (u[0][0][i]);
        refArray2       Fk = reinterpret_cast < refArray2 > (F[0][i]);
        refArray3       one_over_hk =
          reinterpret_cast < refArray3 > (one_over_h[i]);
        Unweighted_Gradient < float, float[16], int[16] > (uk, Fk, one_over_hk);
      }
      for (int __a = 0; __a < 9; __a++)
        for (int __b = 0; __b < 16; __b++)
          if (std::
              abs ((F[__a][__b] -
                    F_reference[__a][__b]) / (F_reference[__a][__b])) > 1)
          {
            std::cerr << "Mismatch detected in SCALAR implementation" << std::
              endl;
            std::cerr << "Variable F:" << std::endl;
            std::
              cerr << "seed=" << seed << ", __a=" << __a << ", __b=" << __b <<
              std::endl;
            std::cerr << "F SCALAR=  " << F[__a][__b] << std::endl;
            std::cerr << "F Reference=  " << F_reference[__a][__b] << std::endl;
            std::cerr << "F Rel Difference=  " << std::
              abs ((F[__a][__b] -
                    F_reference[__a][__b]) /
                   (F_reference[__a][__b])) << std::endl;
            std::cerr << "F Abs Difference=  " << std::abs (F[__a][__b] -
                                                            F_reference[__a]
                                                            [__b]) << std::endl;
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
      typedef         T (&refArray1)[3][8][16];
      typedef         T (&refArray2)[9][16];
      typedef         T (&refArray3)[16];
      for (int __a = 0; __a < 9; __a++)
        for (int __b = 0; __b < 16; __b++)
          F[__a][__b] = F_original[__a][__b];
      for (int i = 0; i < 16; i += 4)
      {
        refArray1       uk = reinterpret_cast < refArray1 > (u[0][0][i]);
        refArray2       Fk = reinterpret_cast < refArray2 > (F[0][i]);
        refArray3       one_over_hk =
          reinterpret_cast < refArray3 > (one_over_h[i]);
        Unweighted_Gradient < __m128, float[16], int[16] > (uk, Fk,
                                                            one_over_hk);
      }
      for (int __a = 0; __a < 9; __a++)
        for (int __b = 0; __b < 16; __b++)
          if (std::
              abs ((F[__a][__b] -
                    F_reference[__a][__b]) / (F_reference[__a][__b])) > 1)
          {
            std::cerr << "Mismatch detected in SSE implementation" << std::endl;
            std::cerr << "Variable F:" << std::endl;
            std::
              cerr << "seed=" << seed << ", __a=" << __a << ", __b=" << __b <<
              std::endl;
            std::cerr << "F SSE=  " << F[__a][__b] << std::endl;
            std::cerr << "F Reference=  " << F_reference[__a][__b] << std::endl;
            std::cerr << "F Rel Difference=  " << std::
              abs ((F[__a][__b] -
                    F_reference[__a][__b]) /
                   (F_reference[__a][__b])) << std::endl;
            std::cerr << "F Abs Difference=  " << std::abs (F[__a][__b] -
                                                            F_reference[__a]
                                                            [__b]) << std::endl;
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
      typedef         T (&refArray1)[3][8][16];
      typedef         T (&refArray2)[9][16];
      typedef         T (&refArray3)[16];
      for (int __a = 0; __a < 9; __a++)
        for (int __b = 0; __b < 16; __b++)
          F[__a][__b] = F_original[__a][__b];
      for (int i = 0; i < 16; i += 8)
      {
        refArray1       uk = reinterpret_cast < refArray1 > (u[0][0][i]);
        refArray2       Fk = reinterpret_cast < refArray2 > (F[0][i]);
        refArray3       one_over_hk =
          reinterpret_cast < refArray3 > (one_over_h[i]);
        Unweighted_Gradient < __m256, float[16], int[16] > (uk, Fk,
                                                            one_over_hk);
      }
      for (int __a = 0; __a < 9; __a++)
        for (int __b = 0; __b < 16; __b++)
          if (std::
              abs ((F[__a][__b] -
                    F_reference[__a][__b]) / (F_reference[__a][__b])) > 1)
          {
            std::cerr << "Mismatch detected in AVX implementation" << std::endl;
            std::cerr << "Variable F:" << std::endl;
            std::
              cerr << "seed=" << seed << ", __a=" << __a << ", __b=" << __b <<
              std::endl;
            std::cerr << "F AVX=  " << F[__a][__b] << std::endl;
            std::cerr << "F Reference=  " << F_reference[__a][__b] << std::endl;
            std::cerr << "F Rel Difference=  " << std::
              abs ((F[__a][__b] -
                    F_reference[__a][__b]) /
                   (F_reference[__a][__b])) << std::endl;
            std::cerr << "F Abs Difference=  " << std::abs (F[__a][__b] -
                                                            F_reference[__a]
                                                            [__b]) << std::endl;
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
      typedef         T (&refArray1)[3][8][16];
      typedef         T (&refArray2)[9][16];
      typedef         T (&refArray3)[16];
      for (int __a = 0; __a < 9; __a++)
        for (int __b = 0; __b < 16; __b++)
          F[__a][__b] = F_original[__a][__b];
      for (int i = 0; i < 16; i += 4)
      {
        refArray1       uk = reinterpret_cast < refArray1 > (u[0][0][i]);
        refArray2       Fk = reinterpret_cast < refArray2 > (F[0][i]);
        refArray3       one_over_hk =
          reinterpret_cast < refArray3 > (one_over_h[i]);
        Unweighted_Gradient < float32x4_t, float[16], int[16] > (uk, Fk,
                                                                 one_over_hk);
      }
      for (int __a = 0; __a < 9; __a++)
        for (int __b = 0; __b < 16; __b++)
          if (std::
              abs ((F[__a][__b] -
                    F_reference[__a][__b]) / (F_reference[__a][__b])) > 1)
          {
            std::cerr << "Mismatch detected in NEON implementation" << std::
              endl;
            std::cerr << "Variable F:" << std::endl;
            std::
              cerr << "seed=" << seed << ", __a=" << __a << ", __b=" << __b <<
              std::endl;
            std::cerr << "F NEON=  " << F[__a][__b] << std::endl;
            std::cerr << "F Reference=  " << F_reference[__a][__b] << std::endl;
            std::cerr << "F Rel Difference=  " << std::
              abs ((F[__a][__b] -
                    F_reference[__a][__b]) /
                   (F_reference[__a][__b])) << std::endl;
            std::cerr << "F Abs Difference=  " << std::abs (F[__a][__b] -
                                                            F_reference[__a]
                                                            [__b]) << std::endl;
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
      typedef         T (&refArray1)[3][8][16];
      typedef         T (&refArray2)[9][16];
      typedef         T (&refArray3)[16];
      for (int __a = 0; __a < 9; __a++)
        for (int __b = 0; __b < 16; __b++)
          F[__a][__b] = F_original[__a][__b];
      for (int i = 0; i < 16; i += 16)
      {
        refArray1       uk = reinterpret_cast < refArray1 > (u[0][0][i]);
        refArray2       Fk = reinterpret_cast < refArray2 > (F[0][i]);
        refArray3       one_over_hk =
          reinterpret_cast < refArray3 > (one_over_h[i]);
        Unweighted_Gradient < __m512, float[16], int[16] > (uk, Fk,
                                                            one_over_hk);
      }
      for (int __a = 0; __a < 9; __a++)
        for (int __b = 0; __b < 16; __b++)
          if (std::
              abs ((F[__a][__b] -
                    F_reference[__a][__b]) / (F_reference[__a][__b])) > 1)
          {
            std::cerr << "Mismatch detected in MIC implementation" << std::endl;
            std::cerr << "Variable F:" << std::endl;
            std::
              cerr << "seed=" << seed << ", __a=" << __a << ", __b=" << __b <<
              std::endl;
            std::cerr << "F MIC=  " << F[__a][__b] << std::endl;
            std::cerr << "F Reference=  " << F_reference[__a][__b] << std::endl;
            std::cerr << "F Rel Difference=  " << std::
              abs ((F[__a][__b] -
                    F_reference[__a][__b]) /
                   (F_reference[__a][__b])) << std::endl;
            std::cerr << "F Abs Difference=  " << std::abs (F[__a][__b] -
                                                            F_reference[__a]
                                                            [__b]) << std::endl;
            return 1;
          }

    }
#endif

  }



  std::cout << "SIMD check successful!" << std::endl;

  return 0;

}

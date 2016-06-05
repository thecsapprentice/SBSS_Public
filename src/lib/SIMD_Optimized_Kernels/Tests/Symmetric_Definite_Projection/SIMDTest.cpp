
#include <cstdlib>
#include <iostream>
#include "KernelCommon.h"

struct NEOHOOKEAN_TAG;
struct COROTATED_TAG;

#include "Symmetric_Definite_Projection.h"
#include "Symmetric_Definite_Projection_Reference.h"

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
    std::cout << "Running SIMD Test for Symmetric_Definite_Projection " << std::
      endl;


//=======================================================
//
//        DEFINE ALL VARIABLES USED BY KERNEL
//
//=======================================================

    T               A[6][16] __attribute__ ((aligned (64)));
    T               Apd[6][16] __attribute__ ((aligned (64)));
    T               Apd_reference[6][16] __attribute__ ((aligned (64)));
    T               Apd_original[6][16] __attribute__ ((aligned (64)));


    for (int __a = 0; __a < 6; __a++)
      for (int __b = 0; __b < 16; __b++)
        A[__a][__b] = Get_Random < float >();
    for (int __a = 0; __a < 6; __a++)
      for (int __b = 0; __b < 16; __b++)
      {
        Apd_original[__a][__b] = Get_Random < float >();
        Apd[__a][__b] = Apd_original[__a][__b];
        Apd_reference[__a][__b] = Apd_original[__a][__b];
      }


//=======================================================
//
//             COMPUTE REFERENCE RESULTS
//
//=======================================================

    T               __mA[6] __attribute__ ((aligned (4)));
    T               __mApd[6] __attribute__ ((aligned (4)));
    T               __mApd_reference[6] __attribute__ ((aligned (4)));
    T               __mApd_original[6] __attribute__ ((aligned (4)));
    for (int k = 0; k < 16; k++)
    {
      for (int __a = 0; __a < 6; __a++)
        __mA[__a] = A[__a][k];
      for (int __a = 0; __a < 6; __a++)
        __mApd_reference[__a] = Apd_reference[__a][k];
      Symmetric_Definite_Projection < float, float, int >(__mA,
                                                          __mApd_reference);
      for (int __a = 0; __a < 6; __a++)
        Apd_reference[__a][k] = __mApd_reference[__a];
    }

//=======================================================
//
//               COMPUTE SCALAR RESULTS
//
//=======================================================

    {
      typedef         T (&refArray1)[6][16];
      typedef         T (&refArray2)[6][16];
      for (int __a = 0; __a < 6; __a++)
        for (int __b = 0; __b < 16; __b++)
          Apd[__a][__b] = Apd_original[__a][__b];
      for (int i = 0; i < 16; i += 1)
      {
        refArray1       Ak = reinterpret_cast < refArray1 > (A[0][i]);
        refArray2       Apdk = reinterpret_cast < refArray2 > (Apd[0][i]);
        Symmetric_Definite_Projection < float, float[16], int[16] > (Ak, Apdk);
      }
      for (int __a = 0; __a < 6; __a++)
        for (int __b = 0; __b < 16; __b++)
          if (std::
              abs ((Apd[__a][__b] -
                    Apd_reference[__a][__b]) / (Apd_reference[__a][__b])) > 1)
          {
            std::cerr << "Mismatch detected in SCALAR implementation" << std::
              endl;
            std::cerr << "Variable Apd:" << std::endl;
            std::
              cerr << "seed=" << seed << ", __a=" << __a << ", __b=" << __b <<
              std::endl;
            std::cerr << "Apd SCALAR=  " << Apd[__a][__b] << std::endl;
            std::cerr << "Apd Reference=  " << Apd_reference[__a][__b] << std::
              endl;
            std::cerr << "Apd Rel Difference=  " << std::
              abs ((Apd[__a][__b] -
                    Apd_reference[__a][__b]) /
                   (Apd_reference[__a][__b])) << std::endl;
            std::cerr << "Apd Abs Difference=  " << std::abs (Apd[__a][__b] -
                                                              Apd_reference[__a]
                                                              [__b]) << std::
              endl;
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
      typedef         T (&refArray1)[6][16];
      typedef         T (&refArray2)[6][16];
      for (int __a = 0; __a < 6; __a++)
        for (int __b = 0; __b < 16; __b++)
          Apd[__a][__b] = Apd_original[__a][__b];
      for (int i = 0; i < 16; i += 4)
      {
        refArray1       Ak = reinterpret_cast < refArray1 > (A[0][i]);
        refArray2       Apdk = reinterpret_cast < refArray2 > (Apd[0][i]);
        Symmetric_Definite_Projection < __m128, float[16], int[16] > (Ak, Apdk);
      }
      for (int __a = 0; __a < 6; __a++)
        for (int __b = 0; __b < 16; __b++)
          if (std::
              abs ((Apd[__a][__b] -
                    Apd_reference[__a][__b]) / (Apd_reference[__a][__b])) > 1)
          {
            std::cerr << "Mismatch detected in SSE implementation" << std::endl;
            std::cerr << "Variable Apd:" << std::endl;
            std::
              cerr << "seed=" << seed << ", __a=" << __a << ", __b=" << __b <<
              std::endl;
            std::cerr << "Apd SSE=  " << Apd[__a][__b] << std::endl;
            std::cerr << "Apd Reference=  " << Apd_reference[__a][__b] << std::
              endl;
            std::cerr << "Apd Rel Difference=  " << std::
              abs ((Apd[__a][__b] -
                    Apd_reference[__a][__b]) /
                   (Apd_reference[__a][__b])) << std::endl;
            std::cerr << "Apd Abs Difference=  " << std::abs (Apd[__a][__b] -
                                                              Apd_reference[__a]
                                                              [__b]) << std::
              endl;
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
      typedef         T (&refArray1)[6][16];
      typedef         T (&refArray2)[6][16];
      for (int __a = 0; __a < 6; __a++)
        for (int __b = 0; __b < 16; __b++)
          Apd[__a][__b] = Apd_original[__a][__b];
      for (int i = 0; i < 16; i += 8)
      {
        refArray1       Ak = reinterpret_cast < refArray1 > (A[0][i]);
        refArray2       Apdk = reinterpret_cast < refArray2 > (Apd[0][i]);
        Symmetric_Definite_Projection < __m256, float[16], int[16] > (Ak, Apdk);
      }
      for (int __a = 0; __a < 6; __a++)
        for (int __b = 0; __b < 16; __b++)
          if (std::
              abs ((Apd[__a][__b] -
                    Apd_reference[__a][__b]) / (Apd_reference[__a][__b])) > 1)
          {
            std::cerr << "Mismatch detected in AVX implementation" << std::endl;
            std::cerr << "Variable Apd:" << std::endl;
            std::
              cerr << "seed=" << seed << ", __a=" << __a << ", __b=" << __b <<
              std::endl;
            std::cerr << "Apd AVX=  " << Apd[__a][__b] << std::endl;
            std::cerr << "Apd Reference=  " << Apd_reference[__a][__b] << std::
              endl;
            std::cerr << "Apd Rel Difference=  " << std::
              abs ((Apd[__a][__b] -
                    Apd_reference[__a][__b]) /
                   (Apd_reference[__a][__b])) << std::endl;
            std::cerr << "Apd Abs Difference=  " << std::abs (Apd[__a][__b] -
                                                              Apd_reference[__a]
                                                              [__b]) << std::
              endl;
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
      typedef         T (&refArray1)[6][16];
      typedef         T (&refArray2)[6][16];
      for (int __a = 0; __a < 6; __a++)
        for (int __b = 0; __b < 16; __b++)
          Apd[__a][__b] = Apd_original[__a][__b];
      for (int i = 0; i < 16; i += 4)
      {
        refArray1       Ak = reinterpret_cast < refArray1 > (A[0][i]);
        refArray2       Apdk = reinterpret_cast < refArray2 > (Apd[0][i]);
        Symmetric_Definite_Projection < float32x4_t, float[16], int[16] > (Ak,
                                                                           Apdk);
      }
      for (int __a = 0; __a < 6; __a++)
        for (int __b = 0; __b < 16; __b++)
          if (std::
              abs ((Apd[__a][__b] -
                    Apd_reference[__a][__b]) / (Apd_reference[__a][__b])) > 1)
          {
            std::cerr << "Mismatch detected in NEON implementation" << std::
              endl;
            std::cerr << "Variable Apd:" << std::endl;
            std::
              cerr << "seed=" << seed << ", __a=" << __a << ", __b=" << __b <<
              std::endl;
            std::cerr << "Apd NEON=  " << Apd[__a][__b] << std::endl;
            std::cerr << "Apd Reference=  " << Apd_reference[__a][__b] << std::
              endl;
            std::cerr << "Apd Rel Difference=  " << std::
              abs ((Apd[__a][__b] -
                    Apd_reference[__a][__b]) /
                   (Apd_reference[__a][__b])) << std::endl;
            std::cerr << "Apd Abs Difference=  " << std::abs (Apd[__a][__b] -
                                                              Apd_reference[__a]
                                                              [__b]) << std::
              endl;
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
      typedef         T (&refArray1)[6][16];
      typedef         T (&refArray2)[6][16];
      for (int __a = 0; __a < 6; __a++)
        for (int __b = 0; __b < 16; __b++)
          Apd[__a][__b] = Apd_original[__a][__b];
      for (int i = 0; i < 16; i += 16)
      {
        refArray1       Ak = reinterpret_cast < refArray1 > (A[0][i]);
        refArray2       Apdk = reinterpret_cast < refArray2 > (Apd[0][i]);
        Symmetric_Definite_Projection < __m512, float[16], int[16] > (Ak, Apdk);
      }
      for (int __a = 0; __a < 6; __a++)
        for (int __b = 0; __b < 16; __b++)
          if (std::
              abs ((Apd[__a][__b] -
                    Apd_reference[__a][__b]) / (Apd_reference[__a][__b])) > 1)
          {
            std::cerr << "Mismatch detected in MIC implementation" << std::endl;
            std::cerr << "Variable Apd:" << std::endl;
            std::
              cerr << "seed=" << seed << ", __a=" << __a << ", __b=" << __b <<
              std::endl;
            std::cerr << "Apd MIC=  " << Apd[__a][__b] << std::endl;
            std::cerr << "Apd Reference=  " << Apd_reference[__a][__b] << std::
              endl;
            std::cerr << "Apd Rel Difference=  " << std::
              abs ((Apd[__a][__b] -
                    Apd_reference[__a][__b]) /
                   (Apd_reference[__a][__b])) << std::endl;
            std::cerr << "Apd Abs Difference=  " << std::abs (Apd[__a][__b] -
                                                              Apd_reference[__a]
                                                              [__b]) << std::
              endl;
            return 1;
          }

    }
#endif

  }



  std::cout << "SIMD check successful!" << std::endl;

  return 0;

}

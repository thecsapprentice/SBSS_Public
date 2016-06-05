
#include <cstdlib>
#include <iostream>
#include "KernelCommon.h"

struct NEOHOOKEAN_TAG;
struct COROTATED_TAG;

#include "Matrix_Times_Transpose.h"
#include "Matrix_Times_Transpose_Reference.h"

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
    std::cout << "Running SIMD Test for Matrix_Times_Transpose " << std::endl;


//=======================================================
//
//        DEFINE ALL VARIABLES USED BY KERNEL
//
//=======================================================

    T               A[9][16] __attribute__ ((aligned (64)));
    T               B[9][16] __attribute__ ((aligned (64)));
    T               C[9][16] __attribute__ ((aligned (64)));
    T               C_reference[9][16] __attribute__ ((aligned (64)));
    T               C_original[9][16] __attribute__ ((aligned (64)));


    for (int __a = 0; __a < 9; __a++)
      for (int __b = 0; __b < 16; __b++)
        A[__a][__b] = Get_Random < float >();
    for (int __a = 0; __a < 9; __a++)
      for (int __b = 0; __b < 16; __b++)
        B[__a][__b] = Get_Random < float >();
    for (int __a = 0; __a < 9; __a++)
      for (int __b = 0; __b < 16; __b++)
      {
        C_original[__a][__b] = Get_Random < float >();
        C[__a][__b] = C_original[__a][__b];
        C_reference[__a][__b] = C_original[__a][__b];
      }


//=======================================================
//
//             COMPUTE REFERENCE RESULTS
//
//=======================================================

    T               __mA[9] __attribute__ ((aligned (4)));
    T               __mB[9] __attribute__ ((aligned (4)));
    T               __mC[9] __attribute__ ((aligned (4)));
    T               __mC_reference[9] __attribute__ ((aligned (4)));
    T               __mC_original[9] __attribute__ ((aligned (4)));
    for (int k = 0; k < 16; k++)
    {
      for (int __a = 0; __a < 9; __a++)
        __mA[__a] = A[__a][k];
      for (int __a = 0; __a < 9; __a++)
        __mB[__a] = B[__a][k];
      for (int __a = 0; __a < 9; __a++)
        __mC_reference[__a] = C_reference[__a][k];
      Matrix_Times_Transpose < float, float, int >(__mA, __mB, __mC_reference);
      for (int __a = 0; __a < 9; __a++)
        C_reference[__a][k] = __mC_reference[__a];
    }

//=======================================================
//
//               COMPUTE SCALAR RESULTS
//
//=======================================================

    {
      typedef         T (&refArray1)[9][16];
      typedef         T (&refArray2)[9][16];
      typedef         T (&refArray3)[9][16];
      for (int __a = 0; __a < 9; __a++)
        for (int __b = 0; __b < 16; __b++)
          C[__a][__b] = C_original[__a][__b];
      for (int i = 0; i < 16; i += 1)
      {
        refArray1       Ak = reinterpret_cast < refArray1 > (A[0][i]);
        refArray2       Bk = reinterpret_cast < refArray2 > (B[0][i]);
        refArray3       Ck = reinterpret_cast < refArray3 > (C[0][i]);
        Matrix_Times_Transpose < float, float[16], int[16] > (Ak, Bk, Ck);
      }
      for (int __a = 0; __a < 9; __a++)
        for (int __b = 0; __b < 16; __b++)
          if (std::
              abs ((C[__a][__b] -
                    C_reference[__a][__b]) / (C_reference[__a][__b])) > 1)
          {
            std::cerr << "Mismatch detected in SCALAR implementation" << std::
              endl;
            std::cerr << "Variable C:" << std::endl;
            std::
              cerr << "seed=" << seed << ", __a=" << __a << ", __b=" << __b <<
              std::endl;
            std::cerr << "C SCALAR=  " << C[__a][__b] << std::endl;
            std::cerr << "C Reference=  " << C_reference[__a][__b] << std::endl;
            std::cerr << "C Rel Difference=  " << std::
              abs ((C[__a][__b] -
                    C_reference[__a][__b]) /
                   (C_reference[__a][__b])) << std::endl;
            std::cerr << "C Abs Difference=  " << std::abs (C[__a][__b] -
                                                            C_reference[__a]
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
      typedef         T (&refArray1)[9][16];
      typedef         T (&refArray2)[9][16];
      typedef         T (&refArray3)[9][16];
      for (int __a = 0; __a < 9; __a++)
        for (int __b = 0; __b < 16; __b++)
          C[__a][__b] = C_original[__a][__b];
      for (int i = 0; i < 16; i += 4)
      {
        refArray1       Ak = reinterpret_cast < refArray1 > (A[0][i]);
        refArray2       Bk = reinterpret_cast < refArray2 > (B[0][i]);
        refArray3       Ck = reinterpret_cast < refArray3 > (C[0][i]);
        Matrix_Times_Transpose < __m128, float[16], int[16] > (Ak, Bk, Ck);
      }
      for (int __a = 0; __a < 9; __a++)
        for (int __b = 0; __b < 16; __b++)
          if (std::
              abs ((C[__a][__b] -
                    C_reference[__a][__b]) / (C_reference[__a][__b])) > 1)
          {
            std::cerr << "Mismatch detected in SSE implementation" << std::endl;
            std::cerr << "Variable C:" << std::endl;
            std::
              cerr << "seed=" << seed << ", __a=" << __a << ", __b=" << __b <<
              std::endl;
            std::cerr << "C SSE=  " << C[__a][__b] << std::endl;
            std::cerr << "C Reference=  " << C_reference[__a][__b] << std::endl;
            std::cerr << "C Rel Difference=  " << std::
              abs ((C[__a][__b] -
                    C_reference[__a][__b]) /
                   (C_reference[__a][__b])) << std::endl;
            std::cerr << "C Abs Difference=  " << std::abs (C[__a][__b] -
                                                            C_reference[__a]
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
      typedef         T (&refArray1)[9][16];
      typedef         T (&refArray2)[9][16];
      typedef         T (&refArray3)[9][16];
      for (int __a = 0; __a < 9; __a++)
        for (int __b = 0; __b < 16; __b++)
          C[__a][__b] = C_original[__a][__b];
      for (int i = 0; i < 16; i += 8)
      {
        refArray1       Ak = reinterpret_cast < refArray1 > (A[0][i]);
        refArray2       Bk = reinterpret_cast < refArray2 > (B[0][i]);
        refArray3       Ck = reinterpret_cast < refArray3 > (C[0][i]);
        Matrix_Times_Transpose < __m256, float[16], int[16] > (Ak, Bk, Ck);
      }
      for (int __a = 0; __a < 9; __a++)
        for (int __b = 0; __b < 16; __b++)
          if (std::
              abs ((C[__a][__b] -
                    C_reference[__a][__b]) / (C_reference[__a][__b])) > 1)
          {
            std::cerr << "Mismatch detected in AVX implementation" << std::endl;
            std::cerr << "Variable C:" << std::endl;
            std::
              cerr << "seed=" << seed << ", __a=" << __a << ", __b=" << __b <<
              std::endl;
            std::cerr << "C AVX=  " << C[__a][__b] << std::endl;
            std::cerr << "C Reference=  " << C_reference[__a][__b] << std::endl;
            std::cerr << "C Rel Difference=  " << std::
              abs ((C[__a][__b] -
                    C_reference[__a][__b]) /
                   (C_reference[__a][__b])) << std::endl;
            std::cerr << "C Abs Difference=  " << std::abs (C[__a][__b] -
                                                            C_reference[__a]
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
      typedef         T (&refArray1)[9][16];
      typedef         T (&refArray2)[9][16];
      typedef         T (&refArray3)[9][16];
      for (int __a = 0; __a < 9; __a++)
        for (int __b = 0; __b < 16; __b++)
          C[__a][__b] = C_original[__a][__b];
      for (int i = 0; i < 16; i += 4)
      {
        refArray1       Ak = reinterpret_cast < refArray1 > (A[0][i]);
        refArray2       Bk = reinterpret_cast < refArray2 > (B[0][i]);
        refArray3       Ck = reinterpret_cast < refArray3 > (C[0][i]);
        Matrix_Times_Transpose < float32x4_t, float[16], int[16] > (Ak, Bk, Ck);
      }
      for (int __a = 0; __a < 9; __a++)
        for (int __b = 0; __b < 16; __b++)
          if (std::
              abs ((C[__a][__b] -
                    C_reference[__a][__b]) / (C_reference[__a][__b])) > 1)
          {
            std::cerr << "Mismatch detected in NEON implementation" << std::
              endl;
            std::cerr << "Variable C:" << std::endl;
            std::
              cerr << "seed=" << seed << ", __a=" << __a << ", __b=" << __b <<
              std::endl;
            std::cerr << "C NEON=  " << C[__a][__b] << std::endl;
            std::cerr << "C Reference=  " << C_reference[__a][__b] << std::endl;
            std::cerr << "C Rel Difference=  " << std::
              abs ((C[__a][__b] -
                    C_reference[__a][__b]) /
                   (C_reference[__a][__b])) << std::endl;
            std::cerr << "C Abs Difference=  " << std::abs (C[__a][__b] -
                                                            C_reference[__a]
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
      typedef         T (&refArray1)[9][16];
      typedef         T (&refArray2)[9][16];
      typedef         T (&refArray3)[9][16];
      for (int __a = 0; __a < 9; __a++)
        for (int __b = 0; __b < 16; __b++)
          C[__a][__b] = C_original[__a][__b];
      for (int i = 0; i < 16; i += 16)
      {
        refArray1       Ak = reinterpret_cast < refArray1 > (A[0][i]);
        refArray2       Bk = reinterpret_cast < refArray2 > (B[0][i]);
        refArray3       Ck = reinterpret_cast < refArray3 > (C[0][i]);
        Matrix_Times_Transpose < __m512, float[16], int[16] > (Ak, Bk, Ck);
      }
      for (int __a = 0; __a < 9; __a++)
        for (int __b = 0; __b < 16; __b++)
          if (std::
              abs ((C[__a][__b] -
                    C_reference[__a][__b]) / (C_reference[__a][__b])) > 1)
          {
            std::cerr << "Mismatch detected in MIC implementation" << std::endl;
            std::cerr << "Variable C:" << std::endl;
            std::
              cerr << "seed=" << seed << ", __a=" << __a << ", __b=" << __b <<
              std::endl;
            std::cerr << "C MIC=  " << C[__a][__b] << std::endl;
            std::cerr << "C Reference=  " << C_reference[__a][__b] << std::endl;
            std::cerr << "C Rel Difference=  " << std::
              abs ((C[__a][__b] -
                    C_reference[__a][__b]) /
                   (C_reference[__a][__b])) << std::endl;
            std::cerr << "C Abs Difference=  " << std::abs (C[__a][__b] -
                                                            C_reference[__a]
                                                            [__b]) << std::endl;
            return 1;
          }

    }
#endif

  }



  std::cout << "SIMD check successful!" << std::endl;

  return 0;

}

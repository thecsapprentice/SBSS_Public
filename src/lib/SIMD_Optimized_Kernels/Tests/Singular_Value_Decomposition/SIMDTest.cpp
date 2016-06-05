
#include <cstdlib>
#include <iostream>
#include "KernelCommon.h"

struct NEOHOOKEAN_TAG;
struct COROTATED_TAG;

#include "Singular_Value_Decomposition.h"
#include "Singular_Value_Decomposition_Reference.h"

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
    std::cout << "Running SIMD Test for Singular_Value_Decomposition " << std::
      endl;


//=======================================================
//
//        DEFINE ALL VARIABLES USED BY KERNEL
//
//=======================================================

    T               A[9][16] __attribute__ ((aligned (64)));
    T               U[9][16] __attribute__ ((aligned (64)));
    T               U_reference[9][16] __attribute__ ((aligned (64)));
    T               U_original[9][16] __attribute__ ((aligned (64)));
    T               S[3][16] __attribute__ ((aligned (64)));
    T               S_reference[3][16] __attribute__ ((aligned (64)));
    T               S_original[3][16] __attribute__ ((aligned (64)));
    T               V[9][16] __attribute__ ((aligned (64)));
    T               V_reference[9][16] __attribute__ ((aligned (64)));
    T               V_original[9][16] __attribute__ ((aligned (64)));


    for (int __a = 0; __a < 9; __a++)
      for (int __b = 0; __b < 16; __b++)
        A[__a][__b] = Get_Random < float >();
    for (int __a = 0; __a < 9; __a++)
      for (int __b = 0; __b < 16; __b++)
      {
        U_original[__a][__b] = Get_Random < float >();
        U[__a][__b] = U_original[__a][__b];
        U_reference[__a][__b] = U_original[__a][__b];
      }
    for (int __a = 0; __a < 3; __a++)
      for (int __b = 0; __b < 16; __b++)
      {
        S_original[__a][__b] = Get_Random < float >();
        S[__a][__b] = S_original[__a][__b];
        S_reference[__a][__b] = S_original[__a][__b];
      }
    for (int __a = 0; __a < 9; __a++)
      for (int __b = 0; __b < 16; __b++)
      {
        V_original[__a][__b] = Get_Random < float >();
        V[__a][__b] = V_original[__a][__b];
        V_reference[__a][__b] = V_original[__a][__b];
      }


//=======================================================
//
//             COMPUTE REFERENCE RESULTS
//
//=======================================================

    T               __mA[9] __attribute__ ((aligned (4)));
    T               __mU[9] __attribute__ ((aligned (4)));
    T               __mU_reference[9] __attribute__ ((aligned (4)));
    T               __mU_original[9] __attribute__ ((aligned (4)));
    T               __mS[3] __attribute__ ((aligned (4)));
    T               __mS_reference[3] __attribute__ ((aligned (4)));
    T               __mS_original[3] __attribute__ ((aligned (4)));
    T               __mV[9] __attribute__ ((aligned (4)));
    T               __mV_reference[9] __attribute__ ((aligned (4)));
    T               __mV_original[9] __attribute__ ((aligned (4)));
    for (int k = 0; k < 16; k++)
    {
      for (int __a = 0; __a < 9; __a++)
        __mA[__a] = A[__a][k];
      for (int __a = 0; __a < 9; __a++)
        __mU_reference[__a] = U_reference[__a][k];
      for (int __a = 0; __a < 3; __a++)
        __mS_reference[__a] = S_reference[__a][k];
      for (int __a = 0; __a < 9; __a++)
        __mV_reference[__a] = V_reference[__a][k];
      Singular_Value_Decomposition < float, float, int >(__mA, __mU_reference,
                                                         __mS_reference,
                                                         __mV_reference);
      for (int __a = 0; __a < 9; __a++)
        U_reference[__a][k] = __mU_reference[__a];
      for (int __a = 0; __a < 3; __a++)
        S_reference[__a][k] = __mS_reference[__a];
      for (int __a = 0; __a < 9; __a++)
        V_reference[__a][k] = __mV_reference[__a];
    }

//=======================================================
//
//               COMPUTE SCALAR RESULTS
//
//=======================================================

    {
      typedef         T (&refArray1)[9][16];
      typedef         T (&refArray2)[9][16];
      typedef         T (&refArray3)[3][16];
      typedef         T (&refArray4)[9][16];
      for (int __a = 0; __a < 9; __a++)
        for (int __b = 0; __b < 16; __b++)
          U[__a][__b] = U_original[__a][__b];
      for (int __a = 0; __a < 3; __a++)
        for (int __b = 0; __b < 16; __b++)
          S[__a][__b] = S_original[__a][__b];
      for (int __a = 0; __a < 9; __a++)
        for (int __b = 0; __b < 16; __b++)
          V[__a][__b] = V_original[__a][__b];
      for (int i = 0; i < 16; i += 1)
      {
        refArray1       Ak = reinterpret_cast < refArray1 > (A[0][i]);
        refArray2       Uk = reinterpret_cast < refArray2 > (U[0][i]);
        refArray3       Sk = reinterpret_cast < refArray3 > (S[0][i]);
        refArray4       Vk = reinterpret_cast < refArray4 > (V[0][i]);
        Singular_Value_Decomposition < float, float[16], int[16] > (Ak, Uk, Sk,
                                                                    Vk);
      }
      for (int __a = 0; __a < 9; __a++)
        for (int __b = 0; __b < 16; __b++)
          if (std::
              abs ((U[__a][__b] -
                    U_reference[__a][__b]) / (U_reference[__a][__b])) > 1)
          {
            std::cerr << "Mismatch detected in SCALAR implementation" << std::
              endl;
            std::cerr << "Variable U:" << std::endl;
            std::
              cerr << "seed=" << seed << ", __a=" << __a << ", __b=" << __b <<
              std::endl;
            std::cerr << "U SCALAR=  " << U[__a][__b] << std::endl;
            std::cerr << "U Reference=  " << U_reference[__a][__b] << std::endl;
            std::cerr << "U Rel Difference=  " << std::
              abs ((U[__a][__b] -
                    U_reference[__a][__b]) /
                   (U_reference[__a][__b])) << std::endl;
            std::cerr << "U Abs Difference=  " << std::abs (U[__a][__b] -
                                                            U_reference[__a]
                                                            [__b]) << std::endl;
            return 1;
          }
      for (int __a = 0; __a < 3; __a++)
        for (int __b = 0; __b < 16; __b++)
          if (std::
              abs ((S[__a][__b] -
                    S_reference[__a][__b]) / (S_reference[__a][__b])) > 1)
          {
            std::cerr << "Mismatch detected in SCALAR implementation" << std::
              endl;
            std::cerr << "Variable S:" << std::endl;
            std::
              cerr << "seed=" << seed << ", __a=" << __a << ", __b=" << __b <<
              std::endl;
            std::cerr << "S SCALAR=  " << S[__a][__b] << std::endl;
            std::cerr << "S Reference=  " << S_reference[__a][__b] << std::endl;
            std::cerr << "S Rel Difference=  " << std::
              abs ((S[__a][__b] -
                    S_reference[__a][__b]) /
                   (S_reference[__a][__b])) << std::endl;
            std::cerr << "S Abs Difference=  " << std::abs (S[__a][__b] -
                                                            S_reference[__a]
                                                            [__b]) << std::endl;
            return 1;
          }
      for (int __a = 0; __a < 9; __a++)
        for (int __b = 0; __b < 16; __b++)
          if (std::
              abs ((V[__a][__b] -
                    V_reference[__a][__b]) / (V_reference[__a][__b])) > 1)
          {
            std::cerr << "Mismatch detected in SCALAR implementation" << std::
              endl;
            std::cerr << "Variable V:" << std::endl;
            std::
              cerr << "seed=" << seed << ", __a=" << __a << ", __b=" << __b <<
              std::endl;
            std::cerr << "V SCALAR=  " << V[__a][__b] << std::endl;
            std::cerr << "V Reference=  " << V_reference[__a][__b] << std::endl;
            std::cerr << "V Rel Difference=  " << std::
              abs ((V[__a][__b] -
                    V_reference[__a][__b]) /
                   (V_reference[__a][__b])) << std::endl;
            std::cerr << "V Abs Difference=  " << std::abs (V[__a][__b] -
                                                            V_reference[__a]
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
      typedef         T (&refArray3)[3][16];
      typedef         T (&refArray4)[9][16];
      for (int __a = 0; __a < 9; __a++)
        for (int __b = 0; __b < 16; __b++)
          U[__a][__b] = U_original[__a][__b];
      for (int __a = 0; __a < 3; __a++)
        for (int __b = 0; __b < 16; __b++)
          S[__a][__b] = S_original[__a][__b];
      for (int __a = 0; __a < 9; __a++)
        for (int __b = 0; __b < 16; __b++)
          V[__a][__b] = V_original[__a][__b];
      for (int i = 0; i < 16; i += 4)
      {
        refArray1       Ak = reinterpret_cast < refArray1 > (A[0][i]);
        refArray2       Uk = reinterpret_cast < refArray2 > (U[0][i]);
        refArray3       Sk = reinterpret_cast < refArray3 > (S[0][i]);
        refArray4       Vk = reinterpret_cast < refArray4 > (V[0][i]);
        Singular_Value_Decomposition < __m128, float[16], int[16] > (Ak, Uk, Sk,
                                                                     Vk);
      }
      for (int __a = 0; __a < 9; __a++)
        for (int __b = 0; __b < 16; __b++)
          if (std::
              abs ((U[__a][__b] -
                    U_reference[__a][__b]) / (U_reference[__a][__b])) > 1)
          {
            std::cerr << "Mismatch detected in SSE implementation" << std::endl;
            std::cerr << "Variable U:" << std::endl;
            std::
              cerr << "seed=" << seed << ", __a=" << __a << ", __b=" << __b <<
              std::endl;
            std::cerr << "U SSE=  " << U[__a][__b] << std::endl;
            std::cerr << "U Reference=  " << U_reference[__a][__b] << std::endl;
            std::cerr << "U Rel Difference=  " << std::
              abs ((U[__a][__b] -
                    U_reference[__a][__b]) /
                   (U_reference[__a][__b])) << std::endl;
            std::cerr << "U Abs Difference=  " << std::abs (U[__a][__b] -
                                                            U_reference[__a]
                                                            [__b]) << std::endl;
            return 1;
          }
      for (int __a = 0; __a < 3; __a++)
        for (int __b = 0; __b < 16; __b++)
          if (std::
              abs ((S[__a][__b] -
                    S_reference[__a][__b]) / (S_reference[__a][__b])) > 1)
          {
            std::cerr << "Mismatch detected in SSE implementation" << std::endl;
            std::cerr << "Variable S:" << std::endl;
            std::
              cerr << "seed=" << seed << ", __a=" << __a << ", __b=" << __b <<
              std::endl;
            std::cerr << "S SSE=  " << S[__a][__b] << std::endl;
            std::cerr << "S Reference=  " << S_reference[__a][__b] << std::endl;
            std::cerr << "S Rel Difference=  " << std::
              abs ((S[__a][__b] -
                    S_reference[__a][__b]) /
                   (S_reference[__a][__b])) << std::endl;
            std::cerr << "S Abs Difference=  " << std::abs (S[__a][__b] -
                                                            S_reference[__a]
                                                            [__b]) << std::endl;
            return 1;
          }
      for (int __a = 0; __a < 9; __a++)
        for (int __b = 0; __b < 16; __b++)
          if (std::
              abs ((V[__a][__b] -
                    V_reference[__a][__b]) / (V_reference[__a][__b])) > 1)
          {
            std::cerr << "Mismatch detected in SSE implementation" << std::endl;
            std::cerr << "Variable V:" << std::endl;
            std::
              cerr << "seed=" << seed << ", __a=" << __a << ", __b=" << __b <<
              std::endl;
            std::cerr << "V SSE=  " << V[__a][__b] << std::endl;
            std::cerr << "V Reference=  " << V_reference[__a][__b] << std::endl;
            std::cerr << "V Rel Difference=  " << std::
              abs ((V[__a][__b] -
                    V_reference[__a][__b]) /
                   (V_reference[__a][__b])) << std::endl;
            std::cerr << "V Abs Difference=  " << std::abs (V[__a][__b] -
                                                            V_reference[__a]
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
      typedef         T (&refArray3)[3][16];
      typedef         T (&refArray4)[9][16];
      for (int __a = 0; __a < 9; __a++)
        for (int __b = 0; __b < 16; __b++)
          U[__a][__b] = U_original[__a][__b];
      for (int __a = 0; __a < 3; __a++)
        for (int __b = 0; __b < 16; __b++)
          S[__a][__b] = S_original[__a][__b];
      for (int __a = 0; __a < 9; __a++)
        for (int __b = 0; __b < 16; __b++)
          V[__a][__b] = V_original[__a][__b];
      for (int i = 0; i < 16; i += 8)
      {
        refArray1       Ak = reinterpret_cast < refArray1 > (A[0][i]);
        refArray2       Uk = reinterpret_cast < refArray2 > (U[0][i]);
        refArray3       Sk = reinterpret_cast < refArray3 > (S[0][i]);
        refArray4       Vk = reinterpret_cast < refArray4 > (V[0][i]);
        Singular_Value_Decomposition < __m256, float[16], int[16] > (Ak, Uk, Sk,
                                                                     Vk);
      }
      for (int __a = 0; __a < 9; __a++)
        for (int __b = 0; __b < 16; __b++)
          if (std::
              abs ((U[__a][__b] -
                    U_reference[__a][__b]) / (U_reference[__a][__b])) > 1)
          {
            std::cerr << "Mismatch detected in AVX implementation" << std::endl;
            std::cerr << "Variable U:" << std::endl;
            std::
              cerr << "seed=" << seed << ", __a=" << __a << ", __b=" << __b <<
              std::endl;
            std::cerr << "U AVX=  " << U[__a][__b] << std::endl;
            std::cerr << "U Reference=  " << U_reference[__a][__b] << std::endl;
            std::cerr << "U Rel Difference=  " << std::
              abs ((U[__a][__b] -
                    U_reference[__a][__b]) /
                   (U_reference[__a][__b])) << std::endl;
            std::cerr << "U Abs Difference=  " << std::abs (U[__a][__b] -
                                                            U_reference[__a]
                                                            [__b]) << std::endl;
            return 1;
          }
      for (int __a = 0; __a < 3; __a++)
        for (int __b = 0; __b < 16; __b++)
          if (std::
              abs ((S[__a][__b] -
                    S_reference[__a][__b]) / (S_reference[__a][__b])) > 1)
          {
            std::cerr << "Mismatch detected in AVX implementation" << std::endl;
            std::cerr << "Variable S:" << std::endl;
            std::
              cerr << "seed=" << seed << ", __a=" << __a << ", __b=" << __b <<
              std::endl;
            std::cerr << "S AVX=  " << S[__a][__b] << std::endl;
            std::cerr << "S Reference=  " << S_reference[__a][__b] << std::endl;
            std::cerr << "S Rel Difference=  " << std::
              abs ((S[__a][__b] -
                    S_reference[__a][__b]) /
                   (S_reference[__a][__b])) << std::endl;
            std::cerr << "S Abs Difference=  " << std::abs (S[__a][__b] -
                                                            S_reference[__a]
                                                            [__b]) << std::endl;
            return 1;
          }
      for (int __a = 0; __a < 9; __a++)
        for (int __b = 0; __b < 16; __b++)
          if (std::
              abs ((V[__a][__b] -
                    V_reference[__a][__b]) / (V_reference[__a][__b])) > 1)
          {
            std::cerr << "Mismatch detected in AVX implementation" << std::endl;
            std::cerr << "Variable V:" << std::endl;
            std::
              cerr << "seed=" << seed << ", __a=" << __a << ", __b=" << __b <<
              std::endl;
            std::cerr << "V AVX=  " << V[__a][__b] << std::endl;
            std::cerr << "V Reference=  " << V_reference[__a][__b] << std::endl;
            std::cerr << "V Rel Difference=  " << std::
              abs ((V[__a][__b] -
                    V_reference[__a][__b]) /
                   (V_reference[__a][__b])) << std::endl;
            std::cerr << "V Abs Difference=  " << std::abs (V[__a][__b] -
                                                            V_reference[__a]
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
      typedef         T (&refArray3)[3][16];
      typedef         T (&refArray4)[9][16];
      for (int __a = 0; __a < 9; __a++)
        for (int __b = 0; __b < 16; __b++)
          U[__a][__b] = U_original[__a][__b];
      for (int __a = 0; __a < 3; __a++)
        for (int __b = 0; __b < 16; __b++)
          S[__a][__b] = S_original[__a][__b];
      for (int __a = 0; __a < 9; __a++)
        for (int __b = 0; __b < 16; __b++)
          V[__a][__b] = V_original[__a][__b];
      for (int i = 0; i < 16; i += 4)
      {
        refArray1       Ak = reinterpret_cast < refArray1 > (A[0][i]);
        refArray2       Uk = reinterpret_cast < refArray2 > (U[0][i]);
        refArray3       Sk = reinterpret_cast < refArray3 > (S[0][i]);
        refArray4       Vk = reinterpret_cast < refArray4 > (V[0][i]);
        Singular_Value_Decomposition < float32x4_t, float[16], int[16] > (Ak,
                                                                          Uk,
                                                                          Sk,
                                                                          Vk);
      }
      for (int __a = 0; __a < 9; __a++)
        for (int __b = 0; __b < 16; __b++)
          if (std::
              abs ((U[__a][__b] -
                    U_reference[__a][__b]) / (U_reference[__a][__b])) > 1)
          {
            std::cerr << "Mismatch detected in NEON implementation" << std::
              endl;
            std::cerr << "Variable U:" << std::endl;
            std::
              cerr << "seed=" << seed << ", __a=" << __a << ", __b=" << __b <<
              std::endl;
            std::cerr << "U NEON=  " << U[__a][__b] << std::endl;
            std::cerr << "U Reference=  " << U_reference[__a][__b] << std::endl;
            std::cerr << "U Rel Difference=  " << std::
              abs ((U[__a][__b] -
                    U_reference[__a][__b]) /
                   (U_reference[__a][__b])) << std::endl;
            std::cerr << "U Abs Difference=  " << std::abs (U[__a][__b] -
                                                            U_reference[__a]
                                                            [__b]) << std::endl;
            return 1;
          }
      for (int __a = 0; __a < 3; __a++)
        for (int __b = 0; __b < 16; __b++)
          if (std::
              abs ((S[__a][__b] -
                    S_reference[__a][__b]) / (S_reference[__a][__b])) > 1)
          {
            std::cerr << "Mismatch detected in NEON implementation" << std::
              endl;
            std::cerr << "Variable S:" << std::endl;
            std::
              cerr << "seed=" << seed << ", __a=" << __a << ", __b=" << __b <<
              std::endl;
            std::cerr << "S NEON=  " << S[__a][__b] << std::endl;
            std::cerr << "S Reference=  " << S_reference[__a][__b] << std::endl;
            std::cerr << "S Rel Difference=  " << std::
              abs ((S[__a][__b] -
                    S_reference[__a][__b]) /
                   (S_reference[__a][__b])) << std::endl;
            std::cerr << "S Abs Difference=  " << std::abs (S[__a][__b] -
                                                            S_reference[__a]
                                                            [__b]) << std::endl;
            return 1;
          }
      for (int __a = 0; __a < 9; __a++)
        for (int __b = 0; __b < 16; __b++)
          if (std::
              abs ((V[__a][__b] -
                    V_reference[__a][__b]) / (V_reference[__a][__b])) > 1)
          {
            std::cerr << "Mismatch detected in NEON implementation" << std::
              endl;
            std::cerr << "Variable V:" << std::endl;
            std::
              cerr << "seed=" << seed << ", __a=" << __a << ", __b=" << __b <<
              std::endl;
            std::cerr << "V NEON=  " << V[__a][__b] << std::endl;
            std::cerr << "V Reference=  " << V_reference[__a][__b] << std::endl;
            std::cerr << "V Rel Difference=  " << std::
              abs ((V[__a][__b] -
                    V_reference[__a][__b]) /
                   (V_reference[__a][__b])) << std::endl;
            std::cerr << "V Abs Difference=  " << std::abs (V[__a][__b] -
                                                            V_reference[__a]
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
      typedef         T (&refArray3)[3][16];
      typedef         T (&refArray4)[9][16];
      for (int __a = 0; __a < 9; __a++)
        for (int __b = 0; __b < 16; __b++)
          U[__a][__b] = U_original[__a][__b];
      for (int __a = 0; __a < 3; __a++)
        for (int __b = 0; __b < 16; __b++)
          S[__a][__b] = S_original[__a][__b];
      for (int __a = 0; __a < 9; __a++)
        for (int __b = 0; __b < 16; __b++)
          V[__a][__b] = V_original[__a][__b];
      for (int i = 0; i < 16; i += 16)
      {
        refArray1       Ak = reinterpret_cast < refArray1 > (A[0][i]);
        refArray2       Uk = reinterpret_cast < refArray2 > (U[0][i]);
        refArray3       Sk = reinterpret_cast < refArray3 > (S[0][i]);
        refArray4       Vk = reinterpret_cast < refArray4 > (V[0][i]);
        Singular_Value_Decomposition < __m512, float[16], int[16] > (Ak, Uk, Sk,
                                                                     Vk);
      }
      for (int __a = 0; __a < 9; __a++)
        for (int __b = 0; __b < 16; __b++)
          if (std::
              abs ((U[__a][__b] -
                    U_reference[__a][__b]) / (U_reference[__a][__b])) > 1)
          {
            std::cerr << "Mismatch detected in MIC implementation" << std::endl;
            std::cerr << "Variable U:" << std::endl;
            std::
              cerr << "seed=" << seed << ", __a=" << __a << ", __b=" << __b <<
              std::endl;
            std::cerr << "U MIC=  " << U[__a][__b] << std::endl;
            std::cerr << "U Reference=  " << U_reference[__a][__b] << std::endl;
            std::cerr << "U Rel Difference=  " << std::
              abs ((U[__a][__b] -
                    U_reference[__a][__b]) /
                   (U_reference[__a][__b])) << std::endl;
            std::cerr << "U Abs Difference=  " << std::abs (U[__a][__b] -
                                                            U_reference[__a]
                                                            [__b]) << std::endl;
            return 1;
          }
      for (int __a = 0; __a < 3; __a++)
        for (int __b = 0; __b < 16; __b++)
          if (std::
              abs ((S[__a][__b] -
                    S_reference[__a][__b]) / (S_reference[__a][__b])) > 1)
          {
            std::cerr << "Mismatch detected in MIC implementation" << std::endl;
            std::cerr << "Variable S:" << std::endl;
            std::
              cerr << "seed=" << seed << ", __a=" << __a << ", __b=" << __b <<
              std::endl;
            std::cerr << "S MIC=  " << S[__a][__b] << std::endl;
            std::cerr << "S Reference=  " << S_reference[__a][__b] << std::endl;
            std::cerr << "S Rel Difference=  " << std::
              abs ((S[__a][__b] -
                    S_reference[__a][__b]) /
                   (S_reference[__a][__b])) << std::endl;
            std::cerr << "S Abs Difference=  " << std::abs (S[__a][__b] -
                                                            S_reference[__a]
                                                            [__b]) << std::endl;
            return 1;
          }
      for (int __a = 0; __a < 9; __a++)
        for (int __b = 0; __b < 16; __b++)
          if (std::
              abs ((V[__a][__b] -
                    V_reference[__a][__b]) / (V_reference[__a][__b])) > 1)
          {
            std::cerr << "Mismatch detected in MIC implementation" << std::endl;
            std::cerr << "Variable V:" << std::endl;
            std::
              cerr << "seed=" << seed << ", __a=" << __a << ", __b=" << __b <<
              std::endl;
            std::cerr << "V MIC=  " << V[__a][__b] << std::endl;
            std::cerr << "V Reference=  " << V_reference[__a][__b] << std::endl;
            std::cerr << "V Rel Difference=  " << std::
              abs ((V[__a][__b] -
                    V_reference[__a][__b]) /
                   (V_reference[__a][__b])) << std::endl;
            std::cerr << "V Abs Difference=  " << std::abs (V[__a][__b] -
                                                            V_reference[__a]
                                                            [__b]) << std::endl;
            return 1;
          }

    }
#endif

  }



  std::cout << "SIMD check successful!" << std::endl;

  return 0;

}

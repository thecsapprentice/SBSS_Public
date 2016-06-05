
#include <cstdlib>
#include <iostream>
#include "KernelCommon.h"

struct NEOHOOKEAN_TAG;
struct COROTATED_TAG;

#include "Muscle_Differential.h"
#include "Muscle_Differential_Reference.h"

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
    std::cout << "Running SIMD Test for Muscle_Differential " << std::endl;


//=======================================================
//
//        DEFINE ALL VARIABLES USED BY KERNEL
//
//=======================================================

    T               dP_fiber[9][16] __attribute__ ((aligned (64)));
    T               dP_fiber_reference[9][16] __attribute__ ((aligned (64)));
    T               dP_fiber_original[9][16] __attribute__ ((aligned (64)));
    T               dF[9][16] __attribute__ ((aligned (64)));
    T               fiber[3][16] __attribute__ ((aligned (64)));
    T               Ffiber[3][16] __attribute__ ((aligned (64)));
    T               c1[16] __attribute__ ((aligned (64)));
    T               c2[16] __attribute__ ((aligned (64)));


    for (int __a = 0; __a < 9; __a++)
      for (int __b = 0; __b < 16; __b++)
      {
        dP_fiber_original[__a][__b] = Get_Random < float >();
        dP_fiber[__a][__b] = dP_fiber_original[__a][__b];
        dP_fiber_reference[__a][__b] = dP_fiber_original[__a][__b];
      }
    for (int __a = 0; __a < 9; __a++)
      for (int __b = 0; __b < 16; __b++)
        dF[__a][__b] = Get_Random < float >();
    for (int __a = 0; __a < 3; __a++)
      for (int __b = 0; __b < 16; __b++)
        fiber[__a][__b] = Get_Random < float >();
    for (int __a = 0; __a < 3; __a++)
      for (int __b = 0; __b < 16; __b++)
        Ffiber[__a][__b] = Get_Random < float >();
    for (int __a = 0; __a < 16; __a++)
      c1[__a] = Get_Random < float >();
    for (int __a = 0; __a < 16; __a++)
      c2[__a] = Get_Random < float >();

//=======================================================
//
//             COMPUTE REFERENCE RESULTS
//
//=======================================================

    T               __mdP_fiber[9] __attribute__ ((aligned (4)));
    T               __mdP_fiber_reference[9] __attribute__ ((aligned (4)));
    T               __mdP_fiber_original[9] __attribute__ ((aligned (4)));
    T               __mdF[9] __attribute__ ((aligned (4)));
    T               __mfiber[3] __attribute__ ((aligned (4)));
    T               __mFfiber[3] __attribute__ ((aligned (4)));
    T __mc1 __attribute__ ((aligned (4)));
    T __mc2 __attribute__ ((aligned (4)));
    for (int k = 0; k < 16; k++)
    {
      for (int __a = 0; __a < 9; __a++)
        __mdP_fiber_reference[__a] = dP_fiber_reference[__a][k];
      for (int __a = 0; __a < 9; __a++)
        __mdF[__a] = dF[__a][k];
      for (int __a = 0; __a < 3; __a++)
        __mfiber[__a] = fiber[__a][k];
      for (int __a = 0; __a < 3; __a++)
        __mFfiber[__a] = Ffiber[__a][k];
      __mc1 = c1[k];
      __mc2 = c2[k];
      Muscle_Differential < float, float, int >(__mdP_fiber_reference, __mdF,
                                                __mfiber, __mFfiber, __mc1,
                                                __mc2);
      for (int __a = 0; __a < 9; __a++)
        dP_fiber_reference[__a][k] = __mdP_fiber_reference[__a];
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
      typedef         T (&refArray4)[3][16];
      typedef         T (&refArray5)[16];
      typedef         T (&refArray6)[16];
      for (int __a = 0; __a < 9; __a++)
        for (int __b = 0; __b < 16; __b++)
          dP_fiber[__a][__b] = dP_fiber_original[__a][__b];
      for (int i = 0; i < 16; i += 1)
      {
        refArray1       dP_fiberk =
          reinterpret_cast < refArray1 > (dP_fiber[0][i]);
        refArray2       dFk = reinterpret_cast < refArray2 > (dF[0][i]);
        refArray3       fiberk = reinterpret_cast < refArray3 > (fiber[0][i]);
        refArray4       Ffiberk = reinterpret_cast < refArray4 > (Ffiber[0][i]);
        refArray5       c1k = reinterpret_cast < refArray5 > (c1[i]);
        refArray6       c2k = reinterpret_cast < refArray6 > (c2[i]);
        Muscle_Differential < float, float[16], int[16] > (dP_fiberk, dFk,
                                                           fiberk, Ffiberk, c1k,
                                                           c2k);
      }
      for (int __a = 0; __a < 9; __a++)
        for (int __b = 0; __b < 16; __b++)
          if (std::
              abs ((dP_fiber[__a][__b] -
                    dP_fiber_reference[__a][__b]) /
                   (dP_fiber_reference[__a][__b])) > 1)
          {
            std::cerr << "Mismatch detected in SCALAR implementation" << std::
              endl;
            std::cerr << "Variable dP_fiber:" << std::endl;
            std::
              cerr << "seed=" << seed << ", __a=" << __a << ", __b=" << __b <<
              std::endl;
            std::cerr << "dP_fiber SCALAR=  " << dP_fiber[__a][__b] << std::
              endl;
            std::
              cerr << "dP_fiber Reference=  " << dP_fiber_reference[__a][__b] <<
              std::endl;
            std::cerr << "dP_fiber Rel Difference=  " << std::
              abs ((dP_fiber[__a][__b] -
                    dP_fiber_reference[__a][__b]) /
                   (dP_fiber_reference[__a][__b])) << std::endl;
            std::cerr << "dP_fiber Abs Difference=  " << std::
              abs (dP_fiber[__a][__b] -
                   dP_fiber_reference[__a][__b]) << std::endl;
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
      typedef         T (&refArray4)[3][16];
      typedef         T (&refArray5)[16];
      typedef         T (&refArray6)[16];
      for (int __a = 0; __a < 9; __a++)
        for (int __b = 0; __b < 16; __b++)
          dP_fiber[__a][__b] = dP_fiber_original[__a][__b];
      for (int i = 0; i < 16; i += 4)
      {
        refArray1       dP_fiberk =
          reinterpret_cast < refArray1 > (dP_fiber[0][i]);
        refArray2       dFk = reinterpret_cast < refArray2 > (dF[0][i]);
        refArray3       fiberk = reinterpret_cast < refArray3 > (fiber[0][i]);
        refArray4       Ffiberk = reinterpret_cast < refArray4 > (Ffiber[0][i]);
        refArray5       c1k = reinterpret_cast < refArray5 > (c1[i]);
        refArray6       c2k = reinterpret_cast < refArray6 > (c2[i]);
        Muscle_Differential < __m128, float[16], int[16] > (dP_fiberk, dFk,
                                                            fiberk, Ffiberk,
                                                            c1k, c2k);
      }
      for (int __a = 0; __a < 9; __a++)
        for (int __b = 0; __b < 16; __b++)
          if (std::
              abs ((dP_fiber[__a][__b] -
                    dP_fiber_reference[__a][__b]) /
                   (dP_fiber_reference[__a][__b])) > 1)
          {
            std::cerr << "Mismatch detected in SSE implementation" << std::endl;
            std::cerr << "Variable dP_fiber:" << std::endl;
            std::
              cerr << "seed=" << seed << ", __a=" << __a << ", __b=" << __b <<
              std::endl;
            std::cerr << "dP_fiber SSE=  " << dP_fiber[__a][__b] << std::endl;
            std::
              cerr << "dP_fiber Reference=  " << dP_fiber_reference[__a][__b] <<
              std::endl;
            std::cerr << "dP_fiber Rel Difference=  " << std::
              abs ((dP_fiber[__a][__b] -
                    dP_fiber_reference[__a][__b]) /
                   (dP_fiber_reference[__a][__b])) << std::endl;
            std::cerr << "dP_fiber Abs Difference=  " << std::
              abs (dP_fiber[__a][__b] -
                   dP_fiber_reference[__a][__b]) << std::endl;
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
      typedef         T (&refArray4)[3][16];
      typedef         T (&refArray5)[16];
      typedef         T (&refArray6)[16];
      for (int __a = 0; __a < 9; __a++)
        for (int __b = 0; __b < 16; __b++)
          dP_fiber[__a][__b] = dP_fiber_original[__a][__b];
      for (int i = 0; i < 16; i += 8)
      {
        refArray1       dP_fiberk =
          reinterpret_cast < refArray1 > (dP_fiber[0][i]);
        refArray2       dFk = reinterpret_cast < refArray2 > (dF[0][i]);
        refArray3       fiberk = reinterpret_cast < refArray3 > (fiber[0][i]);
        refArray4       Ffiberk = reinterpret_cast < refArray4 > (Ffiber[0][i]);
        refArray5       c1k = reinterpret_cast < refArray5 > (c1[i]);
        refArray6       c2k = reinterpret_cast < refArray6 > (c2[i]);
        Muscle_Differential < __m256, float[16], int[16] > (dP_fiberk, dFk,
                                                            fiberk, Ffiberk,
                                                            c1k, c2k);
      }
      for (int __a = 0; __a < 9; __a++)
        for (int __b = 0; __b < 16; __b++)
          if (std::
              abs ((dP_fiber[__a][__b] -
                    dP_fiber_reference[__a][__b]) /
                   (dP_fiber_reference[__a][__b])) > 1)
          {
            std::cerr << "Mismatch detected in AVX implementation" << std::endl;
            std::cerr << "Variable dP_fiber:" << std::endl;
            std::
              cerr << "seed=" << seed << ", __a=" << __a << ", __b=" << __b <<
              std::endl;
            std::cerr << "dP_fiber AVX=  " << dP_fiber[__a][__b] << std::endl;
            std::
              cerr << "dP_fiber Reference=  " << dP_fiber_reference[__a][__b] <<
              std::endl;
            std::cerr << "dP_fiber Rel Difference=  " << std::
              abs ((dP_fiber[__a][__b] -
                    dP_fiber_reference[__a][__b]) /
                   (dP_fiber_reference[__a][__b])) << std::endl;
            std::cerr << "dP_fiber Abs Difference=  " << std::
              abs (dP_fiber[__a][__b] -
                   dP_fiber_reference[__a][__b]) << std::endl;
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
      typedef         T (&refArray4)[3][16];
      typedef         T (&refArray5)[16];
      typedef         T (&refArray6)[16];
      for (int __a = 0; __a < 9; __a++)
        for (int __b = 0; __b < 16; __b++)
          dP_fiber[__a][__b] = dP_fiber_original[__a][__b];
      for (int i = 0; i < 16; i += 4)
      {
        refArray1       dP_fiberk =
          reinterpret_cast < refArray1 > (dP_fiber[0][i]);
        refArray2       dFk = reinterpret_cast < refArray2 > (dF[0][i]);
        refArray3       fiberk = reinterpret_cast < refArray3 > (fiber[0][i]);
        refArray4       Ffiberk = reinterpret_cast < refArray4 > (Ffiber[0][i]);
        refArray5       c1k = reinterpret_cast < refArray5 > (c1[i]);
        refArray6       c2k = reinterpret_cast < refArray6 > (c2[i]);
        Muscle_Differential < float32x4_t, float[16], int[16] > (dP_fiberk, dFk,
                                                                 fiberk,
                                                                 Ffiberk, c1k,
                                                                 c2k);
      }
      for (int __a = 0; __a < 9; __a++)
        for (int __b = 0; __b < 16; __b++)
          if (std::
              abs ((dP_fiber[__a][__b] -
                    dP_fiber_reference[__a][__b]) /
                   (dP_fiber_reference[__a][__b])) > 1)
          {
            std::cerr << "Mismatch detected in NEON implementation" << std::
              endl;
            std::cerr << "Variable dP_fiber:" << std::endl;
            std::
              cerr << "seed=" << seed << ", __a=" << __a << ", __b=" << __b <<
              std::endl;
            std::cerr << "dP_fiber NEON=  " << dP_fiber[__a][__b] << std::endl;
            std::
              cerr << "dP_fiber Reference=  " << dP_fiber_reference[__a][__b] <<
              std::endl;
            std::cerr << "dP_fiber Rel Difference=  " << std::
              abs ((dP_fiber[__a][__b] -
                    dP_fiber_reference[__a][__b]) /
                   (dP_fiber_reference[__a][__b])) << std::endl;
            std::cerr << "dP_fiber Abs Difference=  " << std::
              abs (dP_fiber[__a][__b] -
                   dP_fiber_reference[__a][__b]) << std::endl;
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
      typedef         T (&refArray4)[3][16];
      typedef         T (&refArray5)[16];
      typedef         T (&refArray6)[16];
      for (int __a = 0; __a < 9; __a++)
        for (int __b = 0; __b < 16; __b++)
          dP_fiber[__a][__b] = dP_fiber_original[__a][__b];
      for (int i = 0; i < 16; i += 16)
      {
        refArray1       dP_fiberk =
          reinterpret_cast < refArray1 > (dP_fiber[0][i]);
        refArray2       dFk = reinterpret_cast < refArray2 > (dF[0][i]);
        refArray3       fiberk = reinterpret_cast < refArray3 > (fiber[0][i]);
        refArray4       Ffiberk = reinterpret_cast < refArray4 > (Ffiber[0][i]);
        refArray5       c1k = reinterpret_cast < refArray5 > (c1[i]);
        refArray6       c2k = reinterpret_cast < refArray6 > (c2[i]);
        Muscle_Differential < __m512, float[16], int[16] > (dP_fiberk, dFk,
                                                            fiberk, Ffiberk,
                                                            c1k, c2k);
      }
      for (int __a = 0; __a < 9; __a++)
        for (int __b = 0; __b < 16; __b++)
          if (std::
              abs ((dP_fiber[__a][__b] -
                    dP_fiber_reference[__a][__b]) /
                   (dP_fiber_reference[__a][__b])) > 1)
          {
            std::cerr << "Mismatch detected in MIC implementation" << std::endl;
            std::cerr << "Variable dP_fiber:" << std::endl;
            std::
              cerr << "seed=" << seed << ", __a=" << __a << ", __b=" << __b <<
              std::endl;
            std::cerr << "dP_fiber MIC=  " << dP_fiber[__a][__b] << std::endl;
            std::
              cerr << "dP_fiber Reference=  " << dP_fiber_reference[__a][__b] <<
              std::endl;
            std::cerr << "dP_fiber Rel Difference=  " << std::
              abs ((dP_fiber[__a][__b] -
                    dP_fiber_reference[__a][__b]) /
                   (dP_fiber_reference[__a][__b])) << std::endl;
            std::cerr << "dP_fiber Abs Difference=  " << std::
              abs (dP_fiber[__a][__b] -
                   dP_fiber_reference[__a][__b]) << std::endl;
            return 1;
          }

    }
#endif

  }



  std::cout << "SIMD check successful!" << std::endl;

  return 0;

}

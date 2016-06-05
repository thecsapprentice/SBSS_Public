
#include <cstdlib>
#include <iostream>
#include "KernelCommon.h"

struct NEOHOOKEAN_TAG;
struct COROTATED_TAG;

#include "Muscle_Update_Position_Based_State.h"
#include "Muscle_Update_Position_Based_State_Reference.h"

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
      cout << "Running SIMD Test for Muscle_Update_Position_Based_State " <<
      std::endl;


//=======================================================
//
//        DEFINE ALL VARIABLES USED BY KERNEL
//
//=======================================================

    T               u[3][8][16] __attribute__ ((aligned (64)));
    int             muscle_id[16] __attribute__ ((aligned (64)));
    T               fiber[3][16] __attribute__ ((aligned (64)));
    T               density[16] __attribute__ ((aligned (64)));
    T               one_over_h[16] __attribute__ ((aligned (64)));
    T               c1[16] __attribute__ ((aligned (64)));
    T               c1_reference[16] __attribute__ ((aligned (64)));
    T               c1_original[16] __attribute__ ((aligned (64)));
    T               c2[16] __attribute__ ((aligned (64)));
    T               c2_reference[16] __attribute__ ((aligned (64)));
    T               c2_original[16] __attribute__ ((aligned (64)));
    T               F_fiber[3][16] __attribute__ ((aligned (64)));
    T               F_fiber_reference[3][16] __attribute__ ((aligned (64)));
    T               F_fiber_original[3][16] __attribute__ ((aligned (64)));
    float          *activations __attribute__ ((aligned (64)));
    float          *fiber_max_stresses __attribute__ ((aligned (64)));


    for (int __a = 0; __a < 3; __a++)
      for (int __b = 0; __b < 8; __b++)
        for (int __c = 0; __c < 16; __c++)
          u[__a][__b][__c] = Get_Random < float >();
    for (int __a = 0; __a < 16; __a++)
      muscle_id[__a] = Get_Random < int >(1, 99);
    for (int __a = 0; __a < 3; __a++)
      for (int __b = 0; __b < 16; __b++)
        fiber[__a][__b] = Get_Random < float >();
    for (int __a = 0; __a < 16; __a++)
      density[__a] = Get_Random < float >();
    for (int __a = 0; __a < 16; __a++)
      one_over_h[__a] = Get_Random < float >();
    for (int __a = 0; __a < 16; __a++)
    {
      c1_original[__a] = Get_Random < float >();
      c1[__a] = c1_original[__a];
      c1_reference[__a] = c1_original[__a];
    }
    for (int __a = 0; __a < 16; __a++)
    {
      c2_original[__a] = Get_Random < float >();
      c2[__a] = c2_original[__a];
      c2_reference[__a] = c2_original[__a];
    }
    for (int __a = 0; __a < 3; __a++)
      for (int __b = 0; __b < 16; __b++)
      {
        F_fiber_original[__a][__b] = Get_Random < float >();
        F_fiber[__a][__b] = F_fiber_original[__a][__b];
        F_fiber_reference[__a][__b] = F_fiber_original[__a][__b];
      }

    activations = (float* )(_mm_malloc(100 * sizeof (T), 64));
    for (int __x__ = 0; __x__ < 100; __x__++)
      activations[__x__] = Get_Random < float >();;
    fiber_max_stresses = (float* )(_mm_malloc(100 * sizeof (T), 64));
    for (int __x__ = 0; __x__ < 100; __x__++)
      fiber_max_stresses[__x__] = Get_Random < float >();;

//=======================================================
//
//             COMPUTE REFERENCE RESULTS
//
//=======================================================

    T               __mu[3][8] __attribute__ ((aligned (4)));
    int __mmuscle_id __attribute__ ((aligned (4)));
    T               __mfiber[3] __attribute__ ((aligned (4)));
    T __mdensity __attribute__ ((aligned (4)));
    T __mone_over_h __attribute__ ((aligned (4)));
    T __mc1 __attribute__ ((aligned (4)));
    T __mc1_reference __attribute__ ((aligned (4)));
    T __mc1_original __attribute__ ((aligned (4)));
    T __mc2 __attribute__ ((aligned (4)));
    T __mc2_reference __attribute__ ((aligned (4)));
    T __mc2_original __attribute__ ((aligned (4)));
    T               __mF_fiber[3] __attribute__ ((aligned (4)));
    T               __mF_fiber_reference[3] __attribute__ ((aligned (4)));
    T               __mF_fiber_original[3] __attribute__ ((aligned (4)));
    float          *__mactivations __attribute__ ((aligned (4)));
    float          *__mfiber_max_stresses __attribute__ ((aligned (4)));
    for (int k = 0; k < 16; k++)
    {
      for (int __a = 0; __a < 3; __a++)
        for (int __b = 0; __b < 8; __b++)
          __mu[__a][__b] = u[__a][__b][k];
      __mmuscle_id = muscle_id[k];
      for (int __a = 0; __a < 3; __a++)
        __mfiber[__a] = fiber[__a][k];
      __mdensity = density[k];
      __mone_over_h = one_over_h[k];
      __mc1_reference = c1_reference[k];

      __mc2_reference = c2_reference[k];
      for (int __a = 0; __a < 3; __a++)
        __mF_fiber_reference[__a] = F_fiber_reference[__a][k];

      __mactivations = activations;
      __mfiber_max_stresses = fiber_max_stresses;
      Muscle_Update_Position_Based_State < float, float, int >(__mu,
                                                               __mmuscle_id,
                                                               __mfiber,
                                                               __mdensity,
                                                               __mone_over_h,
                                                               __mc1_reference,
                                                               __mc2_reference,
                                                               __mF_fiber_reference,
                                                               __mactivations,
                                                               __mfiber_max_stresses);
      c1_reference[k] = __mc1_reference;
      c2_reference[k] = __mc2_reference;
      for (int __a = 0; __a < 3; __a++)
        F_fiber_reference[__a][k] = __mF_fiber_reference[__a];
    }

//=======================================================
//
//               COMPUTE SCALAR RESULTS
//
//=======================================================

    {
      typedef         T (&refArray1)[3][8][16];
      typedef int     (&refArray2)[16];
      typedef         T (&refArray3)[3][16];
      typedef         T (&refArray4)[16];
      typedef         T (&refArray5)[16];
      typedef         T (&refArray6)[16];
      typedef         T (&refArray7)[16];
      typedef         T (&refArray8)[3][16];
      typedef float  *(&refArray9);
      typedef float  *(&refArray10);
      for (int __a = 0; __a < 16; __a++)
        c1[__a] = c1_original[__a];
      for (int __a = 0; __a < 16; __a++)
        c2[__a] = c2_original[__a];
      for (int __a = 0; __a < 3; __a++)
        for (int __b = 0; __b < 16; __b++)
          F_fiber[__a][__b] = F_fiber_original[__a][__b];
      for (int i = 0; i < 16; i += 1)
      {
        refArray1       uk = reinterpret_cast < refArray1 > (u[0][0][i]);
        refArray2       muscle_idk =
          reinterpret_cast < refArray2 > (muscle_id[i]);
        refArray3       fiberk = reinterpret_cast < refArray3 > (fiber[0][i]);
        refArray4       densityk = reinterpret_cast < refArray4 > (density[i]);
        refArray5       one_over_hk =
          reinterpret_cast < refArray5 > (one_over_h[i]);
        refArray6       c1k = reinterpret_cast < refArray6 > (c1[i]);
        refArray7       c2k = reinterpret_cast < refArray7 > (c2[i]);
        refArray8       F_fiberk =
          reinterpret_cast < refArray8 > (F_fiber[0][i]);
        refArray9       activationsk =
          reinterpret_cast < refArray9 > (activations);
        refArray10      fiber_max_stressesk =
          reinterpret_cast < refArray10 > (fiber_max_stresses);
        Muscle_Update_Position_Based_State < float, float[16], int[16] > (uk,
                                                                          muscle_idk,
                                                                          fiberk,
                                                                          densityk,
                                                                          one_over_hk,
                                                                          c1k,
                                                                          c2k,
                                                                          F_fiberk,
                                                                          activationsk,
                                                                          fiber_max_stressesk);
      }
      for (int __a = 0; __a < 16; __a++)
        if (std::abs ((c1[__a] - c1_reference[__a]) / (c1_reference[__a])) > 1)
        {
          std::cerr << "Mismatch detected in SCALAR implementation" << std::
            endl;
          std::cerr << "Variable c1:" << std::endl;
          std::cerr << "seed=" << seed << ", __a=" << __a << std::endl;
          std::cerr << "c1 SCALAR=  " << c1[__a] << std::endl;
          std::cerr << "c1 Reference=  " << c1_reference[__a] << std::endl;
          std::cerr << "c1 Rel Difference=  " << std::
            abs ((c1[__a] -
                  c1_reference[__a]) / (c1_reference[__a])) << std::endl;
          std::cerr << "c1 Abs Difference=  " << std::abs (c1[__a] -
                                                           c1_reference[__a]) <<
            std::endl;
          return 1;
        }
      for (int __a = 0; __a < 16; __a++)
        if (std::abs ((c2[__a] - c2_reference[__a]) / (c2_reference[__a])) > 1)
        {
          std::cerr << "Mismatch detected in SCALAR implementation" << std::
            endl;
          std::cerr << "Variable c2:" << std::endl;
          std::cerr << "seed=" << seed << ", __a=" << __a << std::endl;
          std::cerr << "c2 SCALAR=  " << c2[__a] << std::endl;
          std::cerr << "c2 Reference=  " << c2_reference[__a] << std::endl;
          std::cerr << "c2 Rel Difference=  " << std::
            abs ((c2[__a] -
                  c2_reference[__a]) / (c2_reference[__a])) << std::endl;
          std::cerr << "c2 Abs Difference=  " << std::abs (c2[__a] -
                                                           c2_reference[__a]) <<
            std::endl;
          return 1;
        }
      for (int __a = 0; __a < 3; __a++)
        for (int __b = 0; __b < 16; __b++)
          if (std::
              abs ((F_fiber[__a][__b] -
                    F_fiber_reference[__a][__b]) /
                   (F_fiber_reference[__a][__b])) > 1)
          {
            std::cerr << "Mismatch detected in SCALAR implementation" << std::
              endl;
            std::cerr << "Variable F_fiber:" << std::endl;
            std::
              cerr << "seed=" << seed << ", __a=" << __a << ", __b=" << __b <<
              std::endl;
            std::cerr << "F_fiber SCALAR=  " << F_fiber[__a][__b] << std::endl;
            std::
              cerr << "F_fiber Reference=  " << F_fiber_reference[__a][__b] <<
              std::endl;
            std::cerr << "F_fiber Rel Difference=  " << std::
              abs ((F_fiber[__a][__b] -
                    F_fiber_reference[__a][__b]) /
                   (F_fiber_reference[__a][__b])) << std::endl;
            std::cerr << "F_fiber Abs Difference=  " << std::
              abs (F_fiber[__a][__b] -
                   F_fiber_reference[__a][__b]) << std::endl;
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
      typedef int     (&refArray2)[16];
      typedef         T (&refArray3)[3][16];
      typedef         T (&refArray4)[16];
      typedef         T (&refArray5)[16];
      typedef         T (&refArray6)[16];
      typedef         T (&refArray7)[16];
      typedef         T (&refArray8)[3][16];
      typedef float  *(&refArray9);
      typedef float  *(&refArray10);
      for (int __a = 0; __a < 16; __a++)
        c1[__a] = c1_original[__a];
      for (int __a = 0; __a < 16; __a++)
        c2[__a] = c2_original[__a];
      for (int __a = 0; __a < 3; __a++)
        for (int __b = 0; __b < 16; __b++)
          F_fiber[__a][__b] = F_fiber_original[__a][__b];
      for (int i = 0; i < 16; i += 4)
      {
        refArray1       uk = reinterpret_cast < refArray1 > (u[0][0][i]);
        refArray2       muscle_idk =
          reinterpret_cast < refArray2 > (muscle_id[i]);
        refArray3       fiberk = reinterpret_cast < refArray3 > (fiber[0][i]);
        refArray4       densityk = reinterpret_cast < refArray4 > (density[i]);
        refArray5       one_over_hk =
          reinterpret_cast < refArray5 > (one_over_h[i]);
        refArray6       c1k = reinterpret_cast < refArray6 > (c1[i]);
        refArray7       c2k = reinterpret_cast < refArray7 > (c2[i]);
        refArray8       F_fiberk =
          reinterpret_cast < refArray8 > (F_fiber[0][i]);
        refArray9       activationsk =
          reinterpret_cast < refArray9 > (activations);
        refArray10      fiber_max_stressesk =
          reinterpret_cast < refArray10 > (fiber_max_stresses);
        Muscle_Update_Position_Based_State < __m128, float[16], int[16] > (uk,
                                                                           muscle_idk,
                                                                           fiberk,
                                                                           densityk,
                                                                           one_over_hk,
                                                                           c1k,
                                                                           c2k,
                                                                           F_fiberk,
                                                                           activationsk,
                                                                           fiber_max_stressesk);
      }
      for (int __a = 0; __a < 16; __a++)
        if (std::abs ((c1[__a] - c1_reference[__a]) / (c1_reference[__a])) > 1)
        {
          std::cerr << "Mismatch detected in SSE implementation" << std::endl;
          std::cerr << "Variable c1:" << std::endl;
          std::cerr << "seed=" << seed << ", __a=" << __a << std::endl;
          std::cerr << "c1 SSE=  " << c1[__a] << std::endl;
          std::cerr << "c1 Reference=  " << c1_reference[__a] << std::endl;
          std::cerr << "c1 Rel Difference=  " << std::
            abs ((c1[__a] -
                  c1_reference[__a]) / (c1_reference[__a])) << std::endl;
          std::cerr << "c1 Abs Difference=  " << std::abs (c1[__a] -
                                                           c1_reference[__a]) <<
            std::endl;
          return 1;
        }
      for (int __a = 0; __a < 16; __a++)
        if (std::abs ((c2[__a] - c2_reference[__a]) / (c2_reference[__a])) > 1)
        {
          std::cerr << "Mismatch detected in SSE implementation" << std::endl;
          std::cerr << "Variable c2:" << std::endl;
          std::cerr << "seed=" << seed << ", __a=" << __a << std::endl;
          std::cerr << "c2 SSE=  " << c2[__a] << std::endl;
          std::cerr << "c2 Reference=  " << c2_reference[__a] << std::endl;
          std::cerr << "c2 Rel Difference=  " << std::
            abs ((c2[__a] -
                  c2_reference[__a]) / (c2_reference[__a])) << std::endl;
          std::cerr << "c2 Abs Difference=  " << std::abs (c2[__a] -
                                                           c2_reference[__a]) <<
            std::endl;
          return 1;
        }
      for (int __a = 0; __a < 3; __a++)
        for (int __b = 0; __b < 16; __b++)
          if (std::
              abs ((F_fiber[__a][__b] -
                    F_fiber_reference[__a][__b]) /
                   (F_fiber_reference[__a][__b])) > 1)
          {
            std::cerr << "Mismatch detected in SSE implementation" << std::endl;
            std::cerr << "Variable F_fiber:" << std::endl;
            std::
              cerr << "seed=" << seed << ", __a=" << __a << ", __b=" << __b <<
              std::endl;
            std::cerr << "F_fiber SSE=  " << F_fiber[__a][__b] << std::endl;
            std::
              cerr << "F_fiber Reference=  " << F_fiber_reference[__a][__b] <<
              std::endl;
            std::cerr << "F_fiber Rel Difference=  " << std::
              abs ((F_fiber[__a][__b] -
                    F_fiber_reference[__a][__b]) /
                   (F_fiber_reference[__a][__b])) << std::endl;
            std::cerr << "F_fiber Abs Difference=  " << std::
              abs (F_fiber[__a][__b] -
                   F_fiber_reference[__a][__b]) << std::endl;
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
      typedef int     (&refArray2)[16];
      typedef         T (&refArray3)[3][16];
      typedef         T (&refArray4)[16];
      typedef         T (&refArray5)[16];
      typedef         T (&refArray6)[16];
      typedef         T (&refArray7)[16];
      typedef         T (&refArray8)[3][16];
      typedef float  *(&refArray9);
      typedef float  *(&refArray10);
      for (int __a = 0; __a < 16; __a++)
        c1[__a] = c1_original[__a];
      for (int __a = 0; __a < 16; __a++)
        c2[__a] = c2_original[__a];
      for (int __a = 0; __a < 3; __a++)
        for (int __b = 0; __b < 16; __b++)
          F_fiber[__a][__b] = F_fiber_original[__a][__b];
      for (int i = 0; i < 16; i += 8)
      {
        refArray1       uk = reinterpret_cast < refArray1 > (u[0][0][i]);
        refArray2       muscle_idk =
          reinterpret_cast < refArray2 > (muscle_id[i]);
        refArray3       fiberk = reinterpret_cast < refArray3 > (fiber[0][i]);
        refArray4       densityk = reinterpret_cast < refArray4 > (density[i]);
        refArray5       one_over_hk =
          reinterpret_cast < refArray5 > (one_over_h[i]);
        refArray6       c1k = reinterpret_cast < refArray6 > (c1[i]);
        refArray7       c2k = reinterpret_cast < refArray7 > (c2[i]);
        refArray8       F_fiberk =
          reinterpret_cast < refArray8 > (F_fiber[0][i]);
        refArray9       activationsk =
          reinterpret_cast < refArray9 > (activations);
        refArray10      fiber_max_stressesk =
          reinterpret_cast < refArray10 > (fiber_max_stresses);
        Muscle_Update_Position_Based_State < __m256, float[16], int[16] > (uk,
                                                                           muscle_idk,
                                                                           fiberk,
                                                                           densityk,
                                                                           one_over_hk,
                                                                           c1k,
                                                                           c2k,
                                                                           F_fiberk,
                                                                           activationsk,
                                                                           fiber_max_stressesk);
      }
      for (int __a = 0; __a < 16; __a++)
        if (std::abs ((c1[__a] - c1_reference[__a]) / (c1_reference[__a])) > 1)
        {
          std::cerr << "Mismatch detected in AVX implementation" << std::endl;
          std::cerr << "Variable c1:" << std::endl;
          std::cerr << "seed=" << seed << ", __a=" << __a << std::endl;
          std::cerr << "c1 AVX=  " << c1[__a] << std::endl;
          std::cerr << "c1 Reference=  " << c1_reference[__a] << std::endl;
          std::cerr << "c1 Rel Difference=  " << std::
            abs ((c1[__a] -
                  c1_reference[__a]) / (c1_reference[__a])) << std::endl;
          std::cerr << "c1 Abs Difference=  " << std::abs (c1[__a] -
                                                           c1_reference[__a]) <<
            std::endl;
          return 1;
        }
      for (int __a = 0; __a < 16; __a++)
        if (std::abs ((c2[__a] - c2_reference[__a]) / (c2_reference[__a])) > 1)
        {
          std::cerr << "Mismatch detected in AVX implementation" << std::endl;
          std::cerr << "Variable c2:" << std::endl;
          std::cerr << "seed=" << seed << ", __a=" << __a << std::endl;
          std::cerr << "c2 AVX=  " << c2[__a] << std::endl;
          std::cerr << "c2 Reference=  " << c2_reference[__a] << std::endl;
          std::cerr << "c2 Rel Difference=  " << std::
            abs ((c2[__a] -
                  c2_reference[__a]) / (c2_reference[__a])) << std::endl;
          std::cerr << "c2 Abs Difference=  " << std::abs (c2[__a] -
                                                           c2_reference[__a]) <<
            std::endl;
          return 1;
        }
      for (int __a = 0; __a < 3; __a++)
        for (int __b = 0; __b < 16; __b++)
          if (std::
              abs ((F_fiber[__a][__b] -
                    F_fiber_reference[__a][__b]) /
                   (F_fiber_reference[__a][__b])) > 1)
          {
            std::cerr << "Mismatch detected in AVX implementation" << std::endl;
            std::cerr << "Variable F_fiber:" << std::endl;
            std::
              cerr << "seed=" << seed << ", __a=" << __a << ", __b=" << __b <<
              std::endl;
            std::cerr << "F_fiber AVX=  " << F_fiber[__a][__b] << std::endl;
            std::
              cerr << "F_fiber Reference=  " << F_fiber_reference[__a][__b] <<
              std::endl;
            std::cerr << "F_fiber Rel Difference=  " << std::
              abs ((F_fiber[__a][__b] -
                    F_fiber_reference[__a][__b]) /
                   (F_fiber_reference[__a][__b])) << std::endl;
            std::cerr << "F_fiber Abs Difference=  " << std::
              abs (F_fiber[__a][__b] -
                   F_fiber_reference[__a][__b]) << std::endl;
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
      typedef int     (&refArray2)[16];
      typedef         T (&refArray3)[3][16];
      typedef         T (&refArray4)[16];
      typedef         T (&refArray5)[16];
      typedef         T (&refArray6)[16];
      typedef         T (&refArray7)[16];
      typedef         T (&refArray8)[3][16];
      typedef float  *(&refArray9);
      typedef float  *(&refArray10);
      for (int __a = 0; __a < 16; __a++)
        c1[__a] = c1_original[__a];
      for (int __a = 0; __a < 16; __a++)
        c2[__a] = c2_original[__a];
      for (int __a = 0; __a < 3; __a++)
        for (int __b = 0; __b < 16; __b++)
          F_fiber[__a][__b] = F_fiber_original[__a][__b];
      for (int i = 0; i < 16; i += 4)
      {
        refArray1       uk = reinterpret_cast < refArray1 > (u[0][0][i]);
        refArray2       muscle_idk =
          reinterpret_cast < refArray2 > (muscle_id[i]);
        refArray3       fiberk = reinterpret_cast < refArray3 > (fiber[0][i]);
        refArray4       densityk = reinterpret_cast < refArray4 > (density[i]);
        refArray5       one_over_hk =
          reinterpret_cast < refArray5 > (one_over_h[i]);
        refArray6       c1k = reinterpret_cast < refArray6 > (c1[i]);
        refArray7       c2k = reinterpret_cast < refArray7 > (c2[i]);
        refArray8       F_fiberk =
          reinterpret_cast < refArray8 > (F_fiber[0][i]);
        refArray9       activationsk =
          reinterpret_cast < refArray9 > (activations);
        refArray10      fiber_max_stressesk =
          reinterpret_cast < refArray10 > (fiber_max_stresses);
        Muscle_Update_Position_Based_State < float32x4_t, float[16],
          int[16] > (uk, muscle_idk, fiberk, densityk, one_over_hk, c1k, c2k,
                     F_fiberk, activationsk, fiber_max_stressesk);
      }
      for (int __a = 0; __a < 16; __a++)
        if (std::abs ((c1[__a] - c1_reference[__a]) / (c1_reference[__a])) > 1)
        {
          std::cerr << "Mismatch detected in NEON implementation" << std::endl;
          std::cerr << "Variable c1:" << std::endl;
          std::cerr << "seed=" << seed << ", __a=" << __a << std::endl;
          std::cerr << "c1 NEON=  " << c1[__a] << std::endl;
          std::cerr << "c1 Reference=  " << c1_reference[__a] << std::endl;
          std::cerr << "c1 Rel Difference=  " << std::
            abs ((c1[__a] -
                  c1_reference[__a]) / (c1_reference[__a])) << std::endl;
          std::cerr << "c1 Abs Difference=  " << std::abs (c1[__a] -
                                                           c1_reference[__a]) <<
            std::endl;
          return 1;
        }
      for (int __a = 0; __a < 16; __a++)
        if (std::abs ((c2[__a] - c2_reference[__a]) / (c2_reference[__a])) > 1)
        {
          std::cerr << "Mismatch detected in NEON implementation" << std::endl;
          std::cerr << "Variable c2:" << std::endl;
          std::cerr << "seed=" << seed << ", __a=" << __a << std::endl;
          std::cerr << "c2 NEON=  " << c2[__a] << std::endl;
          std::cerr << "c2 Reference=  " << c2_reference[__a] << std::endl;
          std::cerr << "c2 Rel Difference=  " << std::
            abs ((c2[__a] -
                  c2_reference[__a]) / (c2_reference[__a])) << std::endl;
          std::cerr << "c2 Abs Difference=  " << std::abs (c2[__a] -
                                                           c2_reference[__a]) <<
            std::endl;
          return 1;
        }
      for (int __a = 0; __a < 3; __a++)
        for (int __b = 0; __b < 16; __b++)
          if (std::
              abs ((F_fiber[__a][__b] -
                    F_fiber_reference[__a][__b]) /
                   (F_fiber_reference[__a][__b])) > 1)
          {
            std::cerr << "Mismatch detected in NEON implementation" << std::
              endl;
            std::cerr << "Variable F_fiber:" << std::endl;
            std::
              cerr << "seed=" << seed << ", __a=" << __a << ", __b=" << __b <<
              std::endl;
            std::cerr << "F_fiber NEON=  " << F_fiber[__a][__b] << std::endl;
            std::
              cerr << "F_fiber Reference=  " << F_fiber_reference[__a][__b] <<
              std::endl;
            std::cerr << "F_fiber Rel Difference=  " << std::
              abs ((F_fiber[__a][__b] -
                    F_fiber_reference[__a][__b]) /
                   (F_fiber_reference[__a][__b])) << std::endl;
            std::cerr << "F_fiber Abs Difference=  " << std::
              abs (F_fiber[__a][__b] -
                   F_fiber_reference[__a][__b]) << std::endl;
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
      typedef int     (&refArray2)[16];
      typedef         T (&refArray3)[3][16];
      typedef         T (&refArray4)[16];
      typedef         T (&refArray5)[16];
      typedef         T (&refArray6)[16];
      typedef         T (&refArray7)[16];
      typedef         T (&refArray8)[3][16];
      typedef float  *(&refArray9);
      typedef float  *(&refArray10);
      for (int __a = 0; __a < 16; __a++)
        c1[__a] = c1_original[__a];
      for (int __a = 0; __a < 16; __a++)
        c2[__a] = c2_original[__a];
      for (int __a = 0; __a < 3; __a++)
        for (int __b = 0; __b < 16; __b++)
          F_fiber[__a][__b] = F_fiber_original[__a][__b];
      for (int i = 0; i < 16; i += 16)
      {
        refArray1       uk = reinterpret_cast < refArray1 > (u[0][0][i]);
        refArray2       muscle_idk =
          reinterpret_cast < refArray2 > (muscle_id[i]);
        refArray3       fiberk = reinterpret_cast < refArray3 > (fiber[0][i]);
        refArray4       densityk = reinterpret_cast < refArray4 > (density[i]);
        refArray5       one_over_hk =
          reinterpret_cast < refArray5 > (one_over_h[i]);
        refArray6       c1k = reinterpret_cast < refArray6 > (c1[i]);
        refArray7       c2k = reinterpret_cast < refArray7 > (c2[i]);
        refArray8       F_fiberk =
          reinterpret_cast < refArray8 > (F_fiber[0][i]);
        refArray9       activationsk =
          reinterpret_cast < refArray9 > (activations);
        refArray10      fiber_max_stressesk =
          reinterpret_cast < refArray10 > (fiber_max_stresses);
        Muscle_Update_Position_Based_State < __m512, float[16], int[16] > (uk,
                                                                           muscle_idk,
                                                                           fiberk,
                                                                           densityk,
                                                                           one_over_hk,
                                                                           c1k,
                                                                           c2k,
                                                                           F_fiberk,
                                                                           activationsk,
                                                                           fiber_max_stressesk);
      }
      for (int __a = 0; __a < 16; __a++)
        if (std::abs ((c1[__a] - c1_reference[__a]) / (c1_reference[__a])) > 1)
        {
          std::cerr << "Mismatch detected in MIC implementation" << std::endl;
          std::cerr << "Variable c1:" << std::endl;
          std::cerr << "seed=" << seed << ", __a=" << __a << std::endl;
          std::cerr << "c1 MIC=  " << c1[__a] << std::endl;
          std::cerr << "c1 Reference=  " << c1_reference[__a] << std::endl;
          std::cerr << "c1 Rel Difference=  " << std::
            abs ((c1[__a] -
                  c1_reference[__a]) / (c1_reference[__a])) << std::endl;
          std::cerr << "c1 Abs Difference=  " << std::abs (c1[__a] -
                                                           c1_reference[__a]) <<
            std::endl;
          return 1;
        }
      for (int __a = 0; __a < 16; __a++)
        if (std::abs ((c2[__a] - c2_reference[__a]) / (c2_reference[__a])) > 1)
        {
          std::cerr << "Mismatch detected in MIC implementation" << std::endl;
          std::cerr << "Variable c2:" << std::endl;
          std::cerr << "seed=" << seed << ", __a=" << __a << std::endl;
          std::cerr << "c2 MIC=  " << c2[__a] << std::endl;
          std::cerr << "c2 Reference=  " << c2_reference[__a] << std::endl;
          std::cerr << "c2 Rel Difference=  " << std::
            abs ((c2[__a] -
                  c2_reference[__a]) / (c2_reference[__a])) << std::endl;
          std::cerr << "c2 Abs Difference=  " << std::abs (c2[__a] -
                                                           c2_reference[__a]) <<
            std::endl;
          return 1;
        }
      for (int __a = 0; __a < 3; __a++)
        for (int __b = 0; __b < 16; __b++)
          if (std::
              abs ((F_fiber[__a][__b] -
                    F_fiber_reference[__a][__b]) /
                   (F_fiber_reference[__a][__b])) > 1)
          {
            std::cerr << "Mismatch detected in MIC implementation" << std::endl;
            std::cerr << "Variable F_fiber:" << std::endl;
            std::
              cerr << "seed=" << seed << ", __a=" << __a << ", __b=" << __b <<
              std::endl;
            std::cerr << "F_fiber MIC=  " << F_fiber[__a][__b] << std::endl;
            std::
              cerr << "F_fiber Reference=  " << F_fiber_reference[__a][__b] <<
              std::endl;
            std::cerr << "F_fiber Rel Difference=  " << std::
              abs ((F_fiber[__a][__b] -
                    F_fiber_reference[__a][__b]) /
                   (F_fiber_reference[__a][__b])) << std::endl;
            std::cerr << "F_fiber Abs Difference=  " << std::
              abs (F_fiber[__a][__b] -
                   F_fiber_reference[__a][__b]) << std::endl;
            return 1;
          }

    }
#endif

  }



  std::cout << "SIMD check successful!" << std::endl;

  return 0;

}

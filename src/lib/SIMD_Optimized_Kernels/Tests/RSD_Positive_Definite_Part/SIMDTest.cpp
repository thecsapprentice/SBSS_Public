
#include <cstdlib>
#include <iostream>
#include "KernelCommon.h"

struct NEOHOOKEAN_TAG;
struct COROTATED_TAG;

#include "RSD_Positive_Definite_Part.h"
#include "RSD_Positive_Definite_Part_Reference.h"

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
    std::cout << "Running SIMD Test for RSD_Positive_Definite_Part " << std::
      endl;


//=======================================================
//
//        DEFINE ALL VARIABLES USED BY KERNEL
//
//=======================================================

    T               RSD[12][16] __attribute__ ((aligned (64)));
    T               RSDpd[12][16] __attribute__ ((aligned (64)));
    T               RSDpd_reference[12][16] __attribute__ ((aligned (64)));
    T               RSDpd_original[12][16] __attribute__ ((aligned (64)));


    for (int __a = 0; __a < 12; __a++)
      for (int __b = 0; __b < 16; __b++)
        RSD[__a][__b] = Get_Random < float >();
    for (int __a = 0; __a < 12; __a++)
      for (int __b = 0; __b < 16; __b++)
      {
        RSDpd_original[__a][__b] = Get_Random < float >();
        RSDpd[__a][__b] = RSDpd_original[__a][__b];
        RSDpd_reference[__a][__b] = RSDpd_original[__a][__b];
      }


//=======================================================
//
//             COMPUTE REFERENCE RESULTS
//
//=======================================================

    T               __mRSD[12] __attribute__ ((aligned (4)));
    T               __mRSDpd[12] __attribute__ ((aligned (4)));
    T               __mRSDpd_reference[12] __attribute__ ((aligned (4)));
    T               __mRSDpd_original[12] __attribute__ ((aligned (4)));
    for (int k = 0; k < 16; k++)
    {
      for (int __a = 0; __a < 12; __a++)
        __mRSD[__a] = RSD[__a][k];
      for (int __a = 0; __a < 12; __a++)
        __mRSDpd_reference[__a] = RSDpd_reference[__a][k];
      RSD_Positive_Definite_Part < float, float, int >(__mRSD,
                                                       __mRSDpd_reference);
      for (int __a = 0; __a < 12; __a++)
        RSDpd_reference[__a][k] = __mRSDpd_reference[__a];
    }

//=======================================================
//
//               COMPUTE SCALAR RESULTS
//
//=======================================================

    {
      typedef         T (&refArray1)[12][16];
      typedef         T (&refArray2)[12][16];
      for (int __a = 0; __a < 12; __a++)
        for (int __b = 0; __b < 16; __b++)
          RSDpd[__a][__b] = RSDpd_original[__a][__b];
      for (int i = 0; i < 16; i += 1)
      {
        refArray1       RSDk = reinterpret_cast < refArray1 > (RSD[0][i]);
        refArray2       RSDpdk = reinterpret_cast < refArray2 > (RSDpd[0][i]);
        RSD_Positive_Definite_Part < float, float[16], int[16] > (RSDk, RSDpdk);
      }
      for (int __a = 0; __a < 12; __a++)
        for (int __b = 0; __b < 16; __b++)
          if (std::
              abs ((RSDpd[__a][__b] -
                    RSDpd_reference[__a][__b]) / (RSDpd_reference[__a][__b])) >
              1)
          {
            std::cerr << "Mismatch detected in SCALAR implementation" << std::
              endl;
            std::cerr << "Variable RSDpd:" << std::endl;
            std::
              cerr << "seed=" << seed << ", __a=" << __a << ", __b=" << __b <<
              std::endl;
            std::cerr << "RSDpd SCALAR=  " << RSDpd[__a][__b] << std::endl;
            std::
              cerr << "RSDpd Reference=  " << RSDpd_reference[__a][__b] << std::
              endl;
            std::cerr << "RSDpd Rel Difference=  " << std::
              abs ((RSDpd[__a][__b] -
                    RSDpd_reference[__a][__b]) /
                   (RSDpd_reference[__a][__b])) << std::endl;
            std::cerr << "RSDpd Abs Difference=  " << std::
              abs (RSDpd[__a][__b] - RSDpd_reference[__a][__b]) << std::endl;
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
      typedef         T (&refArray2)[12][16];
      for (int __a = 0; __a < 12; __a++)
        for (int __b = 0; __b < 16; __b++)
          RSDpd[__a][__b] = RSDpd_original[__a][__b];
      for (int i = 0; i < 16; i += 4)
      {
        refArray1       RSDk = reinterpret_cast < refArray1 > (RSD[0][i]);
        refArray2       RSDpdk = reinterpret_cast < refArray2 > (RSDpd[0][i]);
        RSD_Positive_Definite_Part < __m128, float[16], int[16] > (RSDk,
                                                                   RSDpdk);
      }
      for (int __a = 0; __a < 12; __a++)
        for (int __b = 0; __b < 16; __b++)
          if (std::
              abs ((RSDpd[__a][__b] -
                    RSDpd_reference[__a][__b]) / (RSDpd_reference[__a][__b])) >
              1)
          {
            std::cerr << "Mismatch detected in SSE implementation" << std::endl;
            std::cerr << "Variable RSDpd:" << std::endl;
            std::
              cerr << "seed=" << seed << ", __a=" << __a << ", __b=" << __b <<
              std::endl;
            std::cerr << "RSDpd SSE=  " << RSDpd[__a][__b] << std::endl;
            std::
              cerr << "RSDpd Reference=  " << RSDpd_reference[__a][__b] << std::
              endl;
            std::cerr << "RSDpd Rel Difference=  " << std::
              abs ((RSDpd[__a][__b] -
                    RSDpd_reference[__a][__b]) /
                   (RSDpd_reference[__a][__b])) << std::endl;
            std::cerr << "RSDpd Abs Difference=  " << std::
              abs (RSDpd[__a][__b] - RSDpd_reference[__a][__b]) << std::endl;
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
      typedef         T (&refArray2)[12][16];
      for (int __a = 0; __a < 12; __a++)
        for (int __b = 0; __b < 16; __b++)
          RSDpd[__a][__b] = RSDpd_original[__a][__b];
      for (int i = 0; i < 16; i += 8)
      {
        refArray1       RSDk = reinterpret_cast < refArray1 > (RSD[0][i]);
        refArray2       RSDpdk = reinterpret_cast < refArray2 > (RSDpd[0][i]);
        RSD_Positive_Definite_Part < __m256, float[16], int[16] > (RSDk,
                                                                   RSDpdk);
      }
      for (int __a = 0; __a < 12; __a++)
        for (int __b = 0; __b < 16; __b++)
          if (std::
              abs ((RSDpd[__a][__b] -
                    RSDpd_reference[__a][__b]) / (RSDpd_reference[__a][__b])) >
              1)
          {
            std::cerr << "Mismatch detected in AVX implementation" << std::endl;
            std::cerr << "Variable RSDpd:" << std::endl;
            std::
              cerr << "seed=" << seed << ", __a=" << __a << ", __b=" << __b <<
              std::endl;
            std::cerr << "RSDpd AVX=  " << RSDpd[__a][__b] << std::endl;
            std::
              cerr << "RSDpd Reference=  " << RSDpd_reference[__a][__b] << std::
              endl;
            std::cerr << "RSDpd Rel Difference=  " << std::
              abs ((RSDpd[__a][__b] -
                    RSDpd_reference[__a][__b]) /
                   (RSDpd_reference[__a][__b])) << std::endl;
            std::cerr << "RSDpd Abs Difference=  " << std::
              abs (RSDpd[__a][__b] - RSDpd_reference[__a][__b]) << std::endl;
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
      typedef         T (&refArray2)[12][16];
      for (int __a = 0; __a < 12; __a++)
        for (int __b = 0; __b < 16; __b++)
          RSDpd[__a][__b] = RSDpd_original[__a][__b];
      for (int i = 0; i < 16; i += 4)
      {
        refArray1       RSDk = reinterpret_cast < refArray1 > (RSD[0][i]);
        refArray2       RSDpdk = reinterpret_cast < refArray2 > (RSDpd[0][i]);
        RSD_Positive_Definite_Part < float32x4_t, float[16], int[16] > (RSDk,
                                                                        RSDpdk);
      }
      for (int __a = 0; __a < 12; __a++)
        for (int __b = 0; __b < 16; __b++)
          if (std::
              abs ((RSDpd[__a][__b] -
                    RSDpd_reference[__a][__b]) / (RSDpd_reference[__a][__b])) >
              1)
          {
            std::cerr << "Mismatch detected in NEON implementation" << std::
              endl;
            std::cerr << "Variable RSDpd:" << std::endl;
            std::
              cerr << "seed=" << seed << ", __a=" << __a << ", __b=" << __b <<
              std::endl;
            std::cerr << "RSDpd NEON=  " << RSDpd[__a][__b] << std::endl;
            std::
              cerr << "RSDpd Reference=  " << RSDpd_reference[__a][__b] << std::
              endl;
            std::cerr << "RSDpd Rel Difference=  " << std::
              abs ((RSDpd[__a][__b] -
                    RSDpd_reference[__a][__b]) /
                   (RSDpd_reference[__a][__b])) << std::endl;
            std::cerr << "RSDpd Abs Difference=  " << std::
              abs (RSDpd[__a][__b] - RSDpd_reference[__a][__b]) << std::endl;
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
      typedef         T (&refArray2)[12][16];
      for (int __a = 0; __a < 12; __a++)
        for (int __b = 0; __b < 16; __b++)
          RSDpd[__a][__b] = RSDpd_original[__a][__b];
      for (int i = 0; i < 16; i += 16)
      {
        refArray1       RSDk = reinterpret_cast < refArray1 > (RSD[0][i]);
        refArray2       RSDpdk = reinterpret_cast < refArray2 > (RSDpd[0][i]);
        RSD_Positive_Definite_Part < __m512, float[16], int[16] > (RSDk,
                                                                   RSDpdk);
      }
      for (int __a = 0; __a < 12; __a++)
        for (int __b = 0; __b < 16; __b++)
          if (std::
              abs ((RSDpd[__a][__b] -
                    RSDpd_reference[__a][__b]) / (RSDpd_reference[__a][__b])) >
              1)
          {
            std::cerr << "Mismatch detected in MIC implementation" << std::endl;
            std::cerr << "Variable RSDpd:" << std::endl;
            std::
              cerr << "seed=" << seed << ", __a=" << __a << ", __b=" << __b <<
              std::endl;
            std::cerr << "RSDpd MIC=  " << RSDpd[__a][__b] << std::endl;
            std::
              cerr << "RSDpd Reference=  " << RSDpd_reference[__a][__b] << std::
              endl;
            std::cerr << "RSDpd Rel Difference=  " << std::
              abs ((RSDpd[__a][__b] -
                    RSDpd_reference[__a][__b]) /
                   (RSDpd_reference[__a][__b])) << std::endl;
            std::cerr << "RSDpd Abs Difference=  " << std::
              abs (RSDpd[__a][__b] - RSDpd_reference[__a][__b]) << std::endl;
            return 1;
          }

    }
#endif

  }



  std::cout << "SIMD check successful!" << std::endl;

  return 0;

}

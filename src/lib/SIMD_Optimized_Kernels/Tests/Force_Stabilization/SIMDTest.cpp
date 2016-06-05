
#include <cstdlib>
#include <iostream>
#include "KernelCommon.h"

struct NEOHOOKEAN_TAG;
struct COROTATED_TAG;

#include "Force_Stabilization.h"
#include "Force_Stabilization_Reference.h"

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
    std::cout << "Running SIMD Test for Force_Stabilization " << std::endl;


//=======================================================
//
//        DEFINE ALL VARIABLES USED BY KERNEL
//
//=======================================================

    T               Du[3][8][16] __attribute__ ((aligned (64)));
    T               constant[16] __attribute__ ((aligned (64)));
    T               dH[3][8][16] __attribute__ ((aligned (64)));
    T               dH_reference[3][8][16] __attribute__ ((aligned (64)));
    T               dH_original[3][8][16] __attribute__ ((aligned (64)));


    for (int __a = 0; __a < 3; __a++)
      for (int __b = 0; __b < 8; __b++)
        for (int __c = 0; __c < 16; __c++)
          Du[__a][__b][__c] = Get_Random < float >();
    for (int __a = 0; __a < 16; __a++)
      constant[__a] = Get_Random < float >();
    for (int __a = 0; __a < 3; __a++)
      for (int __b = 0; __b < 8; __b++)
        for (int __c = 0; __c < 16; __c++)
        {
          dH_original[__a][__b][__c] = Get_Random < float >();
          dH[__a][__b][__c] = dH_original[__a][__b][__c];
          dH_reference[__a][__b][__c] = dH_original[__a][__b][__c];
        }


//=======================================================
//
//             COMPUTE REFERENCE RESULTS
//
//=======================================================

    T               __mDu[3][8] __attribute__ ((aligned (4)));
    T __mconstant __attribute__ ((aligned (4)));
    T               __mdH[3][8] __attribute__ ((aligned (4)));
    T               __mdH_reference[3][8] __attribute__ ((aligned (4)));
    T               __mdH_original[3][8] __attribute__ ((aligned (4)));
    for (int k = 0; k < 16; k++)
    {
      for (int __a = 0; __a < 3; __a++)
        for (int __b = 0; __b < 8; __b++)
          __mDu[__a][__b] = Du[__a][__b][k];
      __mconstant = constant[k];
      for (int __a = 0; __a < 3; __a++)
        for (int __b = 0; __b < 8; __b++)
          __mdH_reference[__a][__b] = dH_reference[__a][__b][k];
      Force_Stabilization < float, float, int >(__mDu, __mconstant,
                                                __mdH_reference);
      for (int __a = 0; __a < 3; __a++)
        for (int __b = 0; __b < 8; __b++)
          dH_reference[__a][__b][k] = __mdH_reference[__a][__b];
    }

//=======================================================
//
//               COMPUTE SCALAR RESULTS
//
//=======================================================

    {
      typedef         T (&refArray1)[3][8][16];
      typedef         T (&refArray2)[16];
      typedef         T (&refArray3)[3][8][16];
      for (int __a = 0; __a < 3; __a++)
        for (int __b = 0; __b < 8; __b++)
          for (int __c = 0; __c < 16; __c++)
            dH[__a][__b][__c] = dH_original[__a][__b][__c];
      for (int i = 0; i < 16; i += 1)
      {
        refArray1       Duk = reinterpret_cast < refArray1 > (Du[0][0][i]);
        refArray2       constantk =
          reinterpret_cast < refArray2 > (constant[i]);
        refArray3       dHk = reinterpret_cast < refArray3 > (dH[0][0][i]);
        Force_Stabilization < float, float[16], int[16] > (Duk, constantk, dHk);
      }
      for (int __a = 0; __a < 3; __a++)
        for (int __b = 0; __b < 8; __b++)
          for (int __c = 0; __c < 16; __c++)
            if (std::
                abs ((dH[__a][__b][__c] -
                      dH_reference[__a][__b][__c]) /
                     (dH_reference[__a][__b][__c])) > 1)
            {
              std::cerr << "Mismatch detected in SCALAR implementation" << std::
                endl;
              std::cerr << "Variable dH:" << std::endl;
              std::
                cerr << "seed=" << seed << ", __a=" << __a << ", __b=" << __b <<
                ", __c=" << __c << std::endl;
              std::cerr << "dH SCALAR=  " << dH[__a][__b][__c] << std::endl;
              std::
                cerr << "dH Reference=  " << dH_reference[__a][__b][__c] <<
                std::endl;
              std::cerr << "dH Rel Difference=  " << std::
                abs ((dH[__a][__b][__c] -
                      dH_reference[__a][__b][__c]) /
                     (dH_reference[__a][__b][__c])) << std::endl;
              std::cerr << "dH Abs Difference=  " << std::
                abs (dH[__a][__b][__c] -
                     dH_reference[__a][__b][__c]) << std::endl;
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
      typedef         T (&refArray2)[16];
      typedef         T (&refArray3)[3][8][16];
      for (int __a = 0; __a < 3; __a++)
        for (int __b = 0; __b < 8; __b++)
          for (int __c = 0; __c < 16; __c++)
            dH[__a][__b][__c] = dH_original[__a][__b][__c];
      for (int i = 0; i < 16; i += 4)
      {
        refArray1       Duk = reinterpret_cast < refArray1 > (Du[0][0][i]);
        refArray2       constantk =
          reinterpret_cast < refArray2 > (constant[i]);
        refArray3       dHk = reinterpret_cast < refArray3 > (dH[0][0][i]);
        Force_Stabilization < __m128, float[16], int[16] > (Duk, constantk,
                                                            dHk);
      }
      for (int __a = 0; __a < 3; __a++)
        for (int __b = 0; __b < 8; __b++)
          for (int __c = 0; __c < 16; __c++)
            if (std::
                abs ((dH[__a][__b][__c] -
                      dH_reference[__a][__b][__c]) /
                     (dH_reference[__a][__b][__c])) > 1)
            {
              std::cerr << "Mismatch detected in SSE implementation" << std::
                endl;
              std::cerr << "Variable dH:" << std::endl;
              std::
                cerr << "seed=" << seed << ", __a=" << __a << ", __b=" << __b <<
                ", __c=" << __c << std::endl;
              std::cerr << "dH SSE=  " << dH[__a][__b][__c] << std::endl;
              std::
                cerr << "dH Reference=  " << dH_reference[__a][__b][__c] <<
                std::endl;
              std::cerr << "dH Rel Difference=  " << std::
                abs ((dH[__a][__b][__c] -
                      dH_reference[__a][__b][__c]) /
                     (dH_reference[__a][__b][__c])) << std::endl;
              std::cerr << "dH Abs Difference=  " << std::
                abs (dH[__a][__b][__c] -
                     dH_reference[__a][__b][__c]) << std::endl;
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
      typedef         T (&refArray2)[16];
      typedef         T (&refArray3)[3][8][16];
      for (int __a = 0; __a < 3; __a++)
        for (int __b = 0; __b < 8; __b++)
          for (int __c = 0; __c < 16; __c++)
            dH[__a][__b][__c] = dH_original[__a][__b][__c];
      for (int i = 0; i < 16; i += 8)
      {
        refArray1       Duk = reinterpret_cast < refArray1 > (Du[0][0][i]);
        refArray2       constantk =
          reinterpret_cast < refArray2 > (constant[i]);
        refArray3       dHk = reinterpret_cast < refArray3 > (dH[0][0][i]);
        Force_Stabilization < __m256, float[16], int[16] > (Duk, constantk,
                                                            dHk);
      }
      for (int __a = 0; __a < 3; __a++)
        for (int __b = 0; __b < 8; __b++)
          for (int __c = 0; __c < 16; __c++)
            if (std::
                abs ((dH[__a][__b][__c] -
                      dH_reference[__a][__b][__c]) /
                     (dH_reference[__a][__b][__c])) > 1)
            {
              std::cerr << "Mismatch detected in AVX implementation" << std::
                endl;
              std::cerr << "Variable dH:" << std::endl;
              std::
                cerr << "seed=" << seed << ", __a=" << __a << ", __b=" << __b <<
                ", __c=" << __c << std::endl;
              std::cerr << "dH AVX=  " << dH[__a][__b][__c] << std::endl;
              std::
                cerr << "dH Reference=  " << dH_reference[__a][__b][__c] <<
                std::endl;
              std::cerr << "dH Rel Difference=  " << std::
                abs ((dH[__a][__b][__c] -
                      dH_reference[__a][__b][__c]) /
                     (dH_reference[__a][__b][__c])) << std::endl;
              std::cerr << "dH Abs Difference=  " << std::
                abs (dH[__a][__b][__c] -
                     dH_reference[__a][__b][__c]) << std::endl;
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
      typedef         T (&refArray2)[16];
      typedef         T (&refArray3)[3][8][16];
      for (int __a = 0; __a < 3; __a++)
        for (int __b = 0; __b < 8; __b++)
          for (int __c = 0; __c < 16; __c++)
            dH[__a][__b][__c] = dH_original[__a][__b][__c];
      for (int i = 0; i < 16; i += 4)
      {
        refArray1       Duk = reinterpret_cast < refArray1 > (Du[0][0][i]);
        refArray2       constantk =
          reinterpret_cast < refArray2 > (constant[i]);
        refArray3       dHk = reinterpret_cast < refArray3 > (dH[0][0][i]);
        Force_Stabilization < float32x4_t, float[16], int[16] > (Duk, constantk,
                                                                 dHk);
      }
      for (int __a = 0; __a < 3; __a++)
        for (int __b = 0; __b < 8; __b++)
          for (int __c = 0; __c < 16; __c++)
            if (std::
                abs ((dH[__a][__b][__c] -
                      dH_reference[__a][__b][__c]) /
                     (dH_reference[__a][__b][__c])) > 1)
            {
              std::cerr << "Mismatch detected in NEON implementation" << std::
                endl;
              std::cerr << "Variable dH:" << std::endl;
              std::
                cerr << "seed=" << seed << ", __a=" << __a << ", __b=" << __b <<
                ", __c=" << __c << std::endl;
              std::cerr << "dH NEON=  " << dH[__a][__b][__c] << std::endl;
              std::
                cerr << "dH Reference=  " << dH_reference[__a][__b][__c] <<
                std::endl;
              std::cerr << "dH Rel Difference=  " << std::
                abs ((dH[__a][__b][__c] -
                      dH_reference[__a][__b][__c]) /
                     (dH_reference[__a][__b][__c])) << std::endl;
              std::cerr << "dH Abs Difference=  " << std::
                abs (dH[__a][__b][__c] -
                     dH_reference[__a][__b][__c]) << std::endl;
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
      typedef         T (&refArray2)[16];
      typedef         T (&refArray3)[3][8][16];
      for (int __a = 0; __a < 3; __a++)
        for (int __b = 0; __b < 8; __b++)
          for (int __c = 0; __c < 16; __c++)
            dH[__a][__b][__c] = dH_original[__a][__b][__c];
      for (int i = 0; i < 16; i += 16)
      {
        refArray1       Duk = reinterpret_cast < refArray1 > (Du[0][0][i]);
        refArray2       constantk =
          reinterpret_cast < refArray2 > (constant[i]);
        refArray3       dHk = reinterpret_cast < refArray3 > (dH[0][0][i]);
        Force_Stabilization < __m512, float[16], int[16] > (Duk, constantk,
                                                            dHk);
      }
      for (int __a = 0; __a < 3; __a++)
        for (int __b = 0; __b < 8; __b++)
          for (int __c = 0; __c < 16; __c++)
            if (std::
                abs ((dH[__a][__b][__c] -
                      dH_reference[__a][__b][__c]) /
                     (dH_reference[__a][__b][__c])) > 1)
            {
              std::cerr << "Mismatch detected in MIC implementation" << std::
                endl;
              std::cerr << "Variable dH:" << std::endl;
              std::
                cerr << "seed=" << seed << ", __a=" << __a << ", __b=" << __b <<
                ", __c=" << __c << std::endl;
              std::cerr << "dH MIC=  " << dH[__a][__b][__c] << std::endl;
              std::
                cerr << "dH Reference=  " << dH_reference[__a][__b][__c] <<
                std::endl;
              std::cerr << "dH Rel Difference=  " << std::
                abs ((dH[__a][__b][__c] -
                      dH_reference[__a][__b][__c]) /
                     (dH_reference[__a][__b][__c])) << std::endl;
              std::cerr << "dH Abs Difference=  " << std::
                abs (dH[__a][__b][__c] -
                     dH_reference[__a][__b][__c]) << std::endl;
              return 1;
            }

    }
#endif

  }



  std::cout << "SIMD check successful!" << std::endl;

  return 0;

}

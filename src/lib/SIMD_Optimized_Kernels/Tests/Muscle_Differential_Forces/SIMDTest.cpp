
#include <cstdlib>
#include <iostream>
#include "KernelCommon.h"

struct NEOHOOKEAN_TAG;
struct COROTATED_TAG;

#include "Muscle_Differential_Forces.h"
#include "Muscle_Differential_Forces_Reference.h"

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
    std::cout << "Running SIMD Test for Muscle_Differential_Forces " << std::
      endl;


//=======================================================
//
//        DEFINE ALL VARIABLES USED BY KERNEL
//
//=======================================================

    T               df[3][8][16] __attribute__ ((aligned (64)));
    T               df_reference[3][8][16] __attribute__ ((aligned (64)));
    T               df_original[3][8][16] __attribute__ ((aligned (64)));
    T               du[3][8][16] __attribute__ ((aligned (64)));
    T               fiber[3][16] __attribute__ ((aligned (64)));
    T               Ffiber[3][16] __attribute__ ((aligned (64)));
    T               c1[16] __attribute__ ((aligned (64)));
    T               c2[16] __attribute__ ((aligned (64)));
    T               one_over_h[16] __attribute__ ((aligned (64)));
    T               cell_volume[16] __attribute__ ((aligned (64)));


    for (int __a = 0; __a < 3; __a++)
      for (int __b = 0; __b < 8; __b++)
        for (int __c = 0; __c < 16; __c++)
        {
          df_original[__a][__b][__c] = Get_Random < float >();
          df[__a][__b][__c] = df_original[__a][__b][__c];
          df_reference[__a][__b][__c] = df_original[__a][__b][__c];
        }
    for (int __a = 0; __a < 3; __a++)
      for (int __b = 0; __b < 8; __b++)
        for (int __c = 0; __c < 16; __c++)
          du[__a][__b][__c] = Get_Random < float >();
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
    for (int __a = 0; __a < 16; __a++)
      one_over_h[__a] = Get_Random < float >();
    for (int __a = 0; __a < 16; __a++)
      cell_volume[__a] = Get_Random < float >();

//=======================================================
//
//             COMPUTE REFERENCE RESULTS
//
//=======================================================

    T               __mdf[3][8] __attribute__ ((aligned (4)));
    T               __mdf_reference[3][8] __attribute__ ((aligned (4)));
    T               __mdf_original[3][8] __attribute__ ((aligned (4)));
    T               __mdu[3][8] __attribute__ ((aligned (4)));
    T               __mfiber[3] __attribute__ ((aligned (4)));
    T               __mFfiber[3] __attribute__ ((aligned (4)));
    T __mc1 __attribute__ ((aligned (4)));
    T __mc2 __attribute__ ((aligned (4)));
    T __mone_over_h __attribute__ ((aligned (4)));
    T __mcell_volume __attribute__ ((aligned (4)));
    for (int k = 0; k < 16; k++)
    {
      for (int __a = 0; __a < 3; __a++)
        for (int __b = 0; __b < 8; __b++)
          __mdf_reference[__a][__b] = df_reference[__a][__b][k];
      for (int __a = 0; __a < 3; __a++)
        for (int __b = 0; __b < 8; __b++)
          __mdu[__a][__b] = du[__a][__b][k];
      for (int __a = 0; __a < 3; __a++)
        __mfiber[__a] = fiber[__a][k];
      for (int __a = 0; __a < 3; __a++)
        __mFfiber[__a] = Ffiber[__a][k];
      __mc1 = c1[k];
      __mc2 = c2[k];
      __mone_over_h = one_over_h[k];
      __mcell_volume = cell_volume[k];
      Muscle_Differential_Forces < float, float, int >(__mdf_reference, __mdu,
                                                       __mfiber, __mFfiber,
                                                       __mc1, __mc2,
                                                       __mone_over_h,
                                                       __mcell_volume);
      for (int __a = 0; __a < 3; __a++)
        for (int __b = 0; __b < 8; __b++)
          df_reference[__a][__b][k] = __mdf_reference[__a][__b];
    }

//=======================================================
//
//               COMPUTE SCALAR RESULTS
//
//=======================================================

    {
      typedef         T (&refArray1)[3][8][16];
      typedef         T (&refArray2)[3][8][16];
      typedef         T (&refArray3)[3][16];
      typedef         T (&refArray4)[3][16];
      typedef         T (&refArray5)[16];
      typedef         T (&refArray6)[16];
      typedef         T (&refArray7)[16];
      typedef         T (&refArray8)[16];
      for (int __a = 0; __a < 3; __a++)
        for (int __b = 0; __b < 8; __b++)
          for (int __c = 0; __c < 16; __c++)
            df[__a][__b][__c] = df_original[__a][__b][__c];
      for (int i = 0; i < 16; i += 1)
      {
        refArray1       dfk = reinterpret_cast < refArray1 > (df[0][0][i]);
        refArray2       duk = reinterpret_cast < refArray2 > (du[0][0][i]);
        refArray3       fiberk = reinterpret_cast < refArray3 > (fiber[0][i]);
        refArray4       Ffiberk = reinterpret_cast < refArray4 > (Ffiber[0][i]);
        refArray5       c1k = reinterpret_cast < refArray5 > (c1[i]);
        refArray6       c2k = reinterpret_cast < refArray6 > (c2[i]);
        refArray7       one_over_hk =
          reinterpret_cast < refArray7 > (one_over_h[i]);
        refArray8       cell_volumek =
          reinterpret_cast < refArray8 > (cell_volume[i]);
        Muscle_Differential_Forces < float, float[16], int[16] > (dfk, duk,
                                                                  fiberk,
                                                                  Ffiberk, c1k,
                                                                  c2k,
                                                                  one_over_hk,
                                                                  cell_volumek);
      }
      for (int __a = 0; __a < 3; __a++)
        for (int __b = 0; __b < 8; __b++)
          for (int __c = 0; __c < 16; __c++)
            if (std::
                abs ((df[__a][__b][__c] -
                      df_reference[__a][__b][__c]) /
                     (df_reference[__a][__b][__c])) > 1)
            {
              std::cerr << "Mismatch detected in SCALAR implementation" << std::
                endl;
              std::cerr << "Variable df:" << std::endl;
              std::
                cerr << "seed=" << seed << ", __a=" << __a << ", __b=" << __b <<
                ", __c=" << __c << std::endl;
              std::cerr << "df SCALAR=  " << df[__a][__b][__c] << std::endl;
              std::
                cerr << "df Reference=  " << df_reference[__a][__b][__c] <<
                std::endl;
              std::cerr << "df Rel Difference=  " << std::
                abs ((df[__a][__b][__c] -
                      df_reference[__a][__b][__c]) /
                     (df_reference[__a][__b][__c])) << std::endl;
              std::cerr << "df Abs Difference=  " << std::
                abs (df[__a][__b][__c] -
                     df_reference[__a][__b][__c]) << std::endl;
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
      typedef         T (&refArray2)[3][8][16];
      typedef         T (&refArray3)[3][16];
      typedef         T (&refArray4)[3][16];
      typedef         T (&refArray5)[16];
      typedef         T (&refArray6)[16];
      typedef         T (&refArray7)[16];
      typedef         T (&refArray8)[16];
      for (int __a = 0; __a < 3; __a++)
        for (int __b = 0; __b < 8; __b++)
          for (int __c = 0; __c < 16; __c++)
            df[__a][__b][__c] = df_original[__a][__b][__c];
      for (int i = 0; i < 16; i += 4)
      {
        refArray1       dfk = reinterpret_cast < refArray1 > (df[0][0][i]);
        refArray2       duk = reinterpret_cast < refArray2 > (du[0][0][i]);
        refArray3       fiberk = reinterpret_cast < refArray3 > (fiber[0][i]);
        refArray4       Ffiberk = reinterpret_cast < refArray4 > (Ffiber[0][i]);
        refArray5       c1k = reinterpret_cast < refArray5 > (c1[i]);
        refArray6       c2k = reinterpret_cast < refArray6 > (c2[i]);
        refArray7       one_over_hk =
          reinterpret_cast < refArray7 > (one_over_h[i]);
        refArray8       cell_volumek =
          reinterpret_cast < refArray8 > (cell_volume[i]);
        Muscle_Differential_Forces < __m128, float[16], int[16] > (dfk, duk,
                                                                   fiberk,
                                                                   Ffiberk, c1k,
                                                                   c2k,
                                                                   one_over_hk,
                                                                   cell_volumek);
      }
      for (int __a = 0; __a < 3; __a++)
        for (int __b = 0; __b < 8; __b++)
          for (int __c = 0; __c < 16; __c++)
            if (std::
                abs ((df[__a][__b][__c] -
                      df_reference[__a][__b][__c]) /
                     (df_reference[__a][__b][__c])) > 1)
            {
              std::cerr << "Mismatch detected in SSE implementation" << std::
                endl;
              std::cerr << "Variable df:" << std::endl;
              std::
                cerr << "seed=" << seed << ", __a=" << __a << ", __b=" << __b <<
                ", __c=" << __c << std::endl;
              std::cerr << "df SSE=  " << df[__a][__b][__c] << std::endl;
              std::
                cerr << "df Reference=  " << df_reference[__a][__b][__c] <<
                std::endl;
              std::cerr << "df Rel Difference=  " << std::
                abs ((df[__a][__b][__c] -
                      df_reference[__a][__b][__c]) /
                     (df_reference[__a][__b][__c])) << std::endl;
              std::cerr << "df Abs Difference=  " << std::
                abs (df[__a][__b][__c] -
                     df_reference[__a][__b][__c]) << std::endl;
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
      typedef         T (&refArray2)[3][8][16];
      typedef         T (&refArray3)[3][16];
      typedef         T (&refArray4)[3][16];
      typedef         T (&refArray5)[16];
      typedef         T (&refArray6)[16];
      typedef         T (&refArray7)[16];
      typedef         T (&refArray8)[16];
      for (int __a = 0; __a < 3; __a++)
        for (int __b = 0; __b < 8; __b++)
          for (int __c = 0; __c < 16; __c++)
            df[__a][__b][__c] = df_original[__a][__b][__c];
      for (int i = 0; i < 16; i += 8)
      {
        refArray1       dfk = reinterpret_cast < refArray1 > (df[0][0][i]);
        refArray2       duk = reinterpret_cast < refArray2 > (du[0][0][i]);
        refArray3       fiberk = reinterpret_cast < refArray3 > (fiber[0][i]);
        refArray4       Ffiberk = reinterpret_cast < refArray4 > (Ffiber[0][i]);
        refArray5       c1k = reinterpret_cast < refArray5 > (c1[i]);
        refArray6       c2k = reinterpret_cast < refArray6 > (c2[i]);
        refArray7       one_over_hk =
          reinterpret_cast < refArray7 > (one_over_h[i]);
        refArray8       cell_volumek =
          reinterpret_cast < refArray8 > (cell_volume[i]);
        Muscle_Differential_Forces < __m256, float[16], int[16] > (dfk, duk,
                                                                   fiberk,
                                                                   Ffiberk, c1k,
                                                                   c2k,
                                                                   one_over_hk,
                                                                   cell_volumek);
      }
      for (int __a = 0; __a < 3; __a++)
        for (int __b = 0; __b < 8; __b++)
          for (int __c = 0; __c < 16; __c++)
            if (std::
                abs ((df[__a][__b][__c] -
                      df_reference[__a][__b][__c]) /
                     (df_reference[__a][__b][__c])) > 1)
            {
              std::cerr << "Mismatch detected in AVX implementation" << std::
                endl;
              std::cerr << "Variable df:" << std::endl;
              std::
                cerr << "seed=" << seed << ", __a=" << __a << ", __b=" << __b <<
                ", __c=" << __c << std::endl;
              std::cerr << "df AVX=  " << df[__a][__b][__c] << std::endl;
              std::
                cerr << "df Reference=  " << df_reference[__a][__b][__c] <<
                std::endl;
              std::cerr << "df Rel Difference=  " << std::
                abs ((df[__a][__b][__c] -
                      df_reference[__a][__b][__c]) /
                     (df_reference[__a][__b][__c])) << std::endl;
              std::cerr << "df Abs Difference=  " << std::
                abs (df[__a][__b][__c] -
                     df_reference[__a][__b][__c]) << std::endl;
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
      typedef         T (&refArray2)[3][8][16];
      typedef         T (&refArray3)[3][16];
      typedef         T (&refArray4)[3][16];
      typedef         T (&refArray5)[16];
      typedef         T (&refArray6)[16];
      typedef         T (&refArray7)[16];
      typedef         T (&refArray8)[16];
      for (int __a = 0; __a < 3; __a++)
        for (int __b = 0; __b < 8; __b++)
          for (int __c = 0; __c < 16; __c++)
            df[__a][__b][__c] = df_original[__a][__b][__c];
      for (int i = 0; i < 16; i += 4)
      {
        refArray1       dfk = reinterpret_cast < refArray1 > (df[0][0][i]);
        refArray2       duk = reinterpret_cast < refArray2 > (du[0][0][i]);
        refArray3       fiberk = reinterpret_cast < refArray3 > (fiber[0][i]);
        refArray4       Ffiberk = reinterpret_cast < refArray4 > (Ffiber[0][i]);
        refArray5       c1k = reinterpret_cast < refArray5 > (c1[i]);
        refArray6       c2k = reinterpret_cast < refArray6 > (c2[i]);
        refArray7       one_over_hk =
          reinterpret_cast < refArray7 > (one_over_h[i]);
        refArray8       cell_volumek =
          reinterpret_cast < refArray8 > (cell_volume[i]);
        Muscle_Differential_Forces < float32x4_t, float[16], int[16] > (dfk,
                                                                        duk,
                                                                        fiberk,
                                                                        Ffiberk,
                                                                        c1k,
                                                                        c2k,
                                                                        one_over_hk,
                                                                        cell_volumek);
      }
      for (int __a = 0; __a < 3; __a++)
        for (int __b = 0; __b < 8; __b++)
          for (int __c = 0; __c < 16; __c++)
            if (std::
                abs ((df[__a][__b][__c] -
                      df_reference[__a][__b][__c]) /
                     (df_reference[__a][__b][__c])) > 1)
            {
              std::cerr << "Mismatch detected in NEON implementation" << std::
                endl;
              std::cerr << "Variable df:" << std::endl;
              std::
                cerr << "seed=" << seed << ", __a=" << __a << ", __b=" << __b <<
                ", __c=" << __c << std::endl;
              std::cerr << "df NEON=  " << df[__a][__b][__c] << std::endl;
              std::
                cerr << "df Reference=  " << df_reference[__a][__b][__c] <<
                std::endl;
              std::cerr << "df Rel Difference=  " << std::
                abs ((df[__a][__b][__c] -
                      df_reference[__a][__b][__c]) /
                     (df_reference[__a][__b][__c])) << std::endl;
              std::cerr << "df Abs Difference=  " << std::
                abs (df[__a][__b][__c] -
                     df_reference[__a][__b][__c]) << std::endl;
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
      typedef         T (&refArray2)[3][8][16];
      typedef         T (&refArray3)[3][16];
      typedef         T (&refArray4)[3][16];
      typedef         T (&refArray5)[16];
      typedef         T (&refArray6)[16];
      typedef         T (&refArray7)[16];
      typedef         T (&refArray8)[16];
      for (int __a = 0; __a < 3; __a++)
        for (int __b = 0; __b < 8; __b++)
          for (int __c = 0; __c < 16; __c++)
            df[__a][__b][__c] = df_original[__a][__b][__c];
      for (int i = 0; i < 16; i += 16)
      {
        refArray1       dfk = reinterpret_cast < refArray1 > (df[0][0][i]);
        refArray2       duk = reinterpret_cast < refArray2 > (du[0][0][i]);
        refArray3       fiberk = reinterpret_cast < refArray3 > (fiber[0][i]);
        refArray4       Ffiberk = reinterpret_cast < refArray4 > (Ffiber[0][i]);
        refArray5       c1k = reinterpret_cast < refArray5 > (c1[i]);
        refArray6       c2k = reinterpret_cast < refArray6 > (c2[i]);
        refArray7       one_over_hk =
          reinterpret_cast < refArray7 > (one_over_h[i]);
        refArray8       cell_volumek =
          reinterpret_cast < refArray8 > (cell_volume[i]);
        Muscle_Differential_Forces < __m512, float[16], int[16] > (dfk, duk,
                                                                   fiberk,
                                                                   Ffiberk, c1k,
                                                                   c2k,
                                                                   one_over_hk,
                                                                   cell_volumek);
      }
      for (int __a = 0; __a < 3; __a++)
        for (int __b = 0; __b < 8; __b++)
          for (int __c = 0; __c < 16; __c++)
            if (std::
                abs ((df[__a][__b][__c] -
                      df_reference[__a][__b][__c]) /
                     (df_reference[__a][__b][__c])) > 1)
            {
              std::cerr << "Mismatch detected in MIC implementation" << std::
                endl;
              std::cerr << "Variable df:" << std::endl;
              std::
                cerr << "seed=" << seed << ", __a=" << __a << ", __b=" << __b <<
                ", __c=" << __c << std::endl;
              std::cerr << "df MIC=  " << df[__a][__b][__c] << std::endl;
              std::
                cerr << "df Reference=  " << df_reference[__a][__b][__c] <<
                std::endl;
              std::cerr << "df Rel Difference=  " << std::
                abs ((df[__a][__b][__c] -
                      df_reference[__a][__b][__c]) /
                     (df_reference[__a][__b][__c])) << std::endl;
              std::cerr << "df Abs Difference=  " << std::
                abs (df[__a][__b][__c] -
                     df_reference[__a][__b][__c]) << std::endl;
              return 1;
            }

    }
#endif

  }



  std::cout << "SIMD check successful!" << std::endl;

  return 0;

}

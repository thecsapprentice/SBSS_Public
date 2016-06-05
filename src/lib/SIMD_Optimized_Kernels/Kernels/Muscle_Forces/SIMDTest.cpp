
struct NEOHOOKEAN_TAG;
struct COROTATED_TAG;

#include <cstdlib>
#include <iostream>
#include <immintrin.h>

#include "Muscle_Forces.h"
#include "Muscle_Forces_Reference.h"

template < class T > T Get_Random (const T a = (T) - 1., const T b = (T) 1.)
{
  return ((b - a) * (T) rand ()) / (T) RAND_MAX + a;
}

int
main (int argc, char *argv[])
{
  typedef float T;
  typedef NEOHOOKEAN_TAG T_MATERIAL_TAG;

  int seed = 1;
  if (argc == 2)
    seed = atoi (argv[1]);
  srand (seed);

  std::cout.precision (10);
  std::cout.setf (std::ios::fixed, std::ios::floatfield);



  {
    T f[3][8][8];
    T f_reference[3][8][8];
    T f_original[3][8][8];
    T fiber[3][8];
    T Ffiber[3][8];
    T c1[8];
    T one_over_h[8];
    T cell_volume[8];

    for (int __a = 0; __a < 3; __a++)
      for (int __b = 0; __b < 8; __b++)
        for (int __c = 0; __c < 8; __c++)
          {
            f_original[__a][__b][__c] = Get_Random < float >();
            f[__a][__b][__c] = f_original[__a][__b][__c];
            f_reference[__a][__b][__c] = f_original[__a][__b][__c];
          }
    for (int __a = 0; __a < 3; __a++)
      for (int __b = 0; __b < 8; __b++)
        fiber[__a][__b] = Get_Random < float >();
    for (int __a = 0; __a < 3; __a++)
      for (int __b = 0; __b < 8; __b++)
        Ffiber[__a][__b] = Get_Random < float >();
    for (int __a = 0; __a < 8; __a++)
      c1[__a] = Get_Random < float >();
    for (int __a = 0; __a < 8; __a++)
      one_over_h[__a] = Get_Random < float >();
    for (int __a = 0; __a < 8; __a++)
      cell_volume[__a] = Get_Random < float >();

    typedef T (&refArray1)[3][8][8];
    typedef T (&refArray2)[3][8];
    typedef T (&refArray3)[3][8];
    typedef T (&refArray4)[8];
    typedef T (&refArray5)[8];
    typedef T (&refArray6)[8];


    for (int k = 0; k < 8; k++)
      {
        T __mf_reference[3][8];
        T __mfiber[3];
        T __mFfiber[3];
        T __mc1;
        T __mone_over_h;
        T __mcell_volume;
        for (int __a = 0; __a < 3; __a++)
          for (int __b = 0; __b < 8; __b++)
            __mf_reference[__a][__b] = f_original[__a][__b][k];
        for (int __a = 0; __a < 3; __a++)
          __mfiber[__a] = fiber[__a][k];
        for (int __a = 0; __a < 3; __a++)
          __mFfiber[__a] = Ffiber[__a][k];
        __mc1 = c1[k];
        __mone_over_h = one_over_h[k];
        __mcell_volume = cell_volume[k];
        Muscle_Forces < T, T, 1 > (__mf_reference, __mfiber, __mFfiber, __mc1,
                                   __mone_over_h, __mcell_volume);
        for (int __a = 0; __a < 3; __a++)
          for (int __b = 0; __b < 8; __b++)
            f_reference[__a][__b][k] = __mf_reference[__a][__b];
      }

    {
      for (int __a = 0; __a < 3; __a++)
        for (int __b = 0; __b < 8; __b++)
          for (int __c = 0; __c < 8; __c++)
            f[__a][__b][__c] = f_original[__a][__b][__c];
      for (int i = 0; i < 8; i += 1)
        {
          refArray1 fk = reinterpret_cast < refArray1 > (f[0][0][i]);
          refArray2 fiberk = reinterpret_cast < refArray2 > (fiber[0][i]);
          refArray3 Ffiberk = reinterpret_cast < refArray3 > (Ffiber[0][i]);
          refArray4 c1k = reinterpret_cast < refArray4 > (c1[i]);
          refArray5 one_over_hk =
            reinterpret_cast < refArray5 > (one_over_h[i]);
          refArray6 cell_volumek =
            reinterpret_cast < refArray6 > (cell_volume[i]);
          Muscle_Forces < float, float[8], 1 > (fk, fiberk, Ffiberk, c1k,
                                                one_over_hk, cell_volumek);
        }
      for (int __a = 0; __a < 3; __a++)
        for (int __b = 0; __b < 8; __b++)
          for (int __c = 0; __c < 8; __c++)
            if (f[__a][__b][__c] != f_reference[__a][__b][__c])
              {
                std::
                  cerr << "Mismatch detected in SCALAR implementation" << std::
                  endl;
                std::cerr << "Variable f:" << std::endl;
                std::
                  cerr << "seed=" << seed << ", __a=" << __a << ", __b=" << __b
                  << ", __c=" << __c << std::endl;
                std::cerr << "f SCALAR=  " << f[__a][__b][__c] << std::endl;
                std::
                  cerr << "f Reference=  " << f_reference[__a][__b][__c] <<
                  std::endl;
                return 1;
              }

    }

    {
      for (int __a = 0; __a < 3; __a++)
        for (int __b = 0; __b < 8; __b++)
          for (int __c = 0; __c < 8; __c++)
            f[__a][__b][__c] = f_original[__a][__b][__c];
      for (int i = 0; i < 8; i += 4)
        {
          refArray1 fk = reinterpret_cast < refArray1 > (f[0][0][i]);
          refArray2 fiberk = reinterpret_cast < refArray2 > (fiber[0][i]);
          refArray3 Ffiberk = reinterpret_cast < refArray3 > (Ffiber[0][i]);
          refArray4 c1k = reinterpret_cast < refArray4 > (c1[i]);
          refArray5 one_over_hk =
            reinterpret_cast < refArray5 > (one_over_h[i]);
          refArray6 cell_volumek =
            reinterpret_cast < refArray6 > (cell_volume[i]);
          Muscle_Forces < __m128, float[8], 4 > (fk, fiberk, Ffiberk, c1k,
                                                 one_over_hk, cell_volumek);
        }
      for (int __a = 0; __a < 3; __a++)
        for (int __b = 0; __b < 8; __b++)
          for (int __c = 0; __c < 8; __c++)
            if (f[__a][__b][__c] != f_reference[__a][__b][__c])
              {
                std::cerr << "Mismatch detected in SSE implementation" << std::
                  endl;
                std::cerr << "Variable f:" << std::endl;
                std::
                  cerr << "seed=" << seed << ", __a=" << __a << ", __b=" << __b
                  << ", __c=" << __c << std::endl;
                std::cerr << "f SSE=  " << f[__a][__b][__c] << std::endl;
                std::
                  cerr << "f Reference=  " << f_reference[__a][__b][__c] <<
                  std::endl;
                return 1;
              }

    }

    {
      for (int __a = 0; __a < 3; __a++)
        for (int __b = 0; __b < 8; __b++)
          for (int __c = 0; __c < 8; __c++)
            f[__a][__b][__c] = f_original[__a][__b][__c];
      for (int i = 0; i < 8; i += 8)
        {
          refArray1 fk = reinterpret_cast < refArray1 > (f[0][0][i]);
          refArray2 fiberk = reinterpret_cast < refArray2 > (fiber[0][i]);
          refArray3 Ffiberk = reinterpret_cast < refArray3 > (Ffiber[0][i]);
          refArray4 c1k = reinterpret_cast < refArray4 > (c1[i]);
          refArray5 one_over_hk =
            reinterpret_cast < refArray5 > (one_over_h[i]);
          refArray6 cell_volumek =
            reinterpret_cast < refArray6 > (cell_volume[i]);
          Muscle_Forces < __m256, float[8], 8 > (fk, fiberk, Ffiberk, c1k,
                                                 one_over_hk, cell_volumek);
        }
      for (int __a = 0; __a < 3; __a++)
        for (int __b = 0; __b < 8; __b++)
          for (int __c = 0; __c < 8; __c++)
            if (f[__a][__b][__c] != f_reference[__a][__b][__c])
              {
                std::cerr << "Mismatch detected in AVX implementation" << std::
                  endl;
                std::cerr << "Variable f:" << std::endl;
                std::
                  cerr << "seed=" << seed << ", __a=" << __a << ", __b=" << __b
                  << ", __c=" << __c << std::endl;
                std::cerr << "f AVX=  " << f[__a][__b][__c] << std::endl;
                std::
                  cerr << "f Reference=  " << f_reference[__a][__b][__c] <<
                  std::endl;
                return 1;
              }

    }

  }



  std::cout << "SIMD check successful!" << std::endl;

  return 0;

}

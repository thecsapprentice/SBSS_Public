
struct NEOHOOKEAN_TAG;
struct COROTATED_TAG;

#include <cstdlib>
#include <immintrin.h>
#include <iomanip>
#include <iostream>
#include <sys/time.h>

#include "Muscle_Forces.h"

#define NUM_TRIALS 1000000

template < class T > T Get_Random (const T a = (T) - 1., const T b = (T) 1.)
{
  return ((b - a) * (T) rand ()) / (T) RAND_MAX + a;
}

struct timeval starttime, stoptime;
void
start_timer ()
{
  gettimeofday (&starttime, NULL);
}

void
stop_timer ()
{
  gettimeofday (&stoptime, NULL);
}

double
get_time ()
{
  return (double) stoptime.tv_sec - (double) starttime.tv_sec +
    (double) 1e-6 *(double) stoptime.tv_usec -
    (double) 1e-6 *(double) starttime.tv_usec;
}

int
main (int argc, char *argv[])
{
  typedef float T;
  typedef NEOHOOKEAN_TAG T_MATERIAL_TAG;

  std::cout << "Preparing to Run " << NUM_TRIALS << " of all kernels." << std::
    endl;

  int seed = 1;
  if (argc == 2)
    seed = atoi (argv[1]);
  srand (seed);



  {
    std::cout << "Running Stream Test for Muscle_Forces " << std::endl;

    T f[3][8][8];
    T fiber[3][8];
    T Ffiber[3][8];
    T c1[8];
    T one_over_h[8];
    T cell_volume[8];


    for (int __a = 0; __a < 3; __a++)
      for (int __b = 0; __b < 8; __b++)
        for (int __c = 0; __c < 8; __c++)
          f[__a][__b][__c] = Get_Random < float >();
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

    {
      std::cout << "	Running " << NUM_TRIALS << " of SCALAR :  ";
      start_timer ();
      for (int n = 0; n < NUM_TRIALS; n++)
        {
          for (int i = 0; i < 8; i++)
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
        }
      stop_timer ();
      std::cout << get_time () << "s" << std::endl;
    }

    {
      std::cout << "	Running " << NUM_TRIALS << " of SSE :  ";
      start_timer ();
      for (int n = 0; n < NUM_TRIALS; n++)
        {
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
        }
      stop_timer ();
      std::cout << get_time () << "s" << std::endl;
    }

    {
      std::cout << "	Running " << NUM_TRIALS << " of AVX :  ";
      start_timer ();
      for (int n = 0; n < NUM_TRIALS; n++)
        {
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
        }
      stop_timer ();
      std::cout << get_time () << "s" << std::endl;
    }


  }



  return 0;

}

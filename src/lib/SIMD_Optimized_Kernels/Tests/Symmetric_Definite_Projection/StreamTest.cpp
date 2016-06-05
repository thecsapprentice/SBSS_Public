
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <sys/time.h>
#include "KernelCommon.h"
struct NEOHOOKEAN_TAG;
struct COROTATED_TAG;

#include "Symmetric_Definite_Projection.h"

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

  std::cout << "Preparing to Run " << NUM_TRIALS << " of all kernels." << std::
    endl;

  int seed = 1;
  if (argc == 2)
    seed = atoi (argv[1]);
  srand (seed);



  {
    std::
      cout << "Running Stream Test for Symmetric_Definite_Projection " << std::
      endl;

//=======================================================
//
//        DEFINE ALL VARIABLES USED BY KERNEL
//
//=======================================================

    T A[6][16] __attribute__ ((aligned (64)));
    T Apd[6][16] __attribute__ ((aligned (64)));
    T Apd_reference[6][16] __attribute__ ((aligned (64)));
    T Apd_original[6][16] __attribute__ ((aligned (64)));


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
//             COMPUTE SCALAR RESULTS
//
//=======================================================

    {
      std::cout << "	Running " << NUM_TRIALS << " of SCALAR :  ";
      start_timer ();
      for (int n = 0; n < NUM_TRIALS; n++)
        {
          typedef T (&refArray1)[6][16];
          typedef T (&refArray2)[6][16];
          for (int i = 0; i < 16; i += 1)
            {
              refArray1 Ak = reinterpret_cast < refArray1 > (A[0][i]);
              refArray2 Apdk = reinterpret_cast < refArray2 > (Apd[0][i]);
              Symmetric_Definite_Projection < float, float[16], int[16] > (Ak,
                                                                           Apdk);
            }

        }
      stop_timer ();
      std::cout << get_time () << "s" << std::endl;
    }

//=======================================================
//
//             COMPUTE SSE RESULTS
//
//=======================================================

#ifdef ENABLE_SSE_INSTRUCTION_SET
    {
      std::cout << "	Running " << NUM_TRIALS << " of SSE :  ";
      start_timer ();
      for (int n = 0; n < NUM_TRIALS; n++)
        {
          typedef T (&refArray1)[6][16];
          typedef T (&refArray2)[6][16];
          for (int i = 0; i < 16; i += 4)
            {
              refArray1 Ak = reinterpret_cast < refArray1 > (A[0][i]);
              refArray2 Apdk = reinterpret_cast < refArray2 > (Apd[0][i]);
              Symmetric_Definite_Projection < __m128, float[16], int[16] > (Ak,
                                                                            Apdk);
            }

        }
      stop_timer ();
      std::cout << get_time () << "s" << std::endl;
    }
#endif

//=======================================================
//
//             COMPUTE AVX RESULTS
//
//=======================================================

#ifdef ENABLE_AVX_INSTRUCTION_SET
    {
      std::cout << "	Running " << NUM_TRIALS << " of AVX :  ";
      start_timer ();
      for (int n = 0; n < NUM_TRIALS; n++)
        {
          typedef T (&refArray1)[6][16];
          typedef T (&refArray2)[6][16];
          for (int i = 0; i < 16; i += 8)
            {
              refArray1 Ak = reinterpret_cast < refArray1 > (A[0][i]);
              refArray2 Apdk = reinterpret_cast < refArray2 > (Apd[0][i]);
              Symmetric_Definite_Projection < __m256, float[16], int[16] > (Ak,
                                                                            Apdk);
            }

        }
      stop_timer ();
      std::cout << get_time () << "s" << std::endl;
    }
#endif

//=======================================================
//
//             COMPUTE NEON RESULTS
//
//=======================================================

#ifdef ENABLE_NEON_INSTRUCTION_SET
    {
      std::cout << "	Running " << NUM_TRIALS << " of NEON :  ";
      start_timer ();
      for (int n = 0; n < NUM_TRIALS; n++)
        {
          typedef T (&refArray1)[6][16];
          typedef T (&refArray2)[6][16];
          for (int i = 0; i < 16; i += 4)
            {
              refArray1 Ak = reinterpret_cast < refArray1 > (A[0][i]);
              refArray2 Apdk = reinterpret_cast < refArray2 > (Apd[0][i]);
              Symmetric_Definite_Projection < float32x4_t, float[16],
                int[16] > (Ak, Apdk);
            }

        }
      stop_timer ();
      std::cout << get_time () << "s" << std::endl;
    }
#endif

//=======================================================
//
//             COMPUTE MIC RESULTS
//
//=======================================================

#ifdef ENABLE_MIC_INSTRUCTION_SET
    {
      std::cout << "	Running " << NUM_TRIALS << " of MIC :  ";
      start_timer ();
      for (int n = 0; n < NUM_TRIALS; n++)
        {
          typedef T (&refArray1)[6][16];
          typedef T (&refArray2)[6][16];
          for (int i = 0; i < 16; i += 16)
            {
              refArray1 Ak = reinterpret_cast < refArray1 > (A[0][i]);
              refArray2 Apdk = reinterpret_cast < refArray2 > (Apd[0][i]);
              Symmetric_Definite_Projection < __m512, float[16], int[16] > (Ak,
                                                                            Apdk);
            }

        }
      stop_timer ();
      std::cout << get_time () << "s" << std::endl;
    }
#endif

  }



  return 0;

}

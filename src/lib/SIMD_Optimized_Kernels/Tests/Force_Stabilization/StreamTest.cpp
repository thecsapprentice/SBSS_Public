
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <sys/time.h>
#include "KernelCommon.h"
struct NEOHOOKEAN_TAG;
struct COROTATED_TAG;

#include "Force_Stabilization.h"

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

float Du[3][8][16] __attribute__ ((aligned (64)));
float constant[16] __attribute__ ((aligned (64)));
float dH[3][8][16] __attribute__ ((aligned (64)));
float dH_reference[3][8][16] __attribute__ ((aligned (64)));
float dH_original[3][8][16] __attribute__ ((aligned (64)));

#pragma omp threadprivate(Du, constant, dH, dH_reference, dH_original)

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
    std::cout << "Running Stream Test for Force_Stabilization " << std::endl;

//=======================================================
//
//        DEFINE ALL VARIABLES USED BY KERNEL
//
//=======================================================

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
//             COMPUTE SCALAR RESULTS
//
//=======================================================

    {
      std::cout << "	Running " << NUM_TRIALS << " of SCALAR :  ";
      for( int k= 0; k < 20; k++){
      start_timer ();
#pragma omp parallel for copyin(dH)
      for (int n = 0; n < NUM_TRIALS; n++)
        {
          typedef T (&refArray1)[3][8][16];
          typedef T (&refArray2)[16];
          typedef T (&refArray3)[3][8][16];
          for (int i = 0; i < 16; i += 1)
            {
              refArray1 Duk = reinterpret_cast < refArray1 > (Du[0][0][i]);
              refArray2 constantk =
                reinterpret_cast < refArray2 > (constant[i]);
              refArray3 dHk = reinterpret_cast < refArray3 > (dH[0][0][i]);
              Force_Stabilization < float, float[16], int[16] > (Duk, constantk,
                                                                 dHk);
            }

        }
      stop_timer ();
      std::cout << get_time () << "s" << std::endl;
      }
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
#pragma omp parallel for copyin(dH)
      for (int n = 0; n < NUM_TRIALS; n++)
        {
          typedef T (&refArray1)[3][8][16];
          typedef T (&refArray2)[16];
          typedef T (&refArray3)[3][8][16];
          for (int i = 0; i < 16; i += 4)
            {
              refArray1 Duk = reinterpret_cast < refArray1 > (Du[0][0][i]);
              refArray2 constantk =
                reinterpret_cast < refArray2 > (constant[i]);
              refArray3 dHk = reinterpret_cast < refArray3 > (dH[0][0][i]);
              Force_Stabilization < __m128, float[16], int[16] > (Duk,
                                                                  constantk,
                                                                  dHk);
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
#pragma omp parallel for copyin(dH)
      for (int n = 0; n < NUM_TRIALS; n++)
        {
          typedef T (&refArray1)[3][8][16];
          typedef T (&refArray2)[16];
          typedef T (&refArray3)[3][8][16];
          for (int i = 0; i < 16; i += 8)
            {
              refArray1 Duk = reinterpret_cast < refArray1 > (Du[0][0][i]);
              refArray2 constantk =
                reinterpret_cast < refArray2 > (constant[i]);
              refArray3 dHk = reinterpret_cast < refArray3 > (dH[0][0][i]);
              Force_Stabilization < __m256, float[16], int[16] > (Duk,
                                                                  constantk,
                                                                  dHk);
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
#pragma omp parallel for copyin(dH)
      for (int n = 0; n < NUM_TRIALS; n++)
        {
          typedef T (&refArray1)[3][8][16];
          typedef T (&refArray2)[16];
          typedef T (&refArray3)[3][8][16];
          for (int i = 0; i < 16; i += 4)
            {
              refArray1 Duk = reinterpret_cast < refArray1 > (Du[0][0][i]);
              refArray2 constantk =
                reinterpret_cast < refArray2 > (constant[i]);
              refArray3 dHk = reinterpret_cast < refArray3 > (dH[0][0][i]);
              Force_Stabilization < float32x4_t, float[16], int[16] > (Duk,
                                                                       constantk,
                                                                       dHk);
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
      for( int k= 0; k < 20; k++){
      start_timer ();
#pragma omp parallel for copyin(dH)
      for (int n = 0; n < NUM_TRIALS; n++)
        {
          typedef T (&refArray1)[3][8][16];
          typedef T (&refArray2)[16];
          typedef T (&refArray3)[3][8][16];
          for (int i = 0; i < 16; i += 16)
            {
              refArray1 Duk = reinterpret_cast < refArray1 > (Du[0][0][i]);
              refArray2 constantk =
                reinterpret_cast < refArray2 > (constant[i]);
              refArray3 dHk = reinterpret_cast < refArray3 > (dH[0][0][i]);
              Force_Stabilization < __m512, float[16], int[16] > (Duk,
                                                                  constantk,
                                                                  dHk);
            }

        }
      stop_timer ();
      std::cout << get_time () << "s" << std::endl;
      }
    }
#endif

  }



  return 0;

}

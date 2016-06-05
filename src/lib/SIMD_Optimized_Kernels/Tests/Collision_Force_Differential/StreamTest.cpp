
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <sys/time.h>
#include "KernelCommon.h"
struct NEOHOOKEAN_TAG;
struct COROTATED_TAG;

#include "Collision_Force_Differential.h"

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
      cout << "Running Stream Test for Collision_Force_Differential " << std::
      endl;

//=======================================================
//
//        DEFINE ALL VARIABLES USED BY KERNEL
//
//=======================================================

    T f[3][8][16] __attribute__ ((aligned (64)));
    T f_reference[3][8][16] __attribute__ ((aligned (64)));
    T f_original[3][8][16] __attribute__ ((aligned (64)));
    T u[3][8][16] __attribute__ ((aligned (64)));
    T W[3][16] __attribute__ ((aligned (64)));
    int spring_id[16] __attribute__ ((aligned (64)));
    float *extern_collision_stiffness __attribute__ ((aligned (64)));


    for (int __a = 0; __a < 3; __a++)
      for (int __b = 0; __b < 8; __b++)
        for (int __c = 0; __c < 16; __c++)
          {
            f_original[__a][__b][__c] = Get_Random < float >();
            f[__a][__b][__c] = f_original[__a][__b][__c];
            f_reference[__a][__b][__c] = f_original[__a][__b][__c];
          }
    for (int __a = 0; __a < 3; __a++)
      for (int __b = 0; __b < 8; __b++)
        for (int __c = 0; __c < 16; __c++)
          u[__a][__b][__c] = Get_Random < float >();
    for (int __a = 0; __a < 3; __a++)
      for (int __b = 0; __b < 16; __b++)
        W[__a][__b] = Get_Random < float >();
    for (int __a = 0; __a < 16; __a++)
      spring_id[__a] = Get_Random < int >(1, 99);
    extern_collision_stiffness = new float[100];
    for (int __x__ = 0; __x__ < 100; __x__++)
      extern_collision_stiffness[__x__] = Get_Random < float >();;

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
          typedef T (&refArray1)[3][8][16];
          typedef T (&refArray2)[3][8][16];
          typedef T (&refArray3)[3][16];
          typedef int (&refArray4)[16];
          typedef float *(&refArray5);
          for (int i = 0; i < 16; i += 1)
            {
              refArray1 fk = reinterpret_cast < refArray1 > (f[0][0][i]);
              refArray2 uk = reinterpret_cast < refArray2 > (u[0][0][i]);
              refArray3 Wk = reinterpret_cast < refArray3 > (W[0][i]);
              refArray4 spring_idk =
                reinterpret_cast < refArray4 > (spring_id[i]);
              refArray5 extern_collision_stiffnessk =
                reinterpret_cast < refArray5 > (extern_collision_stiffness);
              Collision_Force_Differential < float, float[16], int[16] > (fk,
                                                                          uk,
                                                                          Wk,
                                                                          spring_idk,
                                                                          extern_collision_stiffnessk);
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
          typedef T (&refArray1)[3][8][16];
          typedef T (&refArray2)[3][8][16];
          typedef T (&refArray3)[3][16];
          typedef int (&refArray4)[16];
          typedef float *(&refArray5);
          for (int i = 0; i < 16; i += 4)
            {
              refArray1 fk = reinterpret_cast < refArray1 > (f[0][0][i]);
              refArray2 uk = reinterpret_cast < refArray2 > (u[0][0][i]);
              refArray3 Wk = reinterpret_cast < refArray3 > (W[0][i]);
              refArray4 spring_idk =
                reinterpret_cast < refArray4 > (spring_id[i]);
              refArray5 extern_collision_stiffnessk =
                reinterpret_cast < refArray5 > (extern_collision_stiffness);
              Collision_Force_Differential < __m128, float[16], int[16] > (fk,
                                                                           uk,
                                                                           Wk,
                                                                           spring_idk,
                                                                           extern_collision_stiffnessk);
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
          typedef T (&refArray1)[3][8][16];
          typedef T (&refArray2)[3][8][16];
          typedef T (&refArray3)[3][16];
          typedef int (&refArray4)[16];
          typedef float *(&refArray5);
          for (int i = 0; i < 16; i += 8)
            {
              refArray1 fk = reinterpret_cast < refArray1 > (f[0][0][i]);
              refArray2 uk = reinterpret_cast < refArray2 > (u[0][0][i]);
              refArray3 Wk = reinterpret_cast < refArray3 > (W[0][i]);
              refArray4 spring_idk =
                reinterpret_cast < refArray4 > (spring_id[i]);
              refArray5 extern_collision_stiffnessk =
                reinterpret_cast < refArray5 > (extern_collision_stiffness);
              Collision_Force_Differential < __m256, float[16], int[16] > (fk,
                                                                           uk,
                                                                           Wk,
                                                                           spring_idk,
                                                                           extern_collision_stiffnessk);
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
          typedef T (&refArray1)[3][8][16];
          typedef T (&refArray2)[3][8][16];
          typedef T (&refArray3)[3][16];
          typedef int (&refArray4)[16];
          typedef float *(&refArray5);
          for (int i = 0; i < 16; i += 4)
            {
              refArray1 fk = reinterpret_cast < refArray1 > (f[0][0][i]);
              refArray2 uk = reinterpret_cast < refArray2 > (u[0][0][i]);
              refArray3 Wk = reinterpret_cast < refArray3 > (W[0][i]);
              refArray4 spring_idk =
                reinterpret_cast < refArray4 > (spring_id[i]);
              refArray5 extern_collision_stiffnessk =
                reinterpret_cast < refArray5 > (extern_collision_stiffness);
              Collision_Force_Differential < float32x4_t, float[16],
                int[16] > (fk, uk, Wk, spring_idk, extern_collision_stiffnessk);
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
          typedef T (&refArray1)[3][8][16];
          typedef T (&refArray2)[3][8][16];
          typedef T (&refArray3)[3][16];
          typedef int (&refArray4)[16];
          typedef float *(&refArray5);
          for (int i = 0; i < 16; i += 16)
            {
              refArray1 fk = reinterpret_cast < refArray1 > (f[0][0][i]);
              refArray2 uk = reinterpret_cast < refArray2 > (u[0][0][i]);
              refArray3 Wk = reinterpret_cast < refArray3 > (W[0][i]);
              refArray4 spring_idk =
                reinterpret_cast < refArray4 > (spring_id[i]);
              refArray5 extern_collision_stiffnessk =
                reinterpret_cast < refArray5 > (extern_collision_stiffness);
              Collision_Force_Differential < __m512, float[16], int[16] > (fk,
                                                                           uk,
                                                                           Wk,
                                                                           spring_idk,
                                                                           extern_collision_stiffnessk);
            }

        }
      stop_timer ();
      std::cout << get_time () << "s" << std::endl;
    }
#endif

  }



  return 0;

}


#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <sys/time.h>
#include "KernelCommon.h"
struct NEOHOOKEAN_TAG;
struct COROTATED_TAG;

#include "Singular_Value_Decomposition.h"

#include <Thread_Queueing/PTHREAD_QUEUE.h>
#include <Kernel_Serial_Base_Helper.h>


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


template < int SIZE > class Singular_Value_Decomposition_SCALAR_None
{
private:
  // Generate Variables Here
  float *_local_A;
  float *_local_U;
  float *_local_S;
  float *_local_V;

public:
    explicit Singular_Value_Decomposition_SCALAR_None (float *A_in, float *U_in,
                                                       float *S_in,
                                                       float
                                                       *V_in):_local_A (A_in),
    _local_U (U_in), _local_S (S_in), _local_V (V_in)
  {
  }
  void Execute (int index)
  {
    // full array typedefs
    //typedef int (&refArray)[SIZE][3][8];
    typedef float (&fullArray1)[SIZE][9][16];
    typedef float (&fullArray2)[SIZE][9][16];
    typedef float (&fullArray3)[SIZE][3][16];
    typedef float (&fullArray4)[SIZE][9][16];

    // full array extractions
    //refArray _rA = reinterpret_cast < refArray >(*_A);
    fullArray1 _rA = reinterpret_cast < fullArray1 > (*_local_A);
    fullArray2 _rU = reinterpret_cast < fullArray2 > (*_local_U);
    fullArray3 _rS = reinterpret_cast < fullArray3 > (*_local_S);
    fullArray4 _rV = reinterpret_cast < fullArray4 > (*_local_V);

    const int ChunkSize = 1;
    // chunk typedef
    //typedef int (&refArray1)[3][8];
    typedef float (&refArray1)[9][16];
    typedef float (&refArray2)[9][16];
    typedef float (&refArray3)[3][16];
    typedef float (&refArray4)[9][16];

    for (int chunk_offset = 0; chunk_offset < 16; chunk_offset += ChunkSize)
      {
        //refArray1 Ak = reinterpret_cast < refArray1 > (_rA[entry][0][chunk_offset]);
        refArray1 Ak =
          reinterpret_cast < refArray1 > (_rA[index][0][chunk_offset]);
        refArray2 Uk =
          reinterpret_cast < refArray2 > (_rU[index][0][chunk_offset]);
        refArray3 Sk =
          reinterpret_cast < refArray3 > (_rS[index][0][chunk_offset]);
        refArray4 Vk =
          reinterpret_cast < refArray4 > (_rV[index][0][chunk_offset]);

        Singular_Value_Decomposition < float, float[16], int[16] > (Ak, Uk, Sk,
                                                                    Vk);
      }

  }
};


#ifdef ENABLE_SSE_INSTRUCTION_SET

template < int SIZE > class Singular_Value_Decomposition_SSE_None
{
private:
  // Generate Variables Here
  float *_local_A;
  float *_local_U;
  float *_local_S;
  float *_local_V;

public:
    explicit Singular_Value_Decomposition_SSE_None (float *A_in, float *U_in,
                                                    float *S_in,
                                                    float
                                                    *V_in):_local_A (A_in),
    _local_U (U_in), _local_S (S_in), _local_V (V_in)
  {
  }
  void Execute (int index)
  {
    // full array typedefs
    //typedef int (&refArray)[SIZE][3][8];
    typedef float (&fullArray1)[SIZE][9][16];
    typedef float (&fullArray2)[SIZE][9][16];
    typedef float (&fullArray3)[SIZE][3][16];
    typedef float (&fullArray4)[SIZE][9][16];

    // full array extractions
    //refArray _rA = reinterpret_cast < refArray >(*_A);
    fullArray1 _rA = reinterpret_cast < fullArray1 > (*_local_A);
    fullArray2 _rU = reinterpret_cast < fullArray2 > (*_local_U);
    fullArray3 _rS = reinterpret_cast < fullArray3 > (*_local_S);
    fullArray4 _rV = reinterpret_cast < fullArray4 > (*_local_V);

    const int ChunkSize = 4;
    // chunk typedef
    //typedef int (&refArray1)[3][8];
    typedef float (&refArray1)[9][16];
    typedef float (&refArray2)[9][16];
    typedef float (&refArray3)[3][16];
    typedef float (&refArray4)[9][16];

    for (int chunk_offset = 0; chunk_offset < 16; chunk_offset += ChunkSize)
      {
        //refArray1 Ak = reinterpret_cast < refArray1 > (_rA[entry][0][chunk_offset]);
        refArray1 Ak =
          reinterpret_cast < refArray1 > (_rA[index][0][chunk_offset]);
        refArray2 Uk =
          reinterpret_cast < refArray2 > (_rU[index][0][chunk_offset]);
        refArray3 Sk =
          reinterpret_cast < refArray3 > (_rS[index][0][chunk_offset]);
        refArray4 Vk =
          reinterpret_cast < refArray4 > (_rV[index][0][chunk_offset]);

        Singular_Value_Decomposition < __m128, float[16], int[16] > (Ak, Uk, Sk,
                                                                     Vk);
      }

  }
};

#endif

#ifdef ENABLE_AVX_INSTRUCTION_SET

template < int SIZE > class Singular_Value_Decomposition_AVX_None
{
private:
  // Generate Variables Here
  float *_local_A;
  float *_local_U;
  float *_local_S;
  float *_local_V;

public:
    explicit Singular_Value_Decomposition_AVX_None (float *A_in, float *U_in,
                                                    float *S_in,
                                                    float
                                                    *V_in):_local_A (A_in),
    _local_U (U_in), _local_S (S_in), _local_V (V_in)
  {
  }
  void Execute (int index)
  {
    // full array typedefs
    //typedef int (&refArray)[SIZE][3][8];
    typedef float (&fullArray1)[SIZE][9][16];
    typedef float (&fullArray2)[SIZE][9][16];
    typedef float (&fullArray3)[SIZE][3][16];
    typedef float (&fullArray4)[SIZE][9][16];

    // full array extractions
    //refArray _rA = reinterpret_cast < refArray >(*_A);
    fullArray1 _rA = reinterpret_cast < fullArray1 > (*_local_A);
    fullArray2 _rU = reinterpret_cast < fullArray2 > (*_local_U);
    fullArray3 _rS = reinterpret_cast < fullArray3 > (*_local_S);
    fullArray4 _rV = reinterpret_cast < fullArray4 > (*_local_V);

    const int ChunkSize = 8;
    // chunk typedef
    //typedef int (&refArray1)[3][8];
    typedef float (&refArray1)[9][16];
    typedef float (&refArray2)[9][16];
    typedef float (&refArray3)[3][16];
    typedef float (&refArray4)[9][16];

    for (int chunk_offset = 0; chunk_offset < 16; chunk_offset += ChunkSize)
      {
        //refArray1 Ak = reinterpret_cast < refArray1 > (_rA[entry][0][chunk_offset]);
        refArray1 Ak =
          reinterpret_cast < refArray1 > (_rA[index][0][chunk_offset]);
        refArray2 Uk =
          reinterpret_cast < refArray2 > (_rU[index][0][chunk_offset]);
        refArray3 Sk =
          reinterpret_cast < refArray3 > (_rS[index][0][chunk_offset]);
        refArray4 Vk =
          reinterpret_cast < refArray4 > (_rV[index][0][chunk_offset]);

        Singular_Value_Decomposition < __m256, float[16], int[16] > (Ak, Uk, Sk,
                                                                     Vk);
      }

  }
};

#endif

#ifdef ENABLE_NEON_INSTRUCTION_SET

template < int SIZE > class Singular_Value_Decomposition_NEON_None
{
private:
  // Generate Variables Here
  float *_local_A;
  float *_local_U;
  float *_local_S;
  float *_local_V;

public:
    explicit Singular_Value_Decomposition_NEON_None (float *A_in, float *U_in,
                                                     float *S_in,
                                                     float
                                                     *V_in):_local_A (A_in),
    _local_U (U_in), _local_S (S_in), _local_V (V_in)
  {
  }
  void Execute (int index)
  {
    // full array typedefs
    //typedef int (&refArray)[SIZE][3][8];
    typedef float (&fullArray1)[SIZE][9][16];
    typedef float (&fullArray2)[SIZE][9][16];
    typedef float (&fullArray3)[SIZE][3][16];
    typedef float (&fullArray4)[SIZE][9][16];

    // full array extractions
    //refArray _rA = reinterpret_cast < refArray >(*_A);
    fullArray1 _rA = reinterpret_cast < fullArray1 > (*_local_A);
    fullArray2 _rU = reinterpret_cast < fullArray2 > (*_local_U);
    fullArray3 _rS = reinterpret_cast < fullArray3 > (*_local_S);
    fullArray4 _rV = reinterpret_cast < fullArray4 > (*_local_V);

    const int ChunkSize = 4;
    // chunk typedef
    //typedef int (&refArray1)[3][8];
    typedef float (&refArray1)[9][16];
    typedef float (&refArray2)[9][16];
    typedef float (&refArray3)[3][16];
    typedef float (&refArray4)[9][16];

    for (int chunk_offset = 0; chunk_offset < 16; chunk_offset += ChunkSize)
      {
        //refArray1 Ak = reinterpret_cast < refArray1 > (_rA[entry][0][chunk_offset]);
        refArray1 Ak =
          reinterpret_cast < refArray1 > (_rA[index][0][chunk_offset]);
        refArray2 Uk =
          reinterpret_cast < refArray2 > (_rU[index][0][chunk_offset]);
        refArray3 Sk =
          reinterpret_cast < refArray3 > (_rS[index][0][chunk_offset]);
        refArray4 Vk =
          reinterpret_cast < refArray4 > (_rV[index][0][chunk_offset]);

        Singular_Value_Decomposition < float32x4_t, float[16], int[16] > (Ak,
                                                                          Uk,
                                                                          Sk,
                                                                          Vk);
      }

  }
};

#endif

#ifdef ENABLE_MIC_INSTRUCTION_SET

template < int SIZE > class Singular_Value_Decomposition_MIC_None
{
private:
  // Generate Variables Here
  float *_local_A;
  float *_local_U;
  float *_local_S;
  float *_local_V;

public:
    explicit Singular_Value_Decomposition_MIC_None (float *A_in, float *U_in,
                                                    float *S_in,
                                                    float
                                                    *V_in):_local_A (A_in),
    _local_U (U_in), _local_S (S_in), _local_V (V_in)
  {
  }
  void Execute (int index)
  {
    // full array typedefs
    //typedef int (&refArray)[SIZE][3][8];
    typedef float (&fullArray1)[SIZE][9][16];
    typedef float (&fullArray2)[SIZE][9][16];
    typedef float (&fullArray3)[SIZE][3][16];
    typedef float (&fullArray4)[SIZE][9][16];

    // full array extractions
    //refArray _rA = reinterpret_cast < refArray >(*_A);
    fullArray1 _rA = reinterpret_cast < fullArray1 > (*_local_A);
    fullArray2 _rU = reinterpret_cast < fullArray2 > (*_local_U);
    fullArray3 _rS = reinterpret_cast < fullArray3 > (*_local_S);
    fullArray4 _rV = reinterpret_cast < fullArray4 > (*_local_V);

    const int ChunkSize = 16;
    // chunk typedef
    //typedef int (&refArray1)[3][8];
    typedef float (&refArray1)[9][16];
    typedef float (&refArray2)[9][16];
    typedef float (&refArray3)[3][16];
    typedef float (&refArray4)[9][16];

    for (int chunk_offset = 0; chunk_offset < 16; chunk_offset += ChunkSize)
      {
        //refArray1 Ak = reinterpret_cast < refArray1 > (_rA[entry][0][chunk_offset]);
        refArray1 Ak =
          reinterpret_cast < refArray1 > (_rA[index][0][chunk_offset]);
        refArray2 Uk =
          reinterpret_cast < refArray2 > (_rU[index][0][chunk_offset]);
        refArray3 Sk =
          reinterpret_cast < refArray3 > (_rS[index][0][chunk_offset]);
        refArray4 Vk =
          reinterpret_cast < refArray4 > (_rV[index][0][chunk_offset]);

        Singular_Value_Decomposition < __m512, float[16], int[16] > (Ak, Uk, Sk,
                                                                     Vk);
      }

  }
};

#endif

int
main (int argc, char *argv[])
{
  typedef float T;

  int seed = 1;
  int threads = 1;
  int threads_max = 1;
  int passes = 1;
  const int data_size = 1000000;
  if (argc >= 2)
    threads = atoi (argv[1]);
  if (argc >= 3)
    threads_max = atoi (argv[2]);
  if (argc >= 4)
    passes = atoi (argv[3]);
  srand (seed);

  pthread_queue = new PTHREAD_QUEUE (threads_max);

  std::
    cout << "Preparing to Run " << data_size << " of all kernels with " <<
    threads << " threads." << std::endl;



  {
    std::
      cout << "Running Thread Test for Singular_Value_Decomposition " << std::
      endl;

//=======================================================
//
//        DEFINE ALL VARIABLES USED BY KERNEL
//
//=======================================================
    std::cout << "\nAllocating all data: ";
    std::cout.flush ();

    start_timer ();
    typedef T (&A_type)[data_size][9][16];
    A_type A =
      reinterpret_cast < A_type >
      (*((T *) (_mm_malloc (data_size * 9 * 16 * sizeof (T), 64))));
    typedef T (&U_type)[data_size][9][16];
    U_type U =
      reinterpret_cast < U_type >
      (*((T *) (_mm_malloc (data_size * 9 * 16 * sizeof (T), 64))));
    typedef T (&S_type)[data_size][3][16];
    S_type S =
      reinterpret_cast < S_type >
      (*((T *) (_mm_malloc (data_size * 3 * 16 * sizeof (T), 64))));
    typedef T (&V_type)[data_size][9][16];
    V_type V =
      reinterpret_cast < V_type >
      (*((T *) (_mm_malloc (data_size * 9 * 16 * sizeof (T), 64))));


    for (int __a = 0; __a < data_size; __a++)
      for (int __b = 0; __b < 9; __b++)
        for (int __c = 0; __c < 16; __c++)
          {
            A[__a][__b][__c] = Get_Random < float >();
    } for (int __a = 0; __a < data_size; __a++)
      for (int __b = 0; __b < 9; __b++)
        for (int __c = 0; __c < 16; __c++)
          {
            U[__a][__b][__c] = Get_Random < float >();
    } for (int __a = 0; __a < data_size; __a++)
      for (int __b = 0; __b < 3; __b++)
        for (int __c = 0; __c < 16; __c++)
          {
            S[__a][__b][__c] = Get_Random < float >();
    } for (int __a = 0; __a < data_size; __a++)
      for (int __b = 0; __b < 9; __b++)
        for (int __c = 0; __c < 16; __c++)
          {
            V[__a][__b][__c] = Get_Random < float >();
          }
    stop_timer ();

    std::cout << get_time () << "s\n\n" << std::endl;


//=======================================================
//
//             COMPUTE SCALAR RESULTS
//
//=======================================================

    {
      std::cout << "	Running " << data_size << " of SCALAR :  " << std::endl;


      Singular_Value_Decomposition_SCALAR_None < data_size > op ((float *) &A,
                                                                 (float *) &U,
                                                                 (float *) &S,
                                                                 (float *) &V);

      for (int t = threads; t <= threads_max;
           t += std::max < int >(((threads_max - threads) / 30), 1))
        {
          std::cout << "Running Test with " << t << " threads." << std::endl;
          MT_Streaming_Kernels::Kernel_Serial_Base_Helper <
            Singular_Value_Decomposition_SCALAR_None < data_size > >helper (op,
                                                                            data_size,
                                                                            t);

          double min_time = 10000000;
          double max_time = -1;
          double avg_time = 0;

          for (int i = 0; i < passes; i++)
            {
              start_timer ();
              helper.Run_Parallel ();
              stop_timer ();
              std::cout << get_time () << "s" << std::endl;
              min_time = std::min < double >(min_time, get_time ());
              max_time = std::max < double >(max_time, get_time ());
              avg_time += get_time ();
            }
          avg_time = avg_time / passes;
          std::cout << "Min pass time: " << min_time << std::endl;
          std::cout << "Max pass time: " << max_time << std::endl;
          std::cout << "Avg pass time: " << avg_time << std::endl;
        }




    }

//=======================================================
//
//             COMPUTE SSE RESULTS
//
//=======================================================
#ifdef ENABLE_SSE_INSTRUCTION_SET
    {
      std::cout << "	Running " << data_size << " of SSE :  " << std::endl;


      Singular_Value_Decomposition_SSE_None < data_size > op ((float *) &A,
                                                              (float *) &U,
                                                              (float *) &S,
                                                              (float *) &V);

      for (int t = threads; t <= threads_max;
           t += std::max < int >(((threads_max - threads) / 30), 1))
        {
          std::cout << "Running Test with " << t << " threads." << std::endl;
          MT_Streaming_Kernels::Kernel_Serial_Base_Helper <
            Singular_Value_Decomposition_SSE_None < data_size > >helper (op,
                                                                         data_size,
                                                                         t);

          double min_time = 10000000;
          double max_time = -1;
          double avg_time = 0;

          for (int i = 0; i < passes; i++)
            {
              start_timer ();
              helper.Run_Parallel ();
              stop_timer ();
              std::cout << get_time () << "s" << std::endl;
              min_time = std::min < double >(min_time, get_time ());
              max_time = std::max < double >(max_time, get_time ());
              avg_time += get_time ();
            }
          avg_time = avg_time / passes;
          std::cout << "Min pass time: " << min_time << std::endl;
          std::cout << "Max pass time: " << max_time << std::endl;
          std::cout << "Avg pass time: " << avg_time << std::endl;
        }




    }
#endif

//=======================================================
//
//             COMPUTE AVX RESULTS
//
//=======================================================
#ifdef ENABLE_AVX_INSTRUCTION_SET
    {
      std::cout << "	Running " << data_size << " of AVX :  " << std::endl;


      Singular_Value_Decomposition_AVX_None < data_size > op ((float *) &A,
                                                              (float *) &U,
                                                              (float *) &S,
                                                              (float *) &V);

      for (int t = threads; t <= threads_max;
           t += std::max < int >(((threads_max - threads) / 30), 1))
        {
          std::cout << "Running Test with " << t << " threads." << std::endl;
          MT_Streaming_Kernels::Kernel_Serial_Base_Helper <
            Singular_Value_Decomposition_AVX_None < data_size > >helper (op,
                                                                         data_size,
                                                                         t);

          double min_time = 10000000;
          double max_time = -1;
          double avg_time = 0;

          for (int i = 0; i < passes; i++)
            {
              start_timer ();
              helper.Run_Parallel ();
              stop_timer ();
              std::cout << get_time () << "s" << std::endl;
              min_time = std::min < double >(min_time, get_time ());
              max_time = std::max < double >(max_time, get_time ());
              avg_time += get_time ();
            }
          avg_time = avg_time / passes;
          std::cout << "Min pass time: " << min_time << std::endl;
          std::cout << "Max pass time: " << max_time << std::endl;
          std::cout << "Avg pass time: " << avg_time << std::endl;
        }




    }
#endif

//=======================================================
//
//             COMPUTE NEON RESULTS
//
//=======================================================
#ifdef ENABLE_NEON_INSTRUCTION_SET
    {
      std::cout << "	Running " << data_size << " of NEON :  " << std::endl;


      Singular_Value_Decomposition_NEON_None < data_size > op ((float *) &A,
                                                               (float *) &U,
                                                               (float *) &S,
                                                               (float *) &V);

      for (int t = threads; t <= threads_max;
           t += std::max < int >(((threads_max - threads) / 30), 1))
        {
          std::cout << "Running Test with " << t << " threads." << std::endl;
          MT_Streaming_Kernels::Kernel_Serial_Base_Helper <
            Singular_Value_Decomposition_NEON_None < data_size > >helper (op,
                                                                          data_size,
                                                                          t);

          double min_time = 10000000;
          double max_time = -1;
          double avg_time = 0;

          for (int i = 0; i < passes; i++)
            {
              start_timer ();
              helper.Run_Parallel ();
              stop_timer ();
              std::cout << get_time () << "s" << std::endl;
              min_time = std::min < double >(min_time, get_time ());
              max_time = std::max < double >(max_time, get_time ());
              avg_time += get_time ();
            }
          avg_time = avg_time / passes;
          std::cout << "Min pass time: " << min_time << std::endl;
          std::cout << "Max pass time: " << max_time << std::endl;
          std::cout << "Avg pass time: " << avg_time << std::endl;
        }




    }
#endif

//=======================================================
//
//             COMPUTE MIC RESULTS
//
//=======================================================
#ifdef ENABLE_MIC_INSTRUCTION_SET
    {
      std::cout << "	Running " << data_size << " of MIC :  " << std::endl;


      Singular_Value_Decomposition_MIC_None < data_size > op ((float *) &A,
                                                              (float *) &U,
                                                              (float *) &S,
                                                              (float *) &V);

      for (int t = threads; t <= threads_max;
           t += std::max < int >(((threads_max - threads) / 30), 1))
        {
          std::cout << "Running Test with " << t << " threads." << std::endl;
          MT_Streaming_Kernels::Kernel_Serial_Base_Helper <
            Singular_Value_Decomposition_MIC_None < data_size > >helper (op,
                                                                         data_size,
                                                                         t);

          double min_time = 10000000;
          double max_time = -1;
          double avg_time = 0;

          for (int i = 0; i < passes; i++)
            {
              start_timer ();
              helper.Run_Parallel ();
              stop_timer ();
              std::cout << get_time () << "s" << std::endl;
              min_time = std::min < double >(min_time, get_time ());
              max_time = std::max < double >(max_time, get_time ());
              avg_time += get_time ();
            }
          avg_time = avg_time / passes;
          std::cout << "Min pass time: " << min_time << std::endl;
          std::cout << "Max pass time: " << max_time << std::endl;
          std::cout << "Avg pass time: " << avg_time << std::endl;
        }




    }
#endif

//=======================================================
//
//        FREE MEMORY USED BY ALL VARIABLES
//
//=======================================================
    std::cout << "\nFreeing all data: " << std::endl;
    std::cout.flush ();

    _mm_free (reinterpret_cast < void *>(A));
    _mm_free (reinterpret_cast < void *>(U));
    _mm_free (reinterpret_cast < void *>(S));
    _mm_free (reinterpret_cast < void *>(V));


  }


  return 0;

}

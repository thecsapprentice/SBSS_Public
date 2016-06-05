
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <sys/time.h>
#include "KernelCommon.h"
struct NEOHOOKEAN_TAG;
struct COROTATED_TAG;

#include "Weighted_Gradient.h"

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


template < int SIZE > class Weighted_Gradient_SCALAR_None
{
private:
  // Generate Variables Here
  float *_local_u;
  float *_local_F;
  float *_local_W;
  float *_local_one_over_h;

public:
    explicit Weighted_Gradient_SCALAR_None (float *u_in, float *F_in,
                                            float *W_in,
                                            float
                                            *one_over_h_in):_local_u (u_in),
    _local_F (F_in), _local_W (W_in), _local_one_over_h (one_over_h_in)
  {
  }
  void Execute (int index)
  {
    // full array typedefs
    //typedef int (&refArray)[SIZE][3][8];
    typedef float (&fullArray1)[SIZE][3][8][16];
    typedef float (&fullArray2)[SIZE][9][16];
    typedef float (&fullArray3)[SIZE][3][16];
    typedef float (&fullArray4)[SIZE][16];

    // full array extractions
    //refArray _rA = reinterpret_cast < refArray >(*_A);
    fullArray1 _ru = reinterpret_cast < fullArray1 > (*_local_u);
    fullArray2 _rF = reinterpret_cast < fullArray2 > (*_local_F);
    fullArray3 _rW = reinterpret_cast < fullArray3 > (*_local_W);
    fullArray4 _rone_over_h =
      reinterpret_cast < fullArray4 > (*_local_one_over_h);

    const int ChunkSize = 1;
    // chunk typedef
    //typedef int (&refArray1)[3][8];
    typedef float (&refArray1)[3][8][16];
    typedef float (&refArray2)[9][16];
    typedef float (&refArray3)[3][16];
    typedef float (&refArray4)[16];

    for (int chunk_offset = 0; chunk_offset < 16; chunk_offset += ChunkSize)
      {
        //refArray1 Ak = reinterpret_cast < refArray1 > (_rA[entry][0][chunk_offset]);
        refArray1 uk =
          reinterpret_cast < refArray1 > (_ru[index][0][0][chunk_offset]);
        refArray2 Fk =
          reinterpret_cast < refArray2 > (_rF[index][0][chunk_offset]);
        refArray3 Wk =
          reinterpret_cast < refArray3 > (_rW[index][0][chunk_offset]);
        refArray4 one_over_hk =
          reinterpret_cast < refArray4 > (_rone_over_h[index][chunk_offset]);

        Weighted_Gradient < float, float[16], int[16] > (uk, Fk, Wk,
                                                         one_over_hk);
      }

  }
};


#ifdef ENABLE_SSE_INSTRUCTION_SET

template < int SIZE > class Weighted_Gradient_SSE_None
{
private:
  // Generate Variables Here
  float *_local_u;
  float *_local_F;
  float *_local_W;
  float *_local_one_over_h;

public:
    explicit Weighted_Gradient_SSE_None (float *u_in, float *F_in, float *W_in,
                                         float *one_over_h_in):_local_u (u_in),
    _local_F (F_in), _local_W (W_in), _local_one_over_h (one_over_h_in)
  {
  }
  void Execute (int index)
  {
    // full array typedefs
    //typedef int (&refArray)[SIZE][3][8];
    typedef float (&fullArray1)[SIZE][3][8][16];
    typedef float (&fullArray2)[SIZE][9][16];
    typedef float (&fullArray3)[SIZE][3][16];
    typedef float (&fullArray4)[SIZE][16];

    // full array extractions
    //refArray _rA = reinterpret_cast < refArray >(*_A);
    fullArray1 _ru = reinterpret_cast < fullArray1 > (*_local_u);
    fullArray2 _rF = reinterpret_cast < fullArray2 > (*_local_F);
    fullArray3 _rW = reinterpret_cast < fullArray3 > (*_local_W);
    fullArray4 _rone_over_h =
      reinterpret_cast < fullArray4 > (*_local_one_over_h);

    const int ChunkSize = 4;
    // chunk typedef
    //typedef int (&refArray1)[3][8];
    typedef float (&refArray1)[3][8][16];
    typedef float (&refArray2)[9][16];
    typedef float (&refArray3)[3][16];
    typedef float (&refArray4)[16];

    for (int chunk_offset = 0; chunk_offset < 16; chunk_offset += ChunkSize)
      {
        //refArray1 Ak = reinterpret_cast < refArray1 > (_rA[entry][0][chunk_offset]);
        refArray1 uk =
          reinterpret_cast < refArray1 > (_ru[index][0][0][chunk_offset]);
        refArray2 Fk =
          reinterpret_cast < refArray2 > (_rF[index][0][chunk_offset]);
        refArray3 Wk =
          reinterpret_cast < refArray3 > (_rW[index][0][chunk_offset]);
        refArray4 one_over_hk =
          reinterpret_cast < refArray4 > (_rone_over_h[index][chunk_offset]);

        Weighted_Gradient < __m128, float[16], int[16] > (uk, Fk, Wk,
                                                          one_over_hk);
      }

  }
};

#endif

#ifdef ENABLE_AVX_INSTRUCTION_SET

template < int SIZE > class Weighted_Gradient_AVX_None
{
private:
  // Generate Variables Here
  float *_local_u;
  float *_local_F;
  float *_local_W;
  float *_local_one_over_h;

public:
    explicit Weighted_Gradient_AVX_None (float *u_in, float *F_in, float *W_in,
                                         float *one_over_h_in):_local_u (u_in),
    _local_F (F_in), _local_W (W_in), _local_one_over_h (one_over_h_in)
  {
  }
  void Execute (int index)
  {
    // full array typedefs
    //typedef int (&refArray)[SIZE][3][8];
    typedef float (&fullArray1)[SIZE][3][8][16];
    typedef float (&fullArray2)[SIZE][9][16];
    typedef float (&fullArray3)[SIZE][3][16];
    typedef float (&fullArray4)[SIZE][16];

    // full array extractions
    //refArray _rA = reinterpret_cast < refArray >(*_A);
    fullArray1 _ru = reinterpret_cast < fullArray1 > (*_local_u);
    fullArray2 _rF = reinterpret_cast < fullArray2 > (*_local_F);
    fullArray3 _rW = reinterpret_cast < fullArray3 > (*_local_W);
    fullArray4 _rone_over_h =
      reinterpret_cast < fullArray4 > (*_local_one_over_h);

    const int ChunkSize = 8;
    // chunk typedef
    //typedef int (&refArray1)[3][8];
    typedef float (&refArray1)[3][8][16];
    typedef float (&refArray2)[9][16];
    typedef float (&refArray3)[3][16];
    typedef float (&refArray4)[16];

    for (int chunk_offset = 0; chunk_offset < 16; chunk_offset += ChunkSize)
      {
        //refArray1 Ak = reinterpret_cast < refArray1 > (_rA[entry][0][chunk_offset]);
        refArray1 uk =
          reinterpret_cast < refArray1 > (_ru[index][0][0][chunk_offset]);
        refArray2 Fk =
          reinterpret_cast < refArray2 > (_rF[index][0][chunk_offset]);
        refArray3 Wk =
          reinterpret_cast < refArray3 > (_rW[index][0][chunk_offset]);
        refArray4 one_over_hk =
          reinterpret_cast < refArray4 > (_rone_over_h[index][chunk_offset]);

        Weighted_Gradient < __m256, float[16], int[16] > (uk, Fk, Wk,
                                                          one_over_hk);
      }

  }
};

#endif

#ifdef ENABLE_NEON_INSTRUCTION_SET

template < int SIZE > class Weighted_Gradient_NEON_None
{
private:
  // Generate Variables Here
  float *_local_u;
  float *_local_F;
  float *_local_W;
  float *_local_one_over_h;

public:
    explicit Weighted_Gradient_NEON_None (float *u_in, float *F_in, float *W_in,
                                          float *one_over_h_in):_local_u (u_in),
    _local_F (F_in), _local_W (W_in), _local_one_over_h (one_over_h_in)
  {
  }
  void Execute (int index)
  {
    // full array typedefs
    //typedef int (&refArray)[SIZE][3][8];
    typedef float (&fullArray1)[SIZE][3][8][16];
    typedef float (&fullArray2)[SIZE][9][16];
    typedef float (&fullArray3)[SIZE][3][16];
    typedef float (&fullArray4)[SIZE][16];

    // full array extractions
    //refArray _rA = reinterpret_cast < refArray >(*_A);
    fullArray1 _ru = reinterpret_cast < fullArray1 > (*_local_u);
    fullArray2 _rF = reinterpret_cast < fullArray2 > (*_local_F);
    fullArray3 _rW = reinterpret_cast < fullArray3 > (*_local_W);
    fullArray4 _rone_over_h =
      reinterpret_cast < fullArray4 > (*_local_one_over_h);

    const int ChunkSize = 4;
    // chunk typedef
    //typedef int (&refArray1)[3][8];
    typedef float (&refArray1)[3][8][16];
    typedef float (&refArray2)[9][16];
    typedef float (&refArray3)[3][16];
    typedef float (&refArray4)[16];

    for (int chunk_offset = 0; chunk_offset < 16; chunk_offset += ChunkSize)
      {
        //refArray1 Ak = reinterpret_cast < refArray1 > (_rA[entry][0][chunk_offset]);
        refArray1 uk =
          reinterpret_cast < refArray1 > (_ru[index][0][0][chunk_offset]);
        refArray2 Fk =
          reinterpret_cast < refArray2 > (_rF[index][0][chunk_offset]);
        refArray3 Wk =
          reinterpret_cast < refArray3 > (_rW[index][0][chunk_offset]);
        refArray4 one_over_hk =
          reinterpret_cast < refArray4 > (_rone_over_h[index][chunk_offset]);

        Weighted_Gradient < float32x4_t, float[16], int[16] > (uk, Fk, Wk,
                                                               one_over_hk);
      }

  }
};

#endif

#ifdef ENABLE_MIC_INSTRUCTION_SET

template < int SIZE > class Weighted_Gradient_MIC_None
{
private:
  // Generate Variables Here
  float *_local_u;
  float *_local_F;
  float *_local_W;
  float *_local_one_over_h;

public:
    explicit Weighted_Gradient_MIC_None (float *u_in, float *F_in, float *W_in,
                                         float *one_over_h_in):_local_u (u_in),
    _local_F (F_in), _local_W (W_in), _local_one_over_h (one_over_h_in)
  {
  }
  void Execute (int index)
  {
    // full array typedefs
    //typedef int (&refArray)[SIZE][3][8];
    typedef float (&fullArray1)[SIZE][3][8][16];
    typedef float (&fullArray2)[SIZE][9][16];
    typedef float (&fullArray3)[SIZE][3][16];
    typedef float (&fullArray4)[SIZE][16];

    // full array extractions
    //refArray _rA = reinterpret_cast < refArray >(*_A);
    fullArray1 _ru = reinterpret_cast < fullArray1 > (*_local_u);
    fullArray2 _rF = reinterpret_cast < fullArray2 > (*_local_F);
    fullArray3 _rW = reinterpret_cast < fullArray3 > (*_local_W);
    fullArray4 _rone_over_h =
      reinterpret_cast < fullArray4 > (*_local_one_over_h);

    const int ChunkSize = 16;
    // chunk typedef
    //typedef int (&refArray1)[3][8];
    typedef float (&refArray1)[3][8][16];
    typedef float (&refArray2)[9][16];
    typedef float (&refArray3)[3][16];
    typedef float (&refArray4)[16];

    for (int chunk_offset = 0; chunk_offset < 16; chunk_offset += ChunkSize)
      {
        //refArray1 Ak = reinterpret_cast < refArray1 > (_rA[entry][0][chunk_offset]);
        refArray1 uk =
          reinterpret_cast < refArray1 > (_ru[index][0][0][chunk_offset]);
        refArray2 Fk =
          reinterpret_cast < refArray2 > (_rF[index][0][chunk_offset]);
        refArray3 Wk =
          reinterpret_cast < refArray3 > (_rW[index][0][chunk_offset]);
        refArray4 one_over_hk =
          reinterpret_cast < refArray4 > (_rone_over_h[index][chunk_offset]);

        Weighted_Gradient < __m512, float[16], int[16] > (uk, Fk, Wk,
                                                          one_over_hk);
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
    std::cout << "Running Thread Test for Weighted_Gradient " << std::endl;

//=======================================================
//
//        DEFINE ALL VARIABLES USED BY KERNEL
//
//=======================================================
    std::cout << "\nAllocating all data: ";
    std::cout.flush ();

    start_timer ();
    typedef T (&u_type)[data_size][3][8][16];
    u_type u =
      reinterpret_cast < u_type >
      (*((T *) (_mm_malloc (data_size * 3 * 8 * 16 * sizeof (T), 64))));
    typedef T (&F_type)[data_size][9][16];
    F_type F =
      reinterpret_cast < F_type >
      (*((T *) (_mm_malloc (data_size * 9 * 16 * sizeof (T), 64))));
    typedef T (&W_type)[data_size][3][16];
    W_type W =
      reinterpret_cast < W_type >
      (*((T *) (_mm_malloc (data_size * 3 * 16 * sizeof (T), 64))));
    typedef T (&one_over_h_type)[data_size][16];
    one_over_h_type one_over_h =
      reinterpret_cast < one_over_h_type >
      (*((T *) (_mm_malloc (data_size * 16 * sizeof (T), 64))));


    for (int __a = 0; __a < data_size; __a++)
      for (int __b = 0; __b < 3; __b++)
        for (int __c = 0; __c < 8; __c++)
          for (int __d = 0; __d < 16; __d++)
            {
              u[__a][__b][__c][__d] = Get_Random < float >();
    } for (int __a = 0; __a < data_size; __a++)
      for (int __b = 0; __b < 9; __b++)
        for (int __c = 0; __c < 16; __c++)
          {
            F[__a][__b][__c] = Get_Random < float >();
    } for (int __a = 0; __a < data_size; __a++)
      for (int __b = 0; __b < 3; __b++)
        for (int __c = 0; __c < 16; __c++)
          {
            W[__a][__b][__c] = Get_Random < float >();
    } for (int __a = 0; __a < data_size; __a++)
      for (int __b = 0; __b < 16; __b++)
        {
          one_over_h[__a][__b] = Get_Random < float >();
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


      Weighted_Gradient_SCALAR_None < data_size > op ((float *) &u,
                                                      (float *) &F,
                                                      (float *) &W,
                                                      (float *) &one_over_h);

      for (int t = threads; t <= threads_max;
           t += std::max < int >(((threads_max - threads) / 30), 1))
        {
          std::cout << "Running Test with " << t << " threads." << std::endl;
          MT_Streaming_Kernels::Kernel_Serial_Base_Helper <
            Weighted_Gradient_SCALAR_None < data_size > >helper (op, data_size,
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


      Weighted_Gradient_SSE_None < data_size > op ((float *) &u, (float *) &F,
                                                   (float *) &W,
                                                   (float *) &one_over_h);

      for (int t = threads; t <= threads_max;
           t += std::max < int >(((threads_max - threads) / 30), 1))
        {
          std::cout << "Running Test with " << t << " threads." << std::endl;
          MT_Streaming_Kernels::Kernel_Serial_Base_Helper <
            Weighted_Gradient_SSE_None < data_size > >helper (op, data_size, t);

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


      Weighted_Gradient_AVX_None < data_size > op ((float *) &u, (float *) &F,
                                                   (float *) &W,
                                                   (float *) &one_over_h);

      for (int t = threads; t <= threads_max;
           t += std::max < int >(((threads_max - threads) / 30), 1))
        {
          std::cout << "Running Test with " << t << " threads." << std::endl;
          MT_Streaming_Kernels::Kernel_Serial_Base_Helper <
            Weighted_Gradient_AVX_None < data_size > >helper (op, data_size, t);

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


      Weighted_Gradient_NEON_None < data_size > op ((float *) &u, (float *) &F,
                                                    (float *) &W,
                                                    (float *) &one_over_h);

      for (int t = threads; t <= threads_max;
           t += std::max < int >(((threads_max - threads) / 30), 1))
        {
          std::cout << "Running Test with " << t << " threads." << std::endl;
          MT_Streaming_Kernels::Kernel_Serial_Base_Helper <
            Weighted_Gradient_NEON_None < data_size > >helper (op, data_size,
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


      Weighted_Gradient_MIC_None < data_size > op ((float *) &u, (float *) &F,
                                                   (float *) &W,
                                                   (float *) &one_over_h);

      for (int t = threads; t <= threads_max;
           t += std::max < int >(((threads_max - threads) / 30), 1))
        {
          std::cout << "Running Test with " << t << " threads." << std::endl;
          MT_Streaming_Kernels::Kernel_Serial_Base_Helper <
            Weighted_Gradient_MIC_None < data_size > >helper (op, data_size, t);

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

    _mm_free (reinterpret_cast < void *>(u));
    _mm_free (reinterpret_cast < void *>(F));
    _mm_free (reinterpret_cast < void *>(W));
    _mm_free (reinterpret_cast < void *>(one_over_h));


  }


  return 0;

}

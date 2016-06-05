
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <sys/time.h>
#include "KernelCommon.h"
struct NEOHOOKEAN_TAG;
struct COROTATED_TAG;

#include "Rotated_Stress_Derivative_Contraction.h"

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


template < int SIZE > class Rotated_Stress_Derivative_Contraction_SCALAR_None
{
private:
  // Generate Variables Here
  float *_local_dPdF;
  float *_local_dF_Hat;
  float *_local_dP_Hat;

public:
    explicit Rotated_Stress_Derivative_Contraction_SCALAR_None (float *dPdF_in,
                                                                float
                                                                *dF_Hat_in,
                                                                float
                                                                *dP_Hat_in):_local_dPdF
    (dPdF_in), _local_dF_Hat (dF_Hat_in), _local_dP_Hat (dP_Hat_in)
  {
  }
  void Execute (int index)
  {
    // full array typedefs
    //typedef int (&refArray)[SIZE][3][8];
    typedef float (&fullArray1)[SIZE][12][16];
    typedef float (&fullArray2)[SIZE][9][16];
    typedef float (&fullArray3)[SIZE][9][16];

    // full array extractions
    //refArray _rA = reinterpret_cast < refArray >(*_A);
    fullArray1 _rdPdF = reinterpret_cast < fullArray1 > (*_local_dPdF);
    fullArray2 _rdF_Hat = reinterpret_cast < fullArray2 > (*_local_dF_Hat);
    fullArray3 _rdP_Hat = reinterpret_cast < fullArray3 > (*_local_dP_Hat);

    const int ChunkSize = 1;
    // chunk typedef
    //typedef int (&refArray1)[3][8];
    typedef float (&refArray1)[12][16];
    typedef float (&refArray2)[9][16];
    typedef float (&refArray3)[9][16];

    for (int chunk_offset = 0; chunk_offset < 16; chunk_offset += ChunkSize)
      {
        //refArray1 Ak = reinterpret_cast < refArray1 > (_rA[entry][0][chunk_offset]);
        refArray1 dPdFk =
          reinterpret_cast < refArray1 > (_rdPdF[index][0][chunk_offset]);
        refArray2 dF_Hatk =
          reinterpret_cast < refArray2 > (_rdF_Hat[index][0][chunk_offset]);
        refArray3 dP_Hatk =
          reinterpret_cast < refArray3 > (_rdP_Hat[index][0][chunk_offset]);

        Rotated_Stress_Derivative_Contraction < float, float[16],
          int[16] > (dPdFk, dF_Hatk, dP_Hatk);
      }

  }
};


#ifdef ENABLE_SSE_INSTRUCTION_SET

template < int SIZE > class Rotated_Stress_Derivative_Contraction_SSE_None
{
private:
  // Generate Variables Here
  float *_local_dPdF;
  float *_local_dF_Hat;
  float *_local_dP_Hat;

public:
    explicit Rotated_Stress_Derivative_Contraction_SSE_None (float *dPdF_in,
                                                             float *dF_Hat_in,
                                                             float
                                                             *dP_Hat_in):_local_dPdF
    (dPdF_in), _local_dF_Hat (dF_Hat_in), _local_dP_Hat (dP_Hat_in)
  {
  }
  void Execute (int index)
  {
    // full array typedefs
    //typedef int (&refArray)[SIZE][3][8];
    typedef float (&fullArray1)[SIZE][12][16];
    typedef float (&fullArray2)[SIZE][9][16];
    typedef float (&fullArray3)[SIZE][9][16];

    // full array extractions
    //refArray _rA = reinterpret_cast < refArray >(*_A);
    fullArray1 _rdPdF = reinterpret_cast < fullArray1 > (*_local_dPdF);
    fullArray2 _rdF_Hat = reinterpret_cast < fullArray2 > (*_local_dF_Hat);
    fullArray3 _rdP_Hat = reinterpret_cast < fullArray3 > (*_local_dP_Hat);

    const int ChunkSize = 4;
    // chunk typedef
    //typedef int (&refArray1)[3][8];
    typedef float (&refArray1)[12][16];
    typedef float (&refArray2)[9][16];
    typedef float (&refArray3)[9][16];

    for (int chunk_offset = 0; chunk_offset < 16; chunk_offset += ChunkSize)
      {
        //refArray1 Ak = reinterpret_cast < refArray1 > (_rA[entry][0][chunk_offset]);
        refArray1 dPdFk =
          reinterpret_cast < refArray1 > (_rdPdF[index][0][chunk_offset]);
        refArray2 dF_Hatk =
          reinterpret_cast < refArray2 > (_rdF_Hat[index][0][chunk_offset]);
        refArray3 dP_Hatk =
          reinterpret_cast < refArray3 > (_rdP_Hat[index][0][chunk_offset]);

        Rotated_Stress_Derivative_Contraction < __m128, float[16],
          int[16] > (dPdFk, dF_Hatk, dP_Hatk);
      }

  }
};

#endif

#ifdef ENABLE_AVX_INSTRUCTION_SET

template < int SIZE > class Rotated_Stress_Derivative_Contraction_AVX_None
{
private:
  // Generate Variables Here
  float *_local_dPdF;
  float *_local_dF_Hat;
  float *_local_dP_Hat;

public:
    explicit Rotated_Stress_Derivative_Contraction_AVX_None (float *dPdF_in,
                                                             float *dF_Hat_in,
                                                             float
                                                             *dP_Hat_in):_local_dPdF
    (dPdF_in), _local_dF_Hat (dF_Hat_in), _local_dP_Hat (dP_Hat_in)
  {
  }
  void Execute (int index)
  {
    // full array typedefs
    //typedef int (&refArray)[SIZE][3][8];
    typedef float (&fullArray1)[SIZE][12][16];
    typedef float (&fullArray2)[SIZE][9][16];
    typedef float (&fullArray3)[SIZE][9][16];

    // full array extractions
    //refArray _rA = reinterpret_cast < refArray >(*_A);
    fullArray1 _rdPdF = reinterpret_cast < fullArray1 > (*_local_dPdF);
    fullArray2 _rdF_Hat = reinterpret_cast < fullArray2 > (*_local_dF_Hat);
    fullArray3 _rdP_Hat = reinterpret_cast < fullArray3 > (*_local_dP_Hat);

    const int ChunkSize = 8;
    // chunk typedef
    //typedef int (&refArray1)[3][8];
    typedef float (&refArray1)[12][16];
    typedef float (&refArray2)[9][16];
    typedef float (&refArray3)[9][16];

    for (int chunk_offset = 0; chunk_offset < 16; chunk_offset += ChunkSize)
      {
        //refArray1 Ak = reinterpret_cast < refArray1 > (_rA[entry][0][chunk_offset]);
        refArray1 dPdFk =
          reinterpret_cast < refArray1 > (_rdPdF[index][0][chunk_offset]);
        refArray2 dF_Hatk =
          reinterpret_cast < refArray2 > (_rdF_Hat[index][0][chunk_offset]);
        refArray3 dP_Hatk =
          reinterpret_cast < refArray3 > (_rdP_Hat[index][0][chunk_offset]);

        Rotated_Stress_Derivative_Contraction < __m256, float[16],
          int[16] > (dPdFk, dF_Hatk, dP_Hatk);
      }

  }
};

#endif

#ifdef ENABLE_NEON_INSTRUCTION_SET

template < int SIZE > class Rotated_Stress_Derivative_Contraction_NEON_None
{
private:
  // Generate Variables Here
  float *_local_dPdF;
  float *_local_dF_Hat;
  float *_local_dP_Hat;

public:
    explicit Rotated_Stress_Derivative_Contraction_NEON_None (float *dPdF_in,
                                                              float *dF_Hat_in,
                                                              float
                                                              *dP_Hat_in):_local_dPdF
    (dPdF_in), _local_dF_Hat (dF_Hat_in), _local_dP_Hat (dP_Hat_in)
  {
  }
  void Execute (int index)
  {
    // full array typedefs
    //typedef int (&refArray)[SIZE][3][8];
    typedef float (&fullArray1)[SIZE][12][16];
    typedef float (&fullArray2)[SIZE][9][16];
    typedef float (&fullArray3)[SIZE][9][16];

    // full array extractions
    //refArray _rA = reinterpret_cast < refArray >(*_A);
    fullArray1 _rdPdF = reinterpret_cast < fullArray1 > (*_local_dPdF);
    fullArray2 _rdF_Hat = reinterpret_cast < fullArray2 > (*_local_dF_Hat);
    fullArray3 _rdP_Hat = reinterpret_cast < fullArray3 > (*_local_dP_Hat);

    const int ChunkSize = 4;
    // chunk typedef
    //typedef int (&refArray1)[3][8];
    typedef float (&refArray1)[12][16];
    typedef float (&refArray2)[9][16];
    typedef float (&refArray3)[9][16];

    for (int chunk_offset = 0; chunk_offset < 16; chunk_offset += ChunkSize)
      {
        //refArray1 Ak = reinterpret_cast < refArray1 > (_rA[entry][0][chunk_offset]);
        refArray1 dPdFk =
          reinterpret_cast < refArray1 > (_rdPdF[index][0][chunk_offset]);
        refArray2 dF_Hatk =
          reinterpret_cast < refArray2 > (_rdF_Hat[index][0][chunk_offset]);
        refArray3 dP_Hatk =
          reinterpret_cast < refArray3 > (_rdP_Hat[index][0][chunk_offset]);

        Rotated_Stress_Derivative_Contraction < float32x4_t, float[16],
          int[16] > (dPdFk, dF_Hatk, dP_Hatk);
      }

  }
};

#endif

#ifdef ENABLE_MIC_INSTRUCTION_SET

template < int SIZE > class Rotated_Stress_Derivative_Contraction_MIC_None
{
private:
  // Generate Variables Here
  float *_local_dPdF;
  float *_local_dF_Hat;
  float *_local_dP_Hat;

public:
    explicit Rotated_Stress_Derivative_Contraction_MIC_None (float *dPdF_in,
                                                             float *dF_Hat_in,
                                                             float
                                                             *dP_Hat_in):_local_dPdF
    (dPdF_in), _local_dF_Hat (dF_Hat_in), _local_dP_Hat (dP_Hat_in)
  {
  }
  void Execute (int index)
  {
    // full array typedefs
    //typedef int (&refArray)[SIZE][3][8];
    typedef float (&fullArray1)[SIZE][12][16];
    typedef float (&fullArray2)[SIZE][9][16];
    typedef float (&fullArray3)[SIZE][9][16];

    // full array extractions
    //refArray _rA = reinterpret_cast < refArray >(*_A);
    fullArray1 _rdPdF = reinterpret_cast < fullArray1 > (*_local_dPdF);
    fullArray2 _rdF_Hat = reinterpret_cast < fullArray2 > (*_local_dF_Hat);
    fullArray3 _rdP_Hat = reinterpret_cast < fullArray3 > (*_local_dP_Hat);

    const int ChunkSize = 16;
    // chunk typedef
    //typedef int (&refArray1)[3][8];
    typedef float (&refArray1)[12][16];
    typedef float (&refArray2)[9][16];
    typedef float (&refArray3)[9][16];

    for (int chunk_offset = 0; chunk_offset < 16; chunk_offset += ChunkSize)
      {
        //refArray1 Ak = reinterpret_cast < refArray1 > (_rA[entry][0][chunk_offset]);
        refArray1 dPdFk =
          reinterpret_cast < refArray1 > (_rdPdF[index][0][chunk_offset]);
        refArray2 dF_Hatk =
          reinterpret_cast < refArray2 > (_rdF_Hat[index][0][chunk_offset]);
        refArray3 dP_Hatk =
          reinterpret_cast < refArray3 > (_rdP_Hat[index][0][chunk_offset]);

        Rotated_Stress_Derivative_Contraction < __m512, float[16],
          int[16] > (dPdFk, dF_Hatk, dP_Hatk);
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
      cout << "Running Thread Test for Rotated_Stress_Derivative_Contraction "
      << std::endl;

//=======================================================
//
//        DEFINE ALL VARIABLES USED BY KERNEL
//
//=======================================================
    std::cout << "\nAllocating all data: ";
    std::cout.flush ();

    start_timer ();
    typedef T (&dPdF_type)[data_size][12][16];
    dPdF_type dPdF =
      reinterpret_cast < dPdF_type >
      (*((T *) (_mm_malloc (data_size * 12 * 16 * sizeof (T), 64))));
    typedef T (&dF_Hat_type)[data_size][9][16];
    dF_Hat_type dF_Hat =
      reinterpret_cast < dF_Hat_type >
      (*((T *) (_mm_malloc (data_size * 9 * 16 * sizeof (T), 64))));
    typedef T (&dP_Hat_type)[data_size][9][16];
    dP_Hat_type dP_Hat =
      reinterpret_cast < dP_Hat_type >
      (*((T *) (_mm_malloc (data_size * 9 * 16 * sizeof (T), 64))));


    for (int __a = 0; __a < data_size; __a++)
      for (int __b = 0; __b < 12; __b++)
        for (int __c = 0; __c < 16; __c++)
          {
            dPdF[__a][__b][__c] = Get_Random < float >();
    } for (int __a = 0; __a < data_size; __a++)
      for (int __b = 0; __b < 9; __b++)
        for (int __c = 0; __c < 16; __c++)
          {
            dF_Hat[__a][__b][__c] = Get_Random < float >();
    } for (int __a = 0; __a < data_size; __a++)
      for (int __b = 0; __b < 9; __b++)
        for (int __c = 0; __c < 16; __c++)
          {
            dP_Hat[__a][__b][__c] = Get_Random < float >();
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


      Rotated_Stress_Derivative_Contraction_SCALAR_None < data_size >
        op ((float *) &dPdF, (float *) &dF_Hat, (float *) &dP_Hat);

      for (int t = threads; t <= threads_max;
           t += std::max < int >(((threads_max - threads) / 30), 1))
        {
          std::cout << "Running Test with " << t << " threads." << std::endl;
          MT_Streaming_Kernels::Kernel_Serial_Base_Helper <
            Rotated_Stress_Derivative_Contraction_SCALAR_None < data_size >
            >helper (op, data_size, t);

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


      Rotated_Stress_Derivative_Contraction_SSE_None < data_size >
        op ((float *) &dPdF, (float *) &dF_Hat, (float *) &dP_Hat);

      for (int t = threads; t <= threads_max;
           t += std::max < int >(((threads_max - threads) / 30), 1))
        {
          std::cout << "Running Test with " << t << " threads." << std::endl;
          MT_Streaming_Kernels::Kernel_Serial_Base_Helper <
            Rotated_Stress_Derivative_Contraction_SSE_None < data_size >
            >helper (op, data_size, t);

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


      Rotated_Stress_Derivative_Contraction_AVX_None < data_size >
        op ((float *) &dPdF, (float *) &dF_Hat, (float *) &dP_Hat);

      for (int t = threads; t <= threads_max;
           t += std::max < int >(((threads_max - threads) / 30), 1))
        {
          std::cout << "Running Test with " << t << " threads." << std::endl;
          MT_Streaming_Kernels::Kernel_Serial_Base_Helper <
            Rotated_Stress_Derivative_Contraction_AVX_None < data_size >
            >helper (op, data_size, t);

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


      Rotated_Stress_Derivative_Contraction_NEON_None < data_size >
        op ((float *) &dPdF, (float *) &dF_Hat, (float *) &dP_Hat);

      for (int t = threads; t <= threads_max;
           t += std::max < int >(((threads_max - threads) / 30), 1))
        {
          std::cout << "Running Test with " << t << " threads." << std::endl;
          MT_Streaming_Kernels::Kernel_Serial_Base_Helper <
            Rotated_Stress_Derivative_Contraction_NEON_None < data_size >
            >helper (op, data_size, t);

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


      Rotated_Stress_Derivative_Contraction_MIC_None < data_size >
        op ((float *) &dPdF, (float *) &dF_Hat, (float *) &dP_Hat);

      for (int t = threads; t <= threads_max;
           t += std::max < int >(((threads_max - threads) / 30), 1))
        {
          std::cout << "Running Test with " << t << " threads." << std::endl;
          MT_Streaming_Kernels::Kernel_Serial_Base_Helper <
            Rotated_Stress_Derivative_Contraction_MIC_None < data_size >
            >helper (op, data_size, t);

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

    _mm_free (reinterpret_cast < void *>(dPdF));
    _mm_free (reinterpret_cast < void *>(dF_Hat));
    _mm_free (reinterpret_cast < void *>(dP_Hat));


  }


  return 0;

}

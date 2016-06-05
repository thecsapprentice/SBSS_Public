
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <sys/time.h>
#include "KernelCommon.h"
struct NEOHOOKEAN_TAG;
struct COROTATED_TAG;

#include "Muscle_Forces.h"

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


template < int SIZE > class Muscle_Forces_SCALAR_None
{
private:
  // Generate Variables Here
  float *_local_f;
  float *_local_fiber;
  float *_local_Ffiber;
  float *_local_c1;
  float *_local_one_over_h;
  float *_local_cell_volume;

public:
    explicit Muscle_Forces_SCALAR_None (float *f_in, float *fiber_in,
                                        float *Ffiber_in, float *c1_in,
                                        float *one_over_h_in,
                                        float *cell_volume_in):_local_f (f_in),
    _local_fiber (fiber_in), _local_Ffiber (Ffiber_in), _local_c1 (c1_in),
    _local_one_over_h (one_over_h_in), _local_cell_volume (cell_volume_in)
  {
  }
  void Execute (int index)
  {
    // full array typedefs
    //typedef int (&refArray)[SIZE][3][8];
    typedef float (&fullArray1)[SIZE][3][8][16];
    typedef float (&fullArray2)[SIZE][3][16];
    typedef float (&fullArray3)[SIZE][3][16];
    typedef float (&fullArray4)[SIZE][16];
    typedef float (&fullArray5)[SIZE][16];
    typedef float (&fullArray6)[SIZE][16];

    // full array extractions
    //refArray _rA = reinterpret_cast < refArray >(*_A);
    fullArray1 _rf = reinterpret_cast < fullArray1 > (*_local_f);
    fullArray2 _rfiber = reinterpret_cast < fullArray2 > (*_local_fiber);
    fullArray3 _rFfiber = reinterpret_cast < fullArray3 > (*_local_Ffiber);
    fullArray4 _rc1 = reinterpret_cast < fullArray4 > (*_local_c1);
    fullArray5 _rone_over_h =
      reinterpret_cast < fullArray5 > (*_local_one_over_h);
    fullArray6 _rcell_volume =
      reinterpret_cast < fullArray6 > (*_local_cell_volume);

    const int ChunkSize = 1;
    // chunk typedef
    //typedef int (&refArray1)[3][8];
    typedef float (&refArray1)[3][8][16];
    typedef float (&refArray2)[3][16];
    typedef float (&refArray3)[3][16];
    typedef float (&refArray4)[16];
    typedef float (&refArray5)[16];
    typedef float (&refArray6)[16];

    for (int chunk_offset = 0; chunk_offset < 16; chunk_offset += ChunkSize)
      {
        //refArray1 Ak = reinterpret_cast < refArray1 > (_rA[entry][0][chunk_offset]);
        refArray1 fk =
          reinterpret_cast < refArray1 > (_rf[index][0][0][chunk_offset]);
        refArray2 fiberk =
          reinterpret_cast < refArray2 > (_rfiber[index][0][chunk_offset]);
        refArray3 Ffiberk =
          reinterpret_cast < refArray3 > (_rFfiber[index][0][chunk_offset]);
        refArray4 c1k =
          reinterpret_cast < refArray4 > (_rc1[index][chunk_offset]);
        refArray5 one_over_hk =
          reinterpret_cast < refArray5 > (_rone_over_h[index][chunk_offset]);
        refArray6 cell_volumek =
          reinterpret_cast < refArray6 > (_rcell_volume[index][chunk_offset]);

        Muscle_Forces < float, float[16], int[16] > (fk, fiberk, Ffiberk, c1k,
                                                     one_over_hk, cell_volumek);
      }

  }
};


#ifdef ENABLE_SSE_INSTRUCTION_SET

template < int SIZE > class Muscle_Forces_SSE_None
{
private:
  // Generate Variables Here
  float *_local_f;
  float *_local_fiber;
  float *_local_Ffiber;
  float *_local_c1;
  float *_local_one_over_h;
  float *_local_cell_volume;

public:
    explicit Muscle_Forces_SSE_None (float *f_in, float *fiber_in,
                                     float *Ffiber_in, float *c1_in,
                                     float *one_over_h_in,
                                     float *cell_volume_in):_local_f (f_in),
    _local_fiber (fiber_in), _local_Ffiber (Ffiber_in), _local_c1 (c1_in),
    _local_one_over_h (one_over_h_in), _local_cell_volume (cell_volume_in)
  {
  }
  void Execute (int index)
  {
    // full array typedefs
    //typedef int (&refArray)[SIZE][3][8];
    typedef float (&fullArray1)[SIZE][3][8][16];
    typedef float (&fullArray2)[SIZE][3][16];
    typedef float (&fullArray3)[SIZE][3][16];
    typedef float (&fullArray4)[SIZE][16];
    typedef float (&fullArray5)[SIZE][16];
    typedef float (&fullArray6)[SIZE][16];

    // full array extractions
    //refArray _rA = reinterpret_cast < refArray >(*_A);
    fullArray1 _rf = reinterpret_cast < fullArray1 > (*_local_f);
    fullArray2 _rfiber = reinterpret_cast < fullArray2 > (*_local_fiber);
    fullArray3 _rFfiber = reinterpret_cast < fullArray3 > (*_local_Ffiber);
    fullArray4 _rc1 = reinterpret_cast < fullArray4 > (*_local_c1);
    fullArray5 _rone_over_h =
      reinterpret_cast < fullArray5 > (*_local_one_over_h);
    fullArray6 _rcell_volume =
      reinterpret_cast < fullArray6 > (*_local_cell_volume);

    const int ChunkSize = 4;
    // chunk typedef
    //typedef int (&refArray1)[3][8];
    typedef float (&refArray1)[3][8][16];
    typedef float (&refArray2)[3][16];
    typedef float (&refArray3)[3][16];
    typedef float (&refArray4)[16];
    typedef float (&refArray5)[16];
    typedef float (&refArray6)[16];

    for (int chunk_offset = 0; chunk_offset < 16; chunk_offset += ChunkSize)
      {
        //refArray1 Ak = reinterpret_cast < refArray1 > (_rA[entry][0][chunk_offset]);
        refArray1 fk =
          reinterpret_cast < refArray1 > (_rf[index][0][0][chunk_offset]);
        refArray2 fiberk =
          reinterpret_cast < refArray2 > (_rfiber[index][0][chunk_offset]);
        refArray3 Ffiberk =
          reinterpret_cast < refArray3 > (_rFfiber[index][0][chunk_offset]);
        refArray4 c1k =
          reinterpret_cast < refArray4 > (_rc1[index][chunk_offset]);
        refArray5 one_over_hk =
          reinterpret_cast < refArray5 > (_rone_over_h[index][chunk_offset]);
        refArray6 cell_volumek =
          reinterpret_cast < refArray6 > (_rcell_volume[index][chunk_offset]);

        Muscle_Forces < __m128, float[16], int[16] > (fk, fiberk, Ffiberk, c1k,
                                                      one_over_hk,
                                                      cell_volumek);
      }

  }
};

#endif

#ifdef ENABLE_AVX_INSTRUCTION_SET

template < int SIZE > class Muscle_Forces_AVX_None
{
private:
  // Generate Variables Here
  float *_local_f;
  float *_local_fiber;
  float *_local_Ffiber;
  float *_local_c1;
  float *_local_one_over_h;
  float *_local_cell_volume;

public:
    explicit Muscle_Forces_AVX_None (float *f_in, float *fiber_in,
                                     float *Ffiber_in, float *c1_in,
                                     float *one_over_h_in,
                                     float *cell_volume_in):_local_f (f_in),
    _local_fiber (fiber_in), _local_Ffiber (Ffiber_in), _local_c1 (c1_in),
    _local_one_over_h (one_over_h_in), _local_cell_volume (cell_volume_in)
  {
  }
  void Execute (int index)
  {
    // full array typedefs
    //typedef int (&refArray)[SIZE][3][8];
    typedef float (&fullArray1)[SIZE][3][8][16];
    typedef float (&fullArray2)[SIZE][3][16];
    typedef float (&fullArray3)[SIZE][3][16];
    typedef float (&fullArray4)[SIZE][16];
    typedef float (&fullArray5)[SIZE][16];
    typedef float (&fullArray6)[SIZE][16];

    // full array extractions
    //refArray _rA = reinterpret_cast < refArray >(*_A);
    fullArray1 _rf = reinterpret_cast < fullArray1 > (*_local_f);
    fullArray2 _rfiber = reinterpret_cast < fullArray2 > (*_local_fiber);
    fullArray3 _rFfiber = reinterpret_cast < fullArray3 > (*_local_Ffiber);
    fullArray4 _rc1 = reinterpret_cast < fullArray4 > (*_local_c1);
    fullArray5 _rone_over_h =
      reinterpret_cast < fullArray5 > (*_local_one_over_h);
    fullArray6 _rcell_volume =
      reinterpret_cast < fullArray6 > (*_local_cell_volume);

    const int ChunkSize = 8;
    // chunk typedef
    //typedef int (&refArray1)[3][8];
    typedef float (&refArray1)[3][8][16];
    typedef float (&refArray2)[3][16];
    typedef float (&refArray3)[3][16];
    typedef float (&refArray4)[16];
    typedef float (&refArray5)[16];
    typedef float (&refArray6)[16];

    for (int chunk_offset = 0; chunk_offset < 16; chunk_offset += ChunkSize)
      {
        //refArray1 Ak = reinterpret_cast < refArray1 > (_rA[entry][0][chunk_offset]);
        refArray1 fk =
          reinterpret_cast < refArray1 > (_rf[index][0][0][chunk_offset]);
        refArray2 fiberk =
          reinterpret_cast < refArray2 > (_rfiber[index][0][chunk_offset]);
        refArray3 Ffiberk =
          reinterpret_cast < refArray3 > (_rFfiber[index][0][chunk_offset]);
        refArray4 c1k =
          reinterpret_cast < refArray4 > (_rc1[index][chunk_offset]);
        refArray5 one_over_hk =
          reinterpret_cast < refArray5 > (_rone_over_h[index][chunk_offset]);
        refArray6 cell_volumek =
          reinterpret_cast < refArray6 > (_rcell_volume[index][chunk_offset]);

        Muscle_Forces < __m256, float[16], int[16] > (fk, fiberk, Ffiberk, c1k,
                                                      one_over_hk,
                                                      cell_volumek);
      }

  }
};

#endif

#ifdef ENABLE_NEON_INSTRUCTION_SET

template < int SIZE > class Muscle_Forces_NEON_None
{
private:
  // Generate Variables Here
  float *_local_f;
  float *_local_fiber;
  float *_local_Ffiber;
  float *_local_c1;
  float *_local_one_over_h;
  float *_local_cell_volume;

public:
    explicit Muscle_Forces_NEON_None (float *f_in, float *fiber_in,
                                      float *Ffiber_in, float *c1_in,
                                      float *one_over_h_in,
                                      float *cell_volume_in):_local_f (f_in),
    _local_fiber (fiber_in), _local_Ffiber (Ffiber_in), _local_c1 (c1_in),
    _local_one_over_h (one_over_h_in), _local_cell_volume (cell_volume_in)
  {
  }
  void Execute (int index)
  {
    // full array typedefs
    //typedef int (&refArray)[SIZE][3][8];
    typedef float (&fullArray1)[SIZE][3][8][16];
    typedef float (&fullArray2)[SIZE][3][16];
    typedef float (&fullArray3)[SIZE][3][16];
    typedef float (&fullArray4)[SIZE][16];
    typedef float (&fullArray5)[SIZE][16];
    typedef float (&fullArray6)[SIZE][16];

    // full array extractions
    //refArray _rA = reinterpret_cast < refArray >(*_A);
    fullArray1 _rf = reinterpret_cast < fullArray1 > (*_local_f);
    fullArray2 _rfiber = reinterpret_cast < fullArray2 > (*_local_fiber);
    fullArray3 _rFfiber = reinterpret_cast < fullArray3 > (*_local_Ffiber);
    fullArray4 _rc1 = reinterpret_cast < fullArray4 > (*_local_c1);
    fullArray5 _rone_over_h =
      reinterpret_cast < fullArray5 > (*_local_one_over_h);
    fullArray6 _rcell_volume =
      reinterpret_cast < fullArray6 > (*_local_cell_volume);

    const int ChunkSize = 4;
    // chunk typedef
    //typedef int (&refArray1)[3][8];
    typedef float (&refArray1)[3][8][16];
    typedef float (&refArray2)[3][16];
    typedef float (&refArray3)[3][16];
    typedef float (&refArray4)[16];
    typedef float (&refArray5)[16];
    typedef float (&refArray6)[16];

    for (int chunk_offset = 0; chunk_offset < 16; chunk_offset += ChunkSize)
      {
        //refArray1 Ak = reinterpret_cast < refArray1 > (_rA[entry][0][chunk_offset]);
        refArray1 fk =
          reinterpret_cast < refArray1 > (_rf[index][0][0][chunk_offset]);
        refArray2 fiberk =
          reinterpret_cast < refArray2 > (_rfiber[index][0][chunk_offset]);
        refArray3 Ffiberk =
          reinterpret_cast < refArray3 > (_rFfiber[index][0][chunk_offset]);
        refArray4 c1k =
          reinterpret_cast < refArray4 > (_rc1[index][chunk_offset]);
        refArray5 one_over_hk =
          reinterpret_cast < refArray5 > (_rone_over_h[index][chunk_offset]);
        refArray6 cell_volumek =
          reinterpret_cast < refArray6 > (_rcell_volume[index][chunk_offset]);

        Muscle_Forces < float32x4_t, float[16], int[16] > (fk, fiberk, Ffiberk,
                                                           c1k, one_over_hk,
                                                           cell_volumek);
      }

  }
};

#endif

#ifdef ENABLE_MIC_INSTRUCTION_SET

template < int SIZE > class Muscle_Forces_MIC_None
{
private:
  // Generate Variables Here
  float *_local_f;
  float *_local_fiber;
  float *_local_Ffiber;
  float *_local_c1;
  float *_local_one_over_h;
  float *_local_cell_volume;

public:
    explicit Muscle_Forces_MIC_None (float *f_in, float *fiber_in,
                                     float *Ffiber_in, float *c1_in,
                                     float *one_over_h_in,
                                     float *cell_volume_in):_local_f (f_in),
    _local_fiber (fiber_in), _local_Ffiber (Ffiber_in), _local_c1 (c1_in),
    _local_one_over_h (one_over_h_in), _local_cell_volume (cell_volume_in)
  {
  }
  void Execute (int index)
  {
    // full array typedefs
    //typedef int (&refArray)[SIZE][3][8];
    typedef float (&fullArray1)[SIZE][3][8][16];
    typedef float (&fullArray2)[SIZE][3][16];
    typedef float (&fullArray3)[SIZE][3][16];
    typedef float (&fullArray4)[SIZE][16];
    typedef float (&fullArray5)[SIZE][16];
    typedef float (&fullArray6)[SIZE][16];

    // full array extractions
    //refArray _rA = reinterpret_cast < refArray >(*_A);
    fullArray1 _rf = reinterpret_cast < fullArray1 > (*_local_f);
    fullArray2 _rfiber = reinterpret_cast < fullArray2 > (*_local_fiber);
    fullArray3 _rFfiber = reinterpret_cast < fullArray3 > (*_local_Ffiber);
    fullArray4 _rc1 = reinterpret_cast < fullArray4 > (*_local_c1);
    fullArray5 _rone_over_h =
      reinterpret_cast < fullArray5 > (*_local_one_over_h);
    fullArray6 _rcell_volume =
      reinterpret_cast < fullArray6 > (*_local_cell_volume);

    const int ChunkSize = 16;
    // chunk typedef
    //typedef int (&refArray1)[3][8];
    typedef float (&refArray1)[3][8][16];
    typedef float (&refArray2)[3][16];
    typedef float (&refArray3)[3][16];
    typedef float (&refArray4)[16];
    typedef float (&refArray5)[16];
    typedef float (&refArray6)[16];

    for (int chunk_offset = 0; chunk_offset < 16; chunk_offset += ChunkSize)
      {
        //refArray1 Ak = reinterpret_cast < refArray1 > (_rA[entry][0][chunk_offset]);
        refArray1 fk =
          reinterpret_cast < refArray1 > (_rf[index][0][0][chunk_offset]);
        refArray2 fiberk =
          reinterpret_cast < refArray2 > (_rfiber[index][0][chunk_offset]);
        refArray3 Ffiberk =
          reinterpret_cast < refArray3 > (_rFfiber[index][0][chunk_offset]);
        refArray4 c1k =
          reinterpret_cast < refArray4 > (_rc1[index][chunk_offset]);
        refArray5 one_over_hk =
          reinterpret_cast < refArray5 > (_rone_over_h[index][chunk_offset]);
        refArray6 cell_volumek =
          reinterpret_cast < refArray6 > (_rcell_volume[index][chunk_offset]);

        Muscle_Forces < __m512, float[16], int[16] > (fk, fiberk, Ffiberk, c1k,
                                                      one_over_hk,
                                                      cell_volumek);
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
    std::cout << "Running Thread Test for Muscle_Forces " << std::endl;

//=======================================================
//
//        DEFINE ALL VARIABLES USED BY KERNEL
//
//=======================================================
    std::cout << "\nAllocating all data: ";
    std::cout.flush ();

    start_timer ();
    typedef T (&f_type)[data_size][3][8][16];
    f_type f =
      reinterpret_cast < f_type >
      (*((T *) (_mm_malloc (data_size * 3 * 8 * 16 * sizeof (T), 64))));
    typedef T (&fiber_type)[data_size][3][16];
    fiber_type fiber =
      reinterpret_cast < fiber_type >
      (*((T *) (_mm_malloc (data_size * 3 * 16 * sizeof (T), 64))));
    typedef T (&Ffiber_type)[data_size][3][16];
    Ffiber_type Ffiber =
      reinterpret_cast < Ffiber_type >
      (*((T *) (_mm_malloc (data_size * 3 * 16 * sizeof (T), 64))));
    typedef T (&c1_type)[data_size][16];
    c1_type c1 =
      reinterpret_cast < c1_type >
      (*((T *) (_mm_malloc (data_size * 16 * sizeof (T), 64))));
    typedef T (&one_over_h_type)[data_size][16];
    one_over_h_type one_over_h =
      reinterpret_cast < one_over_h_type >
      (*((T *) (_mm_malloc (data_size * 16 * sizeof (T), 64))));
    typedef T (&cell_volume_type)[data_size][16];
    cell_volume_type cell_volume =
      reinterpret_cast < cell_volume_type >
      (*((T *) (_mm_malloc (data_size * 16 * sizeof (T), 64))));


    for (int __a = 0; __a < data_size; __a++)
      for (int __b = 0; __b < 3; __b++)
        for (int __c = 0; __c < 8; __c++)
          for (int __d = 0; __d < 16; __d++)
            {
              f[__a][__b][__c][__d] = Get_Random < float >();
    } for (int __a = 0; __a < data_size; __a++)
      for (int __b = 0; __b < 3; __b++)
        for (int __c = 0; __c < 16; __c++)
          {
            fiber[__a][__b][__c] = Get_Random < float >();
    } for (int __a = 0; __a < data_size; __a++)
      for (int __b = 0; __b < 3; __b++)
        for (int __c = 0; __c < 16; __c++)
          {
            Ffiber[__a][__b][__c] = Get_Random < float >();
    } for (int __a = 0; __a < data_size; __a++)
      for (int __b = 0; __b < 16; __b++)
        {
          c1[__a][__b] = Get_Random < float >();
    } for (int __a = 0; __a < data_size; __a++)
      for (int __b = 0; __b < 16; __b++)
        {
          one_over_h[__a][__b] = Get_Random < float >();
    } for (int __a = 0; __a < data_size; __a++)
      for (int __b = 0; __b < 16; __b++)
        {
          cell_volume[__a][__b] = Get_Random < float >();
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


      Muscle_Forces_SCALAR_None < data_size > op ((float *) &f,
                                                  (float *) &fiber,
                                                  (float *) &Ffiber,
                                                  (float *) &c1,
                                                  (float *) &one_over_h,
                                                  (float *) &cell_volume);

      for (int t = threads; t <= threads_max;
           t += std::max < int >(((threads_max - threads) / 30), 1))
        {
          std::cout << "Running Test with " << t << " threads." << std::endl;
          MT_Streaming_Kernels::Kernel_Serial_Base_Helper <
            Muscle_Forces_SCALAR_None < data_size > >helper (op, data_size, t);

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


      Muscle_Forces_SSE_None < data_size > op ((float *) &f, (float *) &fiber,
                                               (float *) &Ffiber, (float *) &c1,
                                               (float *) &one_over_h,
                                               (float *) &cell_volume);

      for (int t = threads; t <= threads_max;
           t += std::max < int >(((threads_max - threads) / 30), 1))
        {
          std::cout << "Running Test with " << t << " threads." << std::endl;
          MT_Streaming_Kernels::Kernel_Serial_Base_Helper <
            Muscle_Forces_SSE_None < data_size > >helper (op, data_size, t);

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


      Muscle_Forces_AVX_None < data_size > op ((float *) &f, (float *) &fiber,
                                               (float *) &Ffiber, (float *) &c1,
                                               (float *) &one_over_h,
                                               (float *) &cell_volume);

      for (int t = threads; t <= threads_max;
           t += std::max < int >(((threads_max - threads) / 30), 1))
        {
          std::cout << "Running Test with " << t << " threads." << std::endl;
          MT_Streaming_Kernels::Kernel_Serial_Base_Helper <
            Muscle_Forces_AVX_None < data_size > >helper (op, data_size, t);

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


      Muscle_Forces_NEON_None < data_size > op ((float *) &f, (float *) &fiber,
                                                (float *) &Ffiber,
                                                (float *) &c1,
                                                (float *) &one_over_h,
                                                (float *) &cell_volume);

      for (int t = threads; t <= threads_max;
           t += std::max < int >(((threads_max - threads) / 30), 1))
        {
          std::cout << "Running Test with " << t << " threads." << std::endl;
          MT_Streaming_Kernels::Kernel_Serial_Base_Helper <
            Muscle_Forces_NEON_None < data_size > >helper (op, data_size, t);

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


      Muscle_Forces_MIC_None < data_size > op ((float *) &f, (float *) &fiber,
                                               (float *) &Ffiber, (float *) &c1,
                                               (float *) &one_over_h,
                                               (float *) &cell_volume);

      for (int t = threads; t <= threads_max;
           t += std::max < int >(((threads_max - threads) / 30), 1))
        {
          std::cout << "Running Test with " << t << " threads." << std::endl;
          MT_Streaming_Kernels::Kernel_Serial_Base_Helper <
            Muscle_Forces_MIC_None < data_size > >helper (op, data_size, t);

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

    _mm_free (reinterpret_cast < void *>(f));
    _mm_free (reinterpret_cast < void *>(fiber));
    _mm_free (reinterpret_cast < void *>(Ffiber));
    _mm_free (reinterpret_cast < void *>(c1));
    _mm_free (reinterpret_cast < void *>(one_over_h));
    _mm_free (reinterpret_cast < void *>(cell_volume));


  }


  return 0;

}

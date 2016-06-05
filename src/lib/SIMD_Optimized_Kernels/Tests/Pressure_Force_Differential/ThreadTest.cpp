
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <sys/time.h>
#include "KernelCommon.h"
struct NEOHOOKEAN_TAG;
struct COROTATED_TAG;

#include "Pressure_Force_Differential.h"

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


template < int SIZE > class Pressure_Force_Differential_SCALAR_None
{
private:
  // Generate Variables Here
  float *_local_dq;
  float *_local_Q_hat;
  float *_local_dF_hat;
  float *_local_dp;
  float *_local_alpha;
  float *_local_alpha_squared_over_kappa;

public:
    explicit Pressure_Force_Differential_SCALAR_None (float *dq_in,
                                                      float *Q_hat_in,
                                                      float *dF_hat_in,
                                                      float *dp_in,
                                                      float *alpha_in,
                                                      float
                                                      *alpha_squared_over_kappa_in):_local_dq
    (dq_in), _local_Q_hat (Q_hat_in), _local_dF_hat (dF_hat_in),
    _local_dp (dp_in), _local_alpha (alpha_in),
    _local_alpha_squared_over_kappa (alpha_squared_over_kappa_in)
  {
  }
  void Execute (int index)
  {
    // full array typedefs
    //typedef int (&refArray)[SIZE][3][8];
    typedef float (&fullArray1)[SIZE][16];
    typedef float (&fullArray2)[SIZE][3][16];
    typedef float (&fullArray3)[SIZE][9][16];
    typedef float (&fullArray4)[SIZE][16];
    typedef float (&fullArray5)[SIZE][16];
    typedef float (&fullArray6)[SIZE][16];

    // full array extractions
    //refArray _rA = reinterpret_cast < refArray >(*_A);
    fullArray1 _rdq = reinterpret_cast < fullArray1 > (*_local_dq);
    fullArray2 _rQ_hat = reinterpret_cast < fullArray2 > (*_local_Q_hat);
    fullArray3 _rdF_hat = reinterpret_cast < fullArray3 > (*_local_dF_hat);
    fullArray4 _rdp = reinterpret_cast < fullArray4 > (*_local_dp);
    fullArray5 _ralpha = reinterpret_cast < fullArray5 > (*_local_alpha);
    fullArray6 _ralpha_squared_over_kappa =
      reinterpret_cast < fullArray6 > (*_local_alpha_squared_over_kappa);

    const int ChunkSize = 1;
    // chunk typedef
    //typedef int (&refArray1)[3][8];
    typedef float (&refArray1)[16];
    typedef float (&refArray2)[3][16];
    typedef float (&refArray3)[9][16];
    typedef float (&refArray4)[16];
    typedef float (&refArray5)[16];
    typedef float (&refArray6)[16];

    for (int chunk_offset = 0; chunk_offset < 16; chunk_offset += ChunkSize)
      {
        //refArray1 Ak = reinterpret_cast < refArray1 > (_rA[entry][0][chunk_offset]);
        refArray1 dqk =
          reinterpret_cast < refArray1 > (_rdq[index][chunk_offset]);
        refArray2 Q_hatk =
          reinterpret_cast < refArray2 > (_rQ_hat[index][0][chunk_offset]);
        refArray3 dF_hatk =
          reinterpret_cast < refArray3 > (_rdF_hat[index][0][chunk_offset]);
        refArray4 dpk =
          reinterpret_cast < refArray4 > (_rdp[index][chunk_offset]);
        refArray5 alphak =
          reinterpret_cast < refArray5 > (_ralpha[index][chunk_offset]);
        refArray6 alpha_squared_over_kappak =
          reinterpret_cast < refArray6 >
          (_ralpha_squared_over_kappa[index][chunk_offset]);

        Pressure_Force_Differential < float, float[16], int[16] > (dqk, Q_hatk,
                                                                   dF_hatk, dpk,
                                                                   alphak,
                                                                   alpha_squared_over_kappak);
      }

  }
};


#ifdef ENABLE_SSE_INSTRUCTION_SET

template < int SIZE > class Pressure_Force_Differential_SSE_None
{
private:
  // Generate Variables Here
  float *_local_dq;
  float *_local_Q_hat;
  float *_local_dF_hat;
  float *_local_dp;
  float *_local_alpha;
  float *_local_alpha_squared_over_kappa;

public:
    explicit Pressure_Force_Differential_SSE_None (float *dq_in,
                                                   float *Q_hat_in,
                                                   float *dF_hat_in,
                                                   float *dp_in,
                                                   float *alpha_in,
                                                   float
                                                   *alpha_squared_over_kappa_in):_local_dq
    (dq_in), _local_Q_hat (Q_hat_in), _local_dF_hat (dF_hat_in),
    _local_dp (dp_in), _local_alpha (alpha_in),
    _local_alpha_squared_over_kappa (alpha_squared_over_kappa_in)
  {
  }
  void Execute (int index)
  {
    // full array typedefs
    //typedef int (&refArray)[SIZE][3][8];
    typedef float (&fullArray1)[SIZE][16];
    typedef float (&fullArray2)[SIZE][3][16];
    typedef float (&fullArray3)[SIZE][9][16];
    typedef float (&fullArray4)[SIZE][16];
    typedef float (&fullArray5)[SIZE][16];
    typedef float (&fullArray6)[SIZE][16];

    // full array extractions
    //refArray _rA = reinterpret_cast < refArray >(*_A);
    fullArray1 _rdq = reinterpret_cast < fullArray1 > (*_local_dq);
    fullArray2 _rQ_hat = reinterpret_cast < fullArray2 > (*_local_Q_hat);
    fullArray3 _rdF_hat = reinterpret_cast < fullArray3 > (*_local_dF_hat);
    fullArray4 _rdp = reinterpret_cast < fullArray4 > (*_local_dp);
    fullArray5 _ralpha = reinterpret_cast < fullArray5 > (*_local_alpha);
    fullArray6 _ralpha_squared_over_kappa =
      reinterpret_cast < fullArray6 > (*_local_alpha_squared_over_kappa);

    const int ChunkSize = 4;
    // chunk typedef
    //typedef int (&refArray1)[3][8];
    typedef float (&refArray1)[16];
    typedef float (&refArray2)[3][16];
    typedef float (&refArray3)[9][16];
    typedef float (&refArray4)[16];
    typedef float (&refArray5)[16];
    typedef float (&refArray6)[16];

    for (int chunk_offset = 0; chunk_offset < 16; chunk_offset += ChunkSize)
      {
        //refArray1 Ak = reinterpret_cast < refArray1 > (_rA[entry][0][chunk_offset]);
        refArray1 dqk =
          reinterpret_cast < refArray1 > (_rdq[index][chunk_offset]);
        refArray2 Q_hatk =
          reinterpret_cast < refArray2 > (_rQ_hat[index][0][chunk_offset]);
        refArray3 dF_hatk =
          reinterpret_cast < refArray3 > (_rdF_hat[index][0][chunk_offset]);
        refArray4 dpk =
          reinterpret_cast < refArray4 > (_rdp[index][chunk_offset]);
        refArray5 alphak =
          reinterpret_cast < refArray5 > (_ralpha[index][chunk_offset]);
        refArray6 alpha_squared_over_kappak =
          reinterpret_cast < refArray6 >
          (_ralpha_squared_over_kappa[index][chunk_offset]);

        Pressure_Force_Differential < __m128, float[16], int[16] > (dqk, Q_hatk,
                                                                    dF_hatk,
                                                                    dpk, alphak,
                                                                    alpha_squared_over_kappak);
      }

  }
};

#endif

#ifdef ENABLE_AVX_INSTRUCTION_SET

template < int SIZE > class Pressure_Force_Differential_AVX_None
{
private:
  // Generate Variables Here
  float *_local_dq;
  float *_local_Q_hat;
  float *_local_dF_hat;
  float *_local_dp;
  float *_local_alpha;
  float *_local_alpha_squared_over_kappa;

public:
    explicit Pressure_Force_Differential_AVX_None (float *dq_in,
                                                   float *Q_hat_in,
                                                   float *dF_hat_in,
                                                   float *dp_in,
                                                   float *alpha_in,
                                                   float
                                                   *alpha_squared_over_kappa_in):_local_dq
    (dq_in), _local_Q_hat (Q_hat_in), _local_dF_hat (dF_hat_in),
    _local_dp (dp_in), _local_alpha (alpha_in),
    _local_alpha_squared_over_kappa (alpha_squared_over_kappa_in)
  {
  }
  void Execute (int index)
  {
    // full array typedefs
    //typedef int (&refArray)[SIZE][3][8];
    typedef float (&fullArray1)[SIZE][16];
    typedef float (&fullArray2)[SIZE][3][16];
    typedef float (&fullArray3)[SIZE][9][16];
    typedef float (&fullArray4)[SIZE][16];
    typedef float (&fullArray5)[SIZE][16];
    typedef float (&fullArray6)[SIZE][16];

    // full array extractions
    //refArray _rA = reinterpret_cast < refArray >(*_A);
    fullArray1 _rdq = reinterpret_cast < fullArray1 > (*_local_dq);
    fullArray2 _rQ_hat = reinterpret_cast < fullArray2 > (*_local_Q_hat);
    fullArray3 _rdF_hat = reinterpret_cast < fullArray3 > (*_local_dF_hat);
    fullArray4 _rdp = reinterpret_cast < fullArray4 > (*_local_dp);
    fullArray5 _ralpha = reinterpret_cast < fullArray5 > (*_local_alpha);
    fullArray6 _ralpha_squared_over_kappa =
      reinterpret_cast < fullArray6 > (*_local_alpha_squared_over_kappa);

    const int ChunkSize = 8;
    // chunk typedef
    //typedef int (&refArray1)[3][8];
    typedef float (&refArray1)[16];
    typedef float (&refArray2)[3][16];
    typedef float (&refArray3)[9][16];
    typedef float (&refArray4)[16];
    typedef float (&refArray5)[16];
    typedef float (&refArray6)[16];

    for (int chunk_offset = 0; chunk_offset < 16; chunk_offset += ChunkSize)
      {
        //refArray1 Ak = reinterpret_cast < refArray1 > (_rA[entry][0][chunk_offset]);
        refArray1 dqk =
          reinterpret_cast < refArray1 > (_rdq[index][chunk_offset]);
        refArray2 Q_hatk =
          reinterpret_cast < refArray2 > (_rQ_hat[index][0][chunk_offset]);
        refArray3 dF_hatk =
          reinterpret_cast < refArray3 > (_rdF_hat[index][0][chunk_offset]);
        refArray4 dpk =
          reinterpret_cast < refArray4 > (_rdp[index][chunk_offset]);
        refArray5 alphak =
          reinterpret_cast < refArray5 > (_ralpha[index][chunk_offset]);
        refArray6 alpha_squared_over_kappak =
          reinterpret_cast < refArray6 >
          (_ralpha_squared_over_kappa[index][chunk_offset]);

        Pressure_Force_Differential < __m256, float[16], int[16] > (dqk, Q_hatk,
                                                                    dF_hatk,
                                                                    dpk, alphak,
                                                                    alpha_squared_over_kappak);
      }

  }
};

#endif

#ifdef ENABLE_NEON_INSTRUCTION_SET

template < int SIZE > class Pressure_Force_Differential_NEON_None
{
private:
  // Generate Variables Here
  float *_local_dq;
  float *_local_Q_hat;
  float *_local_dF_hat;
  float *_local_dp;
  float *_local_alpha;
  float *_local_alpha_squared_over_kappa;

public:
    explicit Pressure_Force_Differential_NEON_None (float *dq_in,
                                                    float *Q_hat_in,
                                                    float *dF_hat_in,
                                                    float *dp_in,
                                                    float *alpha_in,
                                                    float
                                                    *alpha_squared_over_kappa_in):_local_dq
    (dq_in), _local_Q_hat (Q_hat_in), _local_dF_hat (dF_hat_in),
    _local_dp (dp_in), _local_alpha (alpha_in),
    _local_alpha_squared_over_kappa (alpha_squared_over_kappa_in)
  {
  }
  void Execute (int index)
  {
    // full array typedefs
    //typedef int (&refArray)[SIZE][3][8];
    typedef float (&fullArray1)[SIZE][16];
    typedef float (&fullArray2)[SIZE][3][16];
    typedef float (&fullArray3)[SIZE][9][16];
    typedef float (&fullArray4)[SIZE][16];
    typedef float (&fullArray5)[SIZE][16];
    typedef float (&fullArray6)[SIZE][16];

    // full array extractions
    //refArray _rA = reinterpret_cast < refArray >(*_A);
    fullArray1 _rdq = reinterpret_cast < fullArray1 > (*_local_dq);
    fullArray2 _rQ_hat = reinterpret_cast < fullArray2 > (*_local_Q_hat);
    fullArray3 _rdF_hat = reinterpret_cast < fullArray3 > (*_local_dF_hat);
    fullArray4 _rdp = reinterpret_cast < fullArray4 > (*_local_dp);
    fullArray5 _ralpha = reinterpret_cast < fullArray5 > (*_local_alpha);
    fullArray6 _ralpha_squared_over_kappa =
      reinterpret_cast < fullArray6 > (*_local_alpha_squared_over_kappa);

    const int ChunkSize = 4;
    // chunk typedef
    //typedef int (&refArray1)[3][8];
    typedef float (&refArray1)[16];
    typedef float (&refArray2)[3][16];
    typedef float (&refArray3)[9][16];
    typedef float (&refArray4)[16];
    typedef float (&refArray5)[16];
    typedef float (&refArray6)[16];

    for (int chunk_offset = 0; chunk_offset < 16; chunk_offset += ChunkSize)
      {
        //refArray1 Ak = reinterpret_cast < refArray1 > (_rA[entry][0][chunk_offset]);
        refArray1 dqk =
          reinterpret_cast < refArray1 > (_rdq[index][chunk_offset]);
        refArray2 Q_hatk =
          reinterpret_cast < refArray2 > (_rQ_hat[index][0][chunk_offset]);
        refArray3 dF_hatk =
          reinterpret_cast < refArray3 > (_rdF_hat[index][0][chunk_offset]);
        refArray4 dpk =
          reinterpret_cast < refArray4 > (_rdp[index][chunk_offset]);
        refArray5 alphak =
          reinterpret_cast < refArray5 > (_ralpha[index][chunk_offset]);
        refArray6 alpha_squared_over_kappak =
          reinterpret_cast < refArray6 >
          (_ralpha_squared_over_kappa[index][chunk_offset]);

        Pressure_Force_Differential < float32x4_t, float[16], int[16] > (dqk,
                                                                         Q_hatk,
                                                                         dF_hatk,
                                                                         dpk,
                                                                         alphak,
                                                                         alpha_squared_over_kappak);
      }

  }
};

#endif

#ifdef ENABLE_MIC_INSTRUCTION_SET

template < int SIZE > class Pressure_Force_Differential_MIC_None
{
private:
  // Generate Variables Here
  float *_local_dq;
  float *_local_Q_hat;
  float *_local_dF_hat;
  float *_local_dp;
  float *_local_alpha;
  float *_local_alpha_squared_over_kappa;

public:
    explicit Pressure_Force_Differential_MIC_None (float *dq_in,
                                                   float *Q_hat_in,
                                                   float *dF_hat_in,
                                                   float *dp_in,
                                                   float *alpha_in,
                                                   float
                                                   *alpha_squared_over_kappa_in):_local_dq
    (dq_in), _local_Q_hat (Q_hat_in), _local_dF_hat (dF_hat_in),
    _local_dp (dp_in), _local_alpha (alpha_in),
    _local_alpha_squared_over_kappa (alpha_squared_over_kappa_in)
  {
  }
  void Execute (int index)
  {
    // full array typedefs
    //typedef int (&refArray)[SIZE][3][8];
    typedef float (&fullArray1)[SIZE][16];
    typedef float (&fullArray2)[SIZE][3][16];
    typedef float (&fullArray3)[SIZE][9][16];
    typedef float (&fullArray4)[SIZE][16];
    typedef float (&fullArray5)[SIZE][16];
    typedef float (&fullArray6)[SIZE][16];

    // full array extractions
    //refArray _rA = reinterpret_cast < refArray >(*_A);
    fullArray1 _rdq = reinterpret_cast < fullArray1 > (*_local_dq);
    fullArray2 _rQ_hat = reinterpret_cast < fullArray2 > (*_local_Q_hat);
    fullArray3 _rdF_hat = reinterpret_cast < fullArray3 > (*_local_dF_hat);
    fullArray4 _rdp = reinterpret_cast < fullArray4 > (*_local_dp);
    fullArray5 _ralpha = reinterpret_cast < fullArray5 > (*_local_alpha);
    fullArray6 _ralpha_squared_over_kappa =
      reinterpret_cast < fullArray6 > (*_local_alpha_squared_over_kappa);

    const int ChunkSize = 16;
    // chunk typedef
    //typedef int (&refArray1)[3][8];
    typedef float (&refArray1)[16];
    typedef float (&refArray2)[3][16];
    typedef float (&refArray3)[9][16];
    typedef float (&refArray4)[16];
    typedef float (&refArray5)[16];
    typedef float (&refArray6)[16];

    for (int chunk_offset = 0; chunk_offset < 16; chunk_offset += ChunkSize)
      {
        //refArray1 Ak = reinterpret_cast < refArray1 > (_rA[entry][0][chunk_offset]);
        refArray1 dqk =
          reinterpret_cast < refArray1 > (_rdq[index][chunk_offset]);
        refArray2 Q_hatk =
          reinterpret_cast < refArray2 > (_rQ_hat[index][0][chunk_offset]);
        refArray3 dF_hatk =
          reinterpret_cast < refArray3 > (_rdF_hat[index][0][chunk_offset]);
        refArray4 dpk =
          reinterpret_cast < refArray4 > (_rdp[index][chunk_offset]);
        refArray5 alphak =
          reinterpret_cast < refArray5 > (_ralpha[index][chunk_offset]);
        refArray6 alpha_squared_over_kappak =
          reinterpret_cast < refArray6 >
          (_ralpha_squared_over_kappa[index][chunk_offset]);

        Pressure_Force_Differential < __m512, float[16], int[16] > (dqk, Q_hatk,
                                                                    dF_hatk,
                                                                    dpk, alphak,
                                                                    alpha_squared_over_kappak);
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
    std::cout << "Running Thread Test for Pressure_Force_Differential " << std::
      endl;

//=======================================================
//
//        DEFINE ALL VARIABLES USED BY KERNEL
//
//=======================================================
    std::cout << "\nAllocating all data: ";
    std::cout.flush ();

    start_timer ();
    typedef T (&dq_type)[data_size][16];
    dq_type dq =
      reinterpret_cast < dq_type >
      (*((T *) (_mm_malloc (data_size * 16 * sizeof (T), 64))));
    typedef T (&Q_hat_type)[data_size][3][16];
    Q_hat_type Q_hat =
      reinterpret_cast < Q_hat_type >
      (*((T *) (_mm_malloc (data_size * 3 * 16 * sizeof (T), 64))));
    typedef T (&dF_hat_type)[data_size][9][16];
    dF_hat_type dF_hat =
      reinterpret_cast < dF_hat_type >
      (*((T *) (_mm_malloc (data_size * 9 * 16 * sizeof (T), 64))));
    typedef T (&dp_type)[data_size][16];
    dp_type dp =
      reinterpret_cast < dp_type >
      (*((T *) (_mm_malloc (data_size * 16 * sizeof (T), 64))));
    typedef T (&alpha_type)[data_size][16];
    alpha_type alpha =
      reinterpret_cast < alpha_type >
      (*((T *) (_mm_malloc (data_size * 16 * sizeof (T), 64))));
    typedef T (&alpha_squared_over_kappa_type)[data_size][16];
    alpha_squared_over_kappa_type alpha_squared_over_kappa =
      reinterpret_cast < alpha_squared_over_kappa_type >
      (*((T *) (_mm_malloc (data_size * 16 * sizeof (T), 64))));


    for (int __a = 0; __a < data_size; __a++)
      for (int __b = 0; __b < 16; __b++)
        {
          dq[__a][__b] = Get_Random < float >();
    } for (int __a = 0; __a < data_size; __a++)
      for (int __b = 0; __b < 3; __b++)
        for (int __c = 0; __c < 16; __c++)
          {
            Q_hat[__a][__b][__c] = Get_Random < float >();
    } for (int __a = 0; __a < data_size; __a++)
      for (int __b = 0; __b < 9; __b++)
        for (int __c = 0; __c < 16; __c++)
          {
            dF_hat[__a][__b][__c] = Get_Random < float >();
    } for (int __a = 0; __a < data_size; __a++)
      for (int __b = 0; __b < 16; __b++)
        {
          dp[__a][__b] = Get_Random < float >();
    } for (int __a = 0; __a < data_size; __a++)
      for (int __b = 0; __b < 16; __b++)
        {
          alpha[__a][__b] = Get_Random < float >();
    } for (int __a = 0; __a < data_size; __a++)
      for (int __b = 0; __b < 16; __b++)
        {
          alpha_squared_over_kappa[__a][__b] = Get_Random < float >();
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


      Pressure_Force_Differential_SCALAR_None < data_size > op ((float *) &dq,
                                                                (float *)
                                                                &Q_hat,
                                                                (float *)
                                                                &dF_hat,
                                                                (float *) &dp,
                                                                (float *)
                                                                &alpha,
                                                                (float *)
                                                                &alpha_squared_over_kappa);

      for (int t = threads; t <= threads_max;
           t += std::max < int >(((threads_max - threads) / 30), 1))
        {
          std::cout << "Running Test with " << t << " threads." << std::endl;
          MT_Streaming_Kernels::Kernel_Serial_Base_Helper <
            Pressure_Force_Differential_SCALAR_None < data_size > >helper (op,
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


      Pressure_Force_Differential_SSE_None < data_size > op ((float *) &dq,
                                                             (float *) &Q_hat,
                                                             (float *) &dF_hat,
                                                             (float *) &dp,
                                                             (float *) &alpha,
                                                             (float *)
                                                             &alpha_squared_over_kappa);

      for (int t = threads; t <= threads_max;
           t += std::max < int >(((threads_max - threads) / 30), 1))
        {
          std::cout << "Running Test with " << t << " threads." << std::endl;
          MT_Streaming_Kernels::Kernel_Serial_Base_Helper <
            Pressure_Force_Differential_SSE_None < data_size > >helper (op,
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


      Pressure_Force_Differential_AVX_None < data_size > op ((float *) &dq,
                                                             (float *) &Q_hat,
                                                             (float *) &dF_hat,
                                                             (float *) &dp,
                                                             (float *) &alpha,
                                                             (float *)
                                                             &alpha_squared_over_kappa);

      for (int t = threads; t <= threads_max;
           t += std::max < int >(((threads_max - threads) / 30), 1))
        {
          std::cout << "Running Test with " << t << " threads." << std::endl;
          MT_Streaming_Kernels::Kernel_Serial_Base_Helper <
            Pressure_Force_Differential_AVX_None < data_size > >helper (op,
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


      Pressure_Force_Differential_NEON_None < data_size > op ((float *) &dq,
                                                              (float *) &Q_hat,
                                                              (float *) &dF_hat,
                                                              (float *) &dp,
                                                              (float *) &alpha,
                                                              (float *)
                                                              &alpha_squared_over_kappa);

      for (int t = threads; t <= threads_max;
           t += std::max < int >(((threads_max - threads) / 30), 1))
        {
          std::cout << "Running Test with " << t << " threads." << std::endl;
          MT_Streaming_Kernels::Kernel_Serial_Base_Helper <
            Pressure_Force_Differential_NEON_None < data_size > >helper (op,
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


      Pressure_Force_Differential_MIC_None < data_size > op ((float *) &dq,
                                                             (float *) &Q_hat,
                                                             (float *) &dF_hat,
                                                             (float *) &dp,
                                                             (float *) &alpha,
                                                             (float *)
                                                             &alpha_squared_over_kappa);

      for (int t = threads; t <= threads_max;
           t += std::max < int >(((threads_max - threads) / 30), 1))
        {
          std::cout << "Running Test with " << t << " threads." << std::endl;
          MT_Streaming_Kernels::Kernel_Serial_Base_Helper <
            Pressure_Force_Differential_MIC_None < data_size > >helper (op,
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

    _mm_free (reinterpret_cast < void *>(dq));
    _mm_free (reinterpret_cast < void *>(Q_hat));
    _mm_free (reinterpret_cast < void *>(dF_hat));
    _mm_free (reinterpret_cast < void *>(dp));
    _mm_free (reinterpret_cast < void *>(alpha));
    _mm_free (reinterpret_cast < void *>(alpha_squared_over_kappa));


  }


  return 0;

}

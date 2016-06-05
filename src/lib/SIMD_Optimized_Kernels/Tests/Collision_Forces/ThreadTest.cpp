
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <sys/time.h>
#include "KernelCommon.h"
struct NEOHOOKEAN_TAG;
struct COROTATED_TAG;

#include "Collision_Forces.h"

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


template < int SIZE > class Collision_Forces_SCALAR_None
{
private:
  // Generate Variables Here
  float *_local_f;
  float *_local_u;
  float *_local_N;
  float *_local_W;
  float *_local_h;
  int *_local_spring_id;
  int *_local_spring_id_X;
  int *_local_spring_id_Y;
  int *_local_spring_id_Z;
  float **_local_extern_collision_attach;
  float **_local_extern_collision_stiffness;

public:
    explicit Collision_Forces_SCALAR_None (float *f_in, float *u_in,
                                           float *N_in, float *W_in,
                                           float *h_in, int *spring_id_in,
                                           int *spring_id_X_in,
                                           int *spring_id_Y_in,
                                           int *spring_id_Z_in,
                                           float **extern_collision_attach_in,
                                           float
                                           **extern_collision_stiffness_in):_local_f
    (f_in), _local_u (u_in), _local_N (N_in), _local_W (W_in), _local_h (h_in),
    _local_spring_id (spring_id_in), _local_spring_id_X (spring_id_X_in),
    _local_spring_id_Y (spring_id_Y_in), _local_spring_id_Z (spring_id_Z_in),
    _local_extern_collision_attach (extern_collision_attach_in),
    _local_extern_collision_stiffness (extern_collision_stiffness_in)
  {
  }
  void Execute (int index)
  {
    // full array typedefs
    //typedef int (&refArray)[SIZE][3][8];
    typedef float (&fullArray1)[SIZE][3][8][16];
    typedef float (&fullArray2)[SIZE][3][8][16];
    typedef float (&fullArray3)[SIZE][3][16];
    typedef float (&fullArray4)[SIZE][3][16];
    typedef float (&fullArray5)[SIZE][16];
    typedef int (&fullArray6)[SIZE][16];
    typedef int (&fullArray7)[SIZE][16];
    typedef int (&fullArray8)[SIZE][16];
    typedef int (&fullArray9)[SIZE][16];
    typedef float *(&fullArray10)[SIZE];
    typedef float *(&fullArray11)[SIZE];

    // full array extractions
    //refArray _rA = reinterpret_cast < refArray >(*_A);
    fullArray1 _rf = reinterpret_cast < fullArray1 > (*_local_f);
    fullArray2 _ru = reinterpret_cast < fullArray2 > (*_local_u);
    fullArray3 _rN = reinterpret_cast < fullArray3 > (*_local_N);
    fullArray4 _rW = reinterpret_cast < fullArray4 > (*_local_W);
    fullArray5 _rh = reinterpret_cast < fullArray5 > (*_local_h);
    fullArray6 _rspring_id =
      reinterpret_cast < fullArray6 > (*_local_spring_id);
    fullArray7 _rspring_id_X =
      reinterpret_cast < fullArray7 > (*_local_spring_id_X);
    fullArray8 _rspring_id_Y =
      reinterpret_cast < fullArray8 > (*_local_spring_id_Y);
    fullArray9 _rspring_id_Z =
      reinterpret_cast < fullArray9 > (*_local_spring_id_Z);
    fullArray10 _rextern_collision_attach =
      reinterpret_cast < fullArray10 > (*_local_extern_collision_attach);
    fullArray11 _rextern_collision_stiffness =
      reinterpret_cast < fullArray11 > (*_local_extern_collision_stiffness);

    const int ChunkSize = 1;
    // chunk typedef
    //typedef int (&refArray1)[3][8];
    typedef float (&refArray1)[3][8][16];
    typedef float (&refArray2)[3][8][16];
    typedef float (&refArray3)[3][16];
    typedef float (&refArray4)[3][16];
    typedef float (&refArray5)[16];
    typedef int (&refArray6)[16];
    typedef int (&refArray7)[16];
    typedef int (&refArray8)[16];
    typedef int (&refArray9)[16];
    typedef float *(&refArray10);
    typedef float *(&refArray11);

    for (int chunk_offset = 0; chunk_offset < 16; chunk_offset += ChunkSize)
      {
        //refArray1 Ak = reinterpret_cast < refArray1 > (_rA[entry][0][chunk_offset]);
        refArray1 fk =
          reinterpret_cast < refArray1 > (_rf[index][0][0][chunk_offset]);
        refArray2 uk =
          reinterpret_cast < refArray2 > (_ru[index][0][0][chunk_offset]);
        refArray3 Nk =
          reinterpret_cast < refArray3 > (_rN[index][0][chunk_offset]);
        refArray4 Wk =
          reinterpret_cast < refArray4 > (_rW[index][0][chunk_offset]);
        refArray5 hk =
          reinterpret_cast < refArray5 > (_rh[index][chunk_offset]);
        refArray6 spring_idk =
          reinterpret_cast < refArray6 > (_rspring_id[index][chunk_offset]);
        refArray7 spring_id_Xk =
          reinterpret_cast < refArray7 > (_rspring_id_X[index][chunk_offset]);
        refArray8 spring_id_Yk =
          reinterpret_cast < refArray8 > (_rspring_id_Y[index][chunk_offset]);
        refArray9 spring_id_Zk =
          reinterpret_cast < refArray9 > (_rspring_id_Z[index][chunk_offset]);
        refArray10 extern_collision_attachk =
          reinterpret_cast < refArray10 > (_rextern_collision_attach[index]);
        refArray11 extern_collision_stiffnessk =
          reinterpret_cast < refArray11 > (_rextern_collision_stiffness[index]);

        Collision_Forces < float, float[16], int[16] > (fk, uk, Nk, Wk, hk,
                                                        spring_idk,
                                                        spring_id_Xk,
                                                        spring_id_Yk,
                                                        spring_id_Zk,
                                                        extern_collision_attachk,
                                                        extern_collision_stiffnessk);
      }

  }
};


#ifdef ENABLE_SSE_INSTRUCTION_SET

template < int SIZE > class Collision_Forces_SSE_None
{
private:
  // Generate Variables Here
  float *_local_f;
  float *_local_u;
  float *_local_N;
  float *_local_W;
  float *_local_h;
  int *_local_spring_id;
  int *_local_spring_id_X;
  int *_local_spring_id_Y;
  int *_local_spring_id_Z;
  float **_local_extern_collision_attach;
  float **_local_extern_collision_stiffness;

public:
    explicit Collision_Forces_SSE_None (float *f_in, float *u_in, float *N_in,
                                        float *W_in, float *h_in,
                                        int *spring_id_in, int *spring_id_X_in,
                                        int *spring_id_Y_in,
                                        int *spring_id_Z_in,
                                        float **extern_collision_attach_in,
                                        float
                                        **extern_collision_stiffness_in):_local_f
    (f_in), _local_u (u_in), _local_N (N_in), _local_W (W_in), _local_h (h_in),
    _local_spring_id (spring_id_in), _local_spring_id_X (spring_id_X_in),
    _local_spring_id_Y (spring_id_Y_in), _local_spring_id_Z (spring_id_Z_in),
    _local_extern_collision_attach (extern_collision_attach_in),
    _local_extern_collision_stiffness (extern_collision_stiffness_in)
  {
  }
  void Execute (int index)
  {
    // full array typedefs
    //typedef int (&refArray)[SIZE][3][8];
    typedef float (&fullArray1)[SIZE][3][8][16];
    typedef float (&fullArray2)[SIZE][3][8][16];
    typedef float (&fullArray3)[SIZE][3][16];
    typedef float (&fullArray4)[SIZE][3][16];
    typedef float (&fullArray5)[SIZE][16];
    typedef int (&fullArray6)[SIZE][16];
    typedef int (&fullArray7)[SIZE][16];
    typedef int (&fullArray8)[SIZE][16];
    typedef int (&fullArray9)[SIZE][16];
    typedef float *(&fullArray10)[SIZE];
    typedef float *(&fullArray11)[SIZE];

    // full array extractions
    //refArray _rA = reinterpret_cast < refArray >(*_A);
    fullArray1 _rf = reinterpret_cast < fullArray1 > (*_local_f);
    fullArray2 _ru = reinterpret_cast < fullArray2 > (*_local_u);
    fullArray3 _rN = reinterpret_cast < fullArray3 > (*_local_N);
    fullArray4 _rW = reinterpret_cast < fullArray4 > (*_local_W);
    fullArray5 _rh = reinterpret_cast < fullArray5 > (*_local_h);
    fullArray6 _rspring_id =
      reinterpret_cast < fullArray6 > (*_local_spring_id);
    fullArray7 _rspring_id_X =
      reinterpret_cast < fullArray7 > (*_local_spring_id_X);
    fullArray8 _rspring_id_Y =
      reinterpret_cast < fullArray8 > (*_local_spring_id_Y);
    fullArray9 _rspring_id_Z =
      reinterpret_cast < fullArray9 > (*_local_spring_id_Z);
    fullArray10 _rextern_collision_attach =
      reinterpret_cast < fullArray10 > (*_local_extern_collision_attach);
    fullArray11 _rextern_collision_stiffness =
      reinterpret_cast < fullArray11 > (*_local_extern_collision_stiffness);

    const int ChunkSize = 4;
    // chunk typedef
    //typedef int (&refArray1)[3][8];
    typedef float (&refArray1)[3][8][16];
    typedef float (&refArray2)[3][8][16];
    typedef float (&refArray3)[3][16];
    typedef float (&refArray4)[3][16];
    typedef float (&refArray5)[16];
    typedef int (&refArray6)[16];
    typedef int (&refArray7)[16];
    typedef int (&refArray8)[16];
    typedef int (&refArray9)[16];
    typedef float *(&refArray10);
    typedef float *(&refArray11);

    for (int chunk_offset = 0; chunk_offset < 16; chunk_offset += ChunkSize)
      {
        //refArray1 Ak = reinterpret_cast < refArray1 > (_rA[entry][0][chunk_offset]);
        refArray1 fk =
          reinterpret_cast < refArray1 > (_rf[index][0][0][chunk_offset]);
        refArray2 uk =
          reinterpret_cast < refArray2 > (_ru[index][0][0][chunk_offset]);
        refArray3 Nk =
          reinterpret_cast < refArray3 > (_rN[index][0][chunk_offset]);
        refArray4 Wk =
          reinterpret_cast < refArray4 > (_rW[index][0][chunk_offset]);
        refArray5 hk =
          reinterpret_cast < refArray5 > (_rh[index][chunk_offset]);
        refArray6 spring_idk =
          reinterpret_cast < refArray6 > (_rspring_id[index][chunk_offset]);
        refArray7 spring_id_Xk =
          reinterpret_cast < refArray7 > (_rspring_id_X[index][chunk_offset]);
        refArray8 spring_id_Yk =
          reinterpret_cast < refArray8 > (_rspring_id_Y[index][chunk_offset]);
        refArray9 spring_id_Zk =
          reinterpret_cast < refArray9 > (_rspring_id_Z[index][chunk_offset]);
        refArray10 extern_collision_attachk =
          reinterpret_cast < refArray10 > (_rextern_collision_attach[index]);
        refArray11 extern_collision_stiffnessk =
          reinterpret_cast < refArray11 > (_rextern_collision_stiffness[index]);

        Collision_Forces < __m128, float[16], int[16] > (fk, uk, Nk, Wk, hk,
                                                         spring_idk,
                                                         spring_id_Xk,
                                                         spring_id_Yk,
                                                         spring_id_Zk,
                                                         extern_collision_attachk,
                                                         extern_collision_stiffnessk);
      }

  }
};

#endif

#ifdef ENABLE_AVX_INSTRUCTION_SET

template < int SIZE > class Collision_Forces_AVX_None
{
private:
  // Generate Variables Here
  float *_local_f;
  float *_local_u;
  float *_local_N;
  float *_local_W;
  float *_local_h;
  int *_local_spring_id;
  int *_local_spring_id_X;
  int *_local_spring_id_Y;
  int *_local_spring_id_Z;
  float **_local_extern_collision_attach;
  float **_local_extern_collision_stiffness;

public:
    explicit Collision_Forces_AVX_None (float *f_in, float *u_in, float *N_in,
                                        float *W_in, float *h_in,
                                        int *spring_id_in, int *spring_id_X_in,
                                        int *spring_id_Y_in,
                                        int *spring_id_Z_in,
                                        float **extern_collision_attach_in,
                                        float
                                        **extern_collision_stiffness_in):_local_f
    (f_in), _local_u (u_in), _local_N (N_in), _local_W (W_in), _local_h (h_in),
    _local_spring_id (spring_id_in), _local_spring_id_X (spring_id_X_in),
    _local_spring_id_Y (spring_id_Y_in), _local_spring_id_Z (spring_id_Z_in),
    _local_extern_collision_attach (extern_collision_attach_in),
    _local_extern_collision_stiffness (extern_collision_stiffness_in)
  {
  }
  void Execute (int index)
  {
    // full array typedefs
    //typedef int (&refArray)[SIZE][3][8];
    typedef float (&fullArray1)[SIZE][3][8][16];
    typedef float (&fullArray2)[SIZE][3][8][16];
    typedef float (&fullArray3)[SIZE][3][16];
    typedef float (&fullArray4)[SIZE][3][16];
    typedef float (&fullArray5)[SIZE][16];
    typedef int (&fullArray6)[SIZE][16];
    typedef int (&fullArray7)[SIZE][16];
    typedef int (&fullArray8)[SIZE][16];
    typedef int (&fullArray9)[SIZE][16];
    typedef float *(&fullArray10)[SIZE];
    typedef float *(&fullArray11)[SIZE];

    // full array extractions
    //refArray _rA = reinterpret_cast < refArray >(*_A);
    fullArray1 _rf = reinterpret_cast < fullArray1 > (*_local_f);
    fullArray2 _ru = reinterpret_cast < fullArray2 > (*_local_u);
    fullArray3 _rN = reinterpret_cast < fullArray3 > (*_local_N);
    fullArray4 _rW = reinterpret_cast < fullArray4 > (*_local_W);
    fullArray5 _rh = reinterpret_cast < fullArray5 > (*_local_h);
    fullArray6 _rspring_id =
      reinterpret_cast < fullArray6 > (*_local_spring_id);
    fullArray7 _rspring_id_X =
      reinterpret_cast < fullArray7 > (*_local_spring_id_X);
    fullArray8 _rspring_id_Y =
      reinterpret_cast < fullArray8 > (*_local_spring_id_Y);
    fullArray9 _rspring_id_Z =
      reinterpret_cast < fullArray9 > (*_local_spring_id_Z);
    fullArray10 _rextern_collision_attach =
      reinterpret_cast < fullArray10 > (*_local_extern_collision_attach);
    fullArray11 _rextern_collision_stiffness =
      reinterpret_cast < fullArray11 > (*_local_extern_collision_stiffness);

    const int ChunkSize = 8;
    // chunk typedef
    //typedef int (&refArray1)[3][8];
    typedef float (&refArray1)[3][8][16];
    typedef float (&refArray2)[3][8][16];
    typedef float (&refArray3)[3][16];
    typedef float (&refArray4)[3][16];
    typedef float (&refArray5)[16];
    typedef int (&refArray6)[16];
    typedef int (&refArray7)[16];
    typedef int (&refArray8)[16];
    typedef int (&refArray9)[16];
    typedef float *(&refArray10);
    typedef float *(&refArray11);

    for (int chunk_offset = 0; chunk_offset < 16; chunk_offset += ChunkSize)
      {
        //refArray1 Ak = reinterpret_cast < refArray1 > (_rA[entry][0][chunk_offset]);
        refArray1 fk =
          reinterpret_cast < refArray1 > (_rf[index][0][0][chunk_offset]);
        refArray2 uk =
          reinterpret_cast < refArray2 > (_ru[index][0][0][chunk_offset]);
        refArray3 Nk =
          reinterpret_cast < refArray3 > (_rN[index][0][chunk_offset]);
        refArray4 Wk =
          reinterpret_cast < refArray4 > (_rW[index][0][chunk_offset]);
        refArray5 hk =
          reinterpret_cast < refArray5 > (_rh[index][chunk_offset]);
        refArray6 spring_idk =
          reinterpret_cast < refArray6 > (_rspring_id[index][chunk_offset]);
        refArray7 spring_id_Xk =
          reinterpret_cast < refArray7 > (_rspring_id_X[index][chunk_offset]);
        refArray8 spring_id_Yk =
          reinterpret_cast < refArray8 > (_rspring_id_Y[index][chunk_offset]);
        refArray9 spring_id_Zk =
          reinterpret_cast < refArray9 > (_rspring_id_Z[index][chunk_offset]);
        refArray10 extern_collision_attachk =
          reinterpret_cast < refArray10 > (_rextern_collision_attach[index]);
        refArray11 extern_collision_stiffnessk =
          reinterpret_cast < refArray11 > (_rextern_collision_stiffness[index]);

        Collision_Forces < __m256, float[16], int[16] > (fk, uk, Nk, Wk, hk,
                                                         spring_idk,
                                                         spring_id_Xk,
                                                         spring_id_Yk,
                                                         spring_id_Zk,
                                                         extern_collision_attachk,
                                                         extern_collision_stiffnessk);
      }

  }
};

#endif

#ifdef ENABLE_NEON_INSTRUCTION_SET

template < int SIZE > class Collision_Forces_NEON_None
{
private:
  // Generate Variables Here
  float *_local_f;
  float *_local_u;
  float *_local_N;
  float *_local_W;
  float *_local_h;
  int *_local_spring_id;
  int *_local_spring_id_X;
  int *_local_spring_id_Y;
  int *_local_spring_id_Z;
  float **_local_extern_collision_attach;
  float **_local_extern_collision_stiffness;

public:
    explicit Collision_Forces_NEON_None (float *f_in, float *u_in, float *N_in,
                                         float *W_in, float *h_in,
                                         int *spring_id_in, int *spring_id_X_in,
                                         int *spring_id_Y_in,
                                         int *spring_id_Z_in,
                                         float **extern_collision_attach_in,
                                         float
                                         **extern_collision_stiffness_in):_local_f
    (f_in), _local_u (u_in), _local_N (N_in), _local_W (W_in), _local_h (h_in),
    _local_spring_id (spring_id_in), _local_spring_id_X (spring_id_X_in),
    _local_spring_id_Y (spring_id_Y_in), _local_spring_id_Z (spring_id_Z_in),
    _local_extern_collision_attach (extern_collision_attach_in),
    _local_extern_collision_stiffness (extern_collision_stiffness_in)
  {
  }
  void Execute (int index)
  {
    // full array typedefs
    //typedef int (&refArray)[SIZE][3][8];
    typedef float (&fullArray1)[SIZE][3][8][16];
    typedef float (&fullArray2)[SIZE][3][8][16];
    typedef float (&fullArray3)[SIZE][3][16];
    typedef float (&fullArray4)[SIZE][3][16];
    typedef float (&fullArray5)[SIZE][16];
    typedef int (&fullArray6)[SIZE][16];
    typedef int (&fullArray7)[SIZE][16];
    typedef int (&fullArray8)[SIZE][16];
    typedef int (&fullArray9)[SIZE][16];
    typedef float *(&fullArray10)[SIZE];
    typedef float *(&fullArray11)[SIZE];

    // full array extractions
    //refArray _rA = reinterpret_cast < refArray >(*_A);
    fullArray1 _rf = reinterpret_cast < fullArray1 > (*_local_f);
    fullArray2 _ru = reinterpret_cast < fullArray2 > (*_local_u);
    fullArray3 _rN = reinterpret_cast < fullArray3 > (*_local_N);
    fullArray4 _rW = reinterpret_cast < fullArray4 > (*_local_W);
    fullArray5 _rh = reinterpret_cast < fullArray5 > (*_local_h);
    fullArray6 _rspring_id =
      reinterpret_cast < fullArray6 > (*_local_spring_id);
    fullArray7 _rspring_id_X =
      reinterpret_cast < fullArray7 > (*_local_spring_id_X);
    fullArray8 _rspring_id_Y =
      reinterpret_cast < fullArray8 > (*_local_spring_id_Y);
    fullArray9 _rspring_id_Z =
      reinterpret_cast < fullArray9 > (*_local_spring_id_Z);
    fullArray10 _rextern_collision_attach =
      reinterpret_cast < fullArray10 > (*_local_extern_collision_attach);
    fullArray11 _rextern_collision_stiffness =
      reinterpret_cast < fullArray11 > (*_local_extern_collision_stiffness);

    const int ChunkSize = 4;
    // chunk typedef
    //typedef int (&refArray1)[3][8];
    typedef float (&refArray1)[3][8][16];
    typedef float (&refArray2)[3][8][16];
    typedef float (&refArray3)[3][16];
    typedef float (&refArray4)[3][16];
    typedef float (&refArray5)[16];
    typedef int (&refArray6)[16];
    typedef int (&refArray7)[16];
    typedef int (&refArray8)[16];
    typedef int (&refArray9)[16];
    typedef float *(&refArray10);
    typedef float *(&refArray11);

    for (int chunk_offset = 0; chunk_offset < 16; chunk_offset += ChunkSize)
      {
        //refArray1 Ak = reinterpret_cast < refArray1 > (_rA[entry][0][chunk_offset]);
        refArray1 fk =
          reinterpret_cast < refArray1 > (_rf[index][0][0][chunk_offset]);
        refArray2 uk =
          reinterpret_cast < refArray2 > (_ru[index][0][0][chunk_offset]);
        refArray3 Nk =
          reinterpret_cast < refArray3 > (_rN[index][0][chunk_offset]);
        refArray4 Wk =
          reinterpret_cast < refArray4 > (_rW[index][0][chunk_offset]);
        refArray5 hk =
          reinterpret_cast < refArray5 > (_rh[index][chunk_offset]);
        refArray6 spring_idk =
          reinterpret_cast < refArray6 > (_rspring_id[index][chunk_offset]);
        refArray7 spring_id_Xk =
          reinterpret_cast < refArray7 > (_rspring_id_X[index][chunk_offset]);
        refArray8 spring_id_Yk =
          reinterpret_cast < refArray8 > (_rspring_id_Y[index][chunk_offset]);
        refArray9 spring_id_Zk =
          reinterpret_cast < refArray9 > (_rspring_id_Z[index][chunk_offset]);
        refArray10 extern_collision_attachk =
          reinterpret_cast < refArray10 > (_rextern_collision_attach[index]);
        refArray11 extern_collision_stiffnessk =
          reinterpret_cast < refArray11 > (_rextern_collision_stiffness[index]);

        Collision_Forces < float32x4_t, float[16], int[16] > (fk, uk, Nk, Wk,
                                                              hk, spring_idk,
                                                              spring_id_Xk,
                                                              spring_id_Yk,
                                                              spring_id_Zk,
                                                              extern_collision_attachk,
                                                              extern_collision_stiffnessk);
      }

  }
};

#endif

#ifdef ENABLE_MIC_INSTRUCTION_SET

template < int SIZE > class Collision_Forces_MIC_None
{
private:
  // Generate Variables Here
  float *_local_f;
  float *_local_u;
  float *_local_N;
  float *_local_W;
  float *_local_h;
  int *_local_spring_id;
  int *_local_spring_id_X;
  int *_local_spring_id_Y;
  int *_local_spring_id_Z;
  float **_local_extern_collision_attach;
  float **_local_extern_collision_stiffness;

public:
    explicit Collision_Forces_MIC_None (float *f_in, float *u_in, float *N_in,
                                        float *W_in, float *h_in,
                                        int *spring_id_in, int *spring_id_X_in,
                                        int *spring_id_Y_in,
                                        int *spring_id_Z_in,
                                        float **extern_collision_attach_in,
                                        float
                                        **extern_collision_stiffness_in):_local_f
    (f_in), _local_u (u_in), _local_N (N_in), _local_W (W_in), _local_h (h_in),
    _local_spring_id (spring_id_in), _local_spring_id_X (spring_id_X_in),
    _local_spring_id_Y (spring_id_Y_in), _local_spring_id_Z (spring_id_Z_in),
    _local_extern_collision_attach (extern_collision_attach_in),
    _local_extern_collision_stiffness (extern_collision_stiffness_in)
  {
  }
  void Execute (int index)
  {
    // full array typedefs
    //typedef int (&refArray)[SIZE][3][8];
    typedef float (&fullArray1)[SIZE][3][8][16];
    typedef float (&fullArray2)[SIZE][3][8][16];
    typedef float (&fullArray3)[SIZE][3][16];
    typedef float (&fullArray4)[SIZE][3][16];
    typedef float (&fullArray5)[SIZE][16];
    typedef int (&fullArray6)[SIZE][16];
    typedef int (&fullArray7)[SIZE][16];
    typedef int (&fullArray8)[SIZE][16];
    typedef int (&fullArray9)[SIZE][16];
    typedef float *(&fullArray10)[SIZE];
    typedef float *(&fullArray11)[SIZE];

    // full array extractions
    //refArray _rA = reinterpret_cast < refArray >(*_A);
    fullArray1 _rf = reinterpret_cast < fullArray1 > (*_local_f);
    fullArray2 _ru = reinterpret_cast < fullArray2 > (*_local_u);
    fullArray3 _rN = reinterpret_cast < fullArray3 > (*_local_N);
    fullArray4 _rW = reinterpret_cast < fullArray4 > (*_local_W);
    fullArray5 _rh = reinterpret_cast < fullArray5 > (*_local_h);
    fullArray6 _rspring_id =
      reinterpret_cast < fullArray6 > (*_local_spring_id);
    fullArray7 _rspring_id_X =
      reinterpret_cast < fullArray7 > (*_local_spring_id_X);
    fullArray8 _rspring_id_Y =
      reinterpret_cast < fullArray8 > (*_local_spring_id_Y);
    fullArray9 _rspring_id_Z =
      reinterpret_cast < fullArray9 > (*_local_spring_id_Z);
    fullArray10 _rextern_collision_attach =
      reinterpret_cast < fullArray10 > (*_local_extern_collision_attach);
    fullArray11 _rextern_collision_stiffness =
      reinterpret_cast < fullArray11 > (*_local_extern_collision_stiffness);

    const int ChunkSize = 16;
    // chunk typedef
    //typedef int (&refArray1)[3][8];
    typedef float (&refArray1)[3][8][16];
    typedef float (&refArray2)[3][8][16];
    typedef float (&refArray3)[3][16];
    typedef float (&refArray4)[3][16];
    typedef float (&refArray5)[16];
    typedef int (&refArray6)[16];
    typedef int (&refArray7)[16];
    typedef int (&refArray8)[16];
    typedef int (&refArray9)[16];
    typedef float *(&refArray10);
    typedef float *(&refArray11);

    for (int chunk_offset = 0; chunk_offset < 16; chunk_offset += ChunkSize)
      {
        //refArray1 Ak = reinterpret_cast < refArray1 > (_rA[entry][0][chunk_offset]);
        refArray1 fk =
          reinterpret_cast < refArray1 > (_rf[index][0][0][chunk_offset]);
        refArray2 uk =
          reinterpret_cast < refArray2 > (_ru[index][0][0][chunk_offset]);
        refArray3 Nk =
          reinterpret_cast < refArray3 > (_rN[index][0][chunk_offset]);
        refArray4 Wk =
          reinterpret_cast < refArray4 > (_rW[index][0][chunk_offset]);
        refArray5 hk =
          reinterpret_cast < refArray5 > (_rh[index][chunk_offset]);
        refArray6 spring_idk =
          reinterpret_cast < refArray6 > (_rspring_id[index][chunk_offset]);
        refArray7 spring_id_Xk =
          reinterpret_cast < refArray7 > (_rspring_id_X[index][chunk_offset]);
        refArray8 spring_id_Yk =
          reinterpret_cast < refArray8 > (_rspring_id_Y[index][chunk_offset]);
        refArray9 spring_id_Zk =
          reinterpret_cast < refArray9 > (_rspring_id_Z[index][chunk_offset]);
        refArray10 extern_collision_attachk =
          reinterpret_cast < refArray10 > (_rextern_collision_attach[index]);
        refArray11 extern_collision_stiffnessk =
          reinterpret_cast < refArray11 > (_rextern_collision_stiffness[index]);

        Collision_Forces < __m512, float[16], int[16] > (fk, uk, Nk, Wk, hk,
                                                         spring_idk,
                                                         spring_id_Xk,
                                                         spring_id_Yk,
                                                         spring_id_Zk,
                                                         extern_collision_attachk,
                                                         extern_collision_stiffnessk);
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
    std::cout << "Running Thread Test for Collision_Forces " << std::endl;

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
    typedef T (&u_type)[data_size][3][8][16];
    u_type u =
      reinterpret_cast < u_type >
      (*((T *) (_mm_malloc (data_size * 3 * 8 * 16 * sizeof (T), 64))));
    typedef T (&N_type)[data_size][3][16];
    N_type N =
      reinterpret_cast < N_type >
      (*((T *) (_mm_malloc (data_size * 3 * 16 * sizeof (T), 64))));
    typedef T (&W_type)[data_size][3][16];
    W_type W =
      reinterpret_cast < W_type >
      (*((T *) (_mm_malloc (data_size * 3 * 16 * sizeof (T), 64))));
    typedef T (&h_type)[data_size][16];
    h_type h =
      reinterpret_cast < h_type >
      (*((T *) (_mm_malloc (data_size * 16 * sizeof (T), 64))));
    typedef int (&spring_id_type)[data_size][16];
    spring_id_type spring_id =
      reinterpret_cast < spring_id_type >
      (*((T *) (_mm_malloc (data_size * 16 * sizeof (int), 64))));
    typedef int (&spring_id_X_type)[data_size][16];
    spring_id_X_type spring_id_X =
      reinterpret_cast < spring_id_X_type >
      (*((T *) (_mm_malloc (data_size * 16 * sizeof (int), 64))));
    typedef int (&spring_id_Y_type)[data_size][16];
    spring_id_Y_type spring_id_Y =
      reinterpret_cast < spring_id_Y_type >
      (*((T *) (_mm_malloc (data_size * 16 * sizeof (int), 64))));
    typedef int (&spring_id_Z_type)[data_size][16];
    spring_id_Z_type spring_id_Z =
      reinterpret_cast < spring_id_Z_type >
      (*((T *) (_mm_malloc (data_size * 16 * sizeof (int), 64))));
    typedef float *(&extern_collision_attach_type)[data_size];
    extern_collision_attach_type extern_collision_attach =
      reinterpret_cast < extern_collision_attach_type >
      (*((T *) (_mm_malloc (data_size * sizeof (float *), 64))));
    typedef float *(&extern_collision_stiffness_type)[data_size];
    extern_collision_stiffness_type extern_collision_stiffness =
      reinterpret_cast < extern_collision_stiffness_type >
      (*((T *) (_mm_malloc (data_size * sizeof (float *), 64))));


    for (int __a = 0; __a < data_size; __a++)
      for (int __b = 0; __b < 3; __b++)
        for (int __c = 0; __c < 8; __c++)
          for (int __d = 0; __d < 16; __d++)
            {
              f[__a][__b][__c][__d] = Get_Random < float >();
    } for (int __a = 0; __a < data_size; __a++)
      for (int __b = 0; __b < 3; __b++)
        for (int __c = 0; __c < 8; __c++)
          for (int __d = 0; __d < 16; __d++)
            {
              u[__a][__b][__c][__d] = Get_Random < float >();
    } for (int __a = 0; __a < data_size; __a++)
      for (int __b = 0; __b < 3; __b++)
        for (int __c = 0; __c < 16; __c++)
          {
            N[__a][__b][__c] = Get_Random < float >();
    } for (int __a = 0; __a < data_size; __a++)
      for (int __b = 0; __b < 3; __b++)
        for (int __c = 0; __c < 16; __c++)
          {
            W[__a][__b][__c] = Get_Random < float >();
    } for (int __a = 0; __a < data_size; __a++)
      for (int __b = 0; __b < 16; __b++)
        {
          h[__a][__b] = Get_Random < float >();
    } for (int __a = 0; __a < data_size; __a++)
      for (int __b = 0; __b < 16; __b++)
        {
          spring_id[__a][__b] = Get_Random < int >(1, 99);
    } for (int __a = 0; __a < data_size; __a++)
      for (int __b = 0; __b < 16; __b++)
        {
          spring_id_X[__a][__b] = Get_Random < int >(1, 99);
    } for (int __a = 0; __a < data_size; __a++)
      for (int __b = 0; __b < 16; __b++)
        {
          spring_id_Y[__a][__b] = Get_Random < int >(1, 99);
    } for (int __a = 0; __a < data_size; __a++)
      for (int __b = 0; __b < 16; __b++)
        {
          spring_id_Z[__a][__b] = Get_Random < int >(1, 99);
    } for (int __a = 0; __a < data_size; __a++)
      {
        extern_collision_attach[__a] = new float[100];
        for (int __x__ = 0; __x__ < 100; __x__++)
          extern_collision_attach[__a][__x__] = Get_Random < float >();;
    } for (int __a = 0; __a < data_size; __a++)
      {
        extern_collision_stiffness[__a] = new float[100];
        for (int __x__ = 0; __x__ < 100; __x__++)
          extern_collision_stiffness[__a][__x__] = Get_Random < float >();;
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


      Collision_Forces_SCALAR_None < data_size > op ((float *) &f, (float *) &u,
                                                     (float *) &N, (float *) &W,
                                                     (float *) &h,
                                                     (int *) &spring_id,
                                                     (int *) &spring_id_X,
                                                     (int *) &spring_id_Y,
                                                     (int *) &spring_id_Z,
                                                     (float **)
                                                     &extern_collision_attach,
                                                     (float **)
                                                     &extern_collision_stiffness);

      for (int t = threads; t <= threads_max;
           t += std::max < int >(((threads_max - threads) / 30), 1))
        {
          std::cout << "Running Test with " << t << " threads." << std::endl;
          MT_Streaming_Kernels::Kernel_Serial_Base_Helper <
            Collision_Forces_SCALAR_None < data_size > >helper (op, data_size,
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


      Collision_Forces_SSE_None < data_size > op ((float *) &f, (float *) &u,
                                                  (float *) &N, (float *) &W,
                                                  (float *) &h,
                                                  (int *) &spring_id,
                                                  (int *) &spring_id_X,
                                                  (int *) &spring_id_Y,
                                                  (int *) &spring_id_Z,
                                                  (float **)
                                                  &extern_collision_attach,
                                                  (float **)
                                                  &extern_collision_stiffness);

      for (int t = threads; t <= threads_max;
           t += std::max < int >(((threads_max - threads) / 30), 1))
        {
          std::cout << "Running Test with " << t << " threads." << std::endl;
          MT_Streaming_Kernels::Kernel_Serial_Base_Helper <
            Collision_Forces_SSE_None < data_size > >helper (op, data_size, t);

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


      Collision_Forces_AVX_None < data_size > op ((float *) &f, (float *) &u,
                                                  (float *) &N, (float *) &W,
                                                  (float *) &h,
                                                  (int *) &spring_id,
                                                  (int *) &spring_id_X,
                                                  (int *) &spring_id_Y,
                                                  (int *) &spring_id_Z,
                                                  (float **)
                                                  &extern_collision_attach,
                                                  (float **)
                                                  &extern_collision_stiffness);

      for (int t = threads; t <= threads_max;
           t += std::max < int >(((threads_max - threads) / 30), 1))
        {
          std::cout << "Running Test with " << t << " threads." << std::endl;
          MT_Streaming_Kernels::Kernel_Serial_Base_Helper <
            Collision_Forces_AVX_None < data_size > >helper (op, data_size, t);

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


      Collision_Forces_NEON_None < data_size > op ((float *) &f, (float *) &u,
                                                   (float *) &N, (float *) &W,
                                                   (float *) &h,
                                                   (int *) &spring_id,
                                                   (int *) &spring_id_X,
                                                   (int *) &spring_id_Y,
                                                   (int *) &spring_id_Z,
                                                   (float **)
                                                   &extern_collision_attach,
                                                   (float **)
                                                   &extern_collision_stiffness);

      for (int t = threads; t <= threads_max;
           t += std::max < int >(((threads_max - threads) / 30), 1))
        {
          std::cout << "Running Test with " << t << " threads." << std::endl;
          MT_Streaming_Kernels::Kernel_Serial_Base_Helper <
            Collision_Forces_NEON_None < data_size > >helper (op, data_size, t);

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


      Collision_Forces_MIC_None < data_size > op ((float *) &f, (float *) &u,
                                                  (float *) &N, (float *) &W,
                                                  (float *) &h,
                                                  (int *) &spring_id,
                                                  (int *) &spring_id_X,
                                                  (int *) &spring_id_Y,
                                                  (int *) &spring_id_Z,
                                                  (float **)
                                                  &extern_collision_attach,
                                                  (float **)
                                                  &extern_collision_stiffness);

      for (int t = threads; t <= threads_max;
           t += std::max < int >(((threads_max - threads) / 30), 1))
        {
          std::cout << "Running Test with " << t << " threads." << std::endl;
          MT_Streaming_Kernels::Kernel_Serial_Base_Helper <
            Collision_Forces_MIC_None < data_size > >helper (op, data_size, t);

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
    _mm_free (reinterpret_cast < void *>(u));
    _mm_free (reinterpret_cast < void *>(N));
    _mm_free (reinterpret_cast < void *>(W));
    _mm_free (reinterpret_cast < void *>(h));
    _mm_free (reinterpret_cast < void *>(spring_id));
    _mm_free (reinterpret_cast < void *>(spring_id_X));
    _mm_free (reinterpret_cast < void *>(spring_id_Y));
    _mm_free (reinterpret_cast < void *>(spring_id_Z));
    _mm_free (reinterpret_cast < void *>(extern_collision_attach));
    _mm_free (reinterpret_cast < void *>(extern_collision_stiffness));


  }


  return 0;

}

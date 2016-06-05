
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <sys/time.h>
#include "KernelCommon.h"
struct NEOHOOKEAN_TAG;
struct COROTATED_TAG;
#include <omp.h>
#include "Add_Force_Differential.h"
//#include "Force_Stabilization.h"

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

float du[3][8][16] __attribute__ ((aligned (64)));
float dp[16] __attribute__ ((aligned (64)));
float alpha_squared_over_kappa[16] __attribute__ ((aligned (64)));
float alpha[16] __attribute__ ((aligned (64)));
float one_over_h[16] __attribute__ ((aligned (64)));
float cell_volume[16] __attribute__ ((aligned (64)));
float Q_hat[3][16] __attribute__ ((aligned (64)));
float U[9][16] __attribute__ ((aligned (64)));
float V[9][16] __attribute__ ((aligned (64)));
float dPdF[12][16] __attribute__ ((aligned (64)));
float stab_constant[16] __attribute__ ((aligned (64)));
float df[3][8][16] __attribute__ ((aligned (64)));
float df_reference[3][8][16] __attribute__ ((aligned (64)));
float df_original[3][8][16] __attribute__ ((aligned (64)));
float dq[16] __attribute__ ((aligned (64)));
float dq_reference[16] __attribute__ ((aligned (64)));
float dq_original[16] __attribute__ ((aligned (64)));

#pragma omp threadprivate(du, dp, alpha_squared_over_kappa, alpha, one_over_h, cell_volume, Q_hat, U, V, dPdF, stab_constant, df, df_reference, df_original, dq, dq_reference, dq_original)

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

  omp_set_dynamic(0);


  {
    std::cout << "Running Stream Test for Add_Force_Differential " << std::endl;

//=======================================================
//
//        DEFINE ALL VARIABLES USED BY KERNEL
//
//=======================================================
 
    for (int __a = 0; __a < 3; __a++)
      for (int __b = 0; __b < 8; __b++)
        for (int __c = 0; __c < 16; __c++)
          du[__a][__b][__c] = Get_Random < float >();
    for (int __a = 0; __a < 16; __a++)
      dp[__a] = Get_Random < float >();
    for (int __a = 0; __a < 16; __a++)
      alpha_squared_over_kappa[__a] = Get_Random < float >();
    for (int __a = 0; __a < 16; __a++)
      alpha[__a] = Get_Random < float >();
    for (int __a = 0; __a < 16; __a++)
      one_over_h[__a] = Get_Random < float >();
    for (int __a = 0; __a < 16; __a++)
      cell_volume[__a] = Get_Random < float >();
    for (int __a = 0; __a < 3; __a++)
      for (int __b = 0; __b < 16; __b++)
        Q_hat[__a][__b] = Get_Random < float >();
    for (int __a = 0; __a < 9; __a++)
      for (int __b = 0; __b < 16; __b++)
        U[__a][__b] = Get_Random < float >();
    for (int __a = 0; __a < 9; __a++)
      for (int __b = 0; __b < 16; __b++)
        V[__a][__b] = Get_Random < float >();
    for (int __a = 0; __a < 12; __a++)
      for (int __b = 0; __b < 16; __b++)
        dPdF[__a][__b] = Get_Random < float >();
    for (int __a = 0; __a < 3; __a++)
      for (int __b = 0; __b < 8; __b++)
        for (int __c = 0; __c < 16; __c++)
          {
            df_original[__a][__b][__c] = Get_Random < float >();
            df[__a][__b][__c] = df_original[__a][__b][__c];
            df_reference[__a][__b][__c] = df_original[__a][__b][__c];
          }
    for (int __a = 0; __a < 16; __a++)
      {
        dq_original[__a] = Get_Random < float >();
        dq[__a] = dq_original[__a];
        dq_reference[__a] = dq_original[__a];
      }
    for (int __a = 0; __a < 16; __a++)
        stab_constant[__a] = Get_Random < float >();


//=======================================================
//
//             COMPUTE SCALAR RESULTS
//
//=======================================================

    {
      std::cout << "	Running " << NUM_TRIALS << " of SCALAR :  ";
      for( int k= 0; k < 20; k++){
      start_timer ();
#pragma omp parallel for copyin(df,dq)
      for (int n = 0; n < NUM_TRIALS; n++)
        {
          typedef T (&refArray1)[3][8][16];
          typedef T (&refArray2)[16];
          typedef T (&refArray3)[16];
          typedef T (&refArray4)[16];
          typedef T (&refArray5)[16];
          typedef T (&refArray6)[16];
          typedef T (&refArray7)[3][16];
          typedef T (&refArray8)[9][16];
          typedef T (&refArray9)[9][16];
          typedef T (&refArray10)[12][16];
          typedef T (&refArray11)[3][8][16];
          typedef T (&refArray12)[16];
          typedef T (&refArray13)[16];
          for (int i = 0; i < 16; i += 1)
            {
              refArray1 duk = reinterpret_cast < refArray1 > (du[0][0][i]);
              refArray2 dpk = reinterpret_cast < refArray2 > (dp[i]);
              refArray3 alpha_squared_over_kappak =
                reinterpret_cast < refArray3 > (alpha_squared_over_kappa[i]);
              refArray4 alphak = reinterpret_cast < refArray4 > (alpha[i]);
              refArray5 one_over_hk =
                reinterpret_cast < refArray5 > (one_over_h[i]);
              refArray6 cell_volumek =
                reinterpret_cast < refArray6 > (cell_volume[i]);
              refArray7 Q_hatk = reinterpret_cast < refArray7 > (Q_hat[0][i]);
              refArray8 Uk = reinterpret_cast < refArray8 > (U[0][i]);
              refArray9 Vk = reinterpret_cast < refArray9 > (V[0][i]);
              refArray10 dPdFk = reinterpret_cast < refArray10 > (dPdF[0][i]);
              refArray11 dfk = reinterpret_cast < refArray11 > (df[0][0][i]);
              refArray12 dqk = reinterpret_cast < refArray12 > (dq[i]);
              refArray13 stab_constantk = reinterpret_cast < refArray13 > (stab_constant[i]);
#if 1
              Add_Force_Differential < float, float[16], int[16] > (duk, dpk,
                                                                    alpha_squared_over_kappak,
                                                                    alphak,
                                                                    one_over_hk,
                                                                    cell_volumek,
                                                                    Q_hatk, Uk,
                                                                    Vk, dPdFk,
                                                                    dfk, dqk);
#endif
#if 0
              Add_Force_Differential < float, float[16], int[16] > (duk, dpk,
                                                                    alpha_squared_over_kappak,
                                                                    alphak,
                                                                    one_over_hk,
                                                                    cell_volumek,
                                                                    Q_hatk, Uk,
                                                                    Vk, dPdFk,
                                                                    dfk, dqk);
#endif
//              Force_Stabilization < float, float[16], int[16] > (duk, stab_constantk, dfk);

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
#pragma omp parallel for copyin(df,dq)
      for (int n = 0; n < NUM_TRIALS; n++)
        {
          typedef T (&refArray1)[3][8][16];
          typedef T (&refArray2)[16];
          typedef T (&refArray3)[16];
          typedef T (&refArray4)[16];
          typedef T (&refArray5)[16];
          typedef T (&refArray6)[16];
          typedef T (&refArray7)[3][16];
          typedef T (&refArray8)[9][16];
          typedef T (&refArray9)[9][16];
          typedef T (&refArray10)[12][16];
          typedef T (&refArray11)[3][8][16];
          typedef T (&refArray12)[16];
          for (int i = 0; i < 16; i += 4)
            {
              refArray1 duk = reinterpret_cast < refArray1 > (du[0][0][i]);
              refArray2 dpk = reinterpret_cast < refArray2 > (dp[i]);
              refArray3 alpha_squared_over_kappak =
                reinterpret_cast < refArray3 > (alpha_squared_over_kappa[i]);
              refArray4 alphak = reinterpret_cast < refArray4 > (alpha[i]);
              refArray5 one_over_hk =
                reinterpret_cast < refArray5 > (one_over_h[i]);
              refArray6 cell_volumek =
                reinterpret_cast < refArray6 > (cell_volume[i]);
              refArray7 Q_hatk = reinterpret_cast < refArray7 > (Q_hat[0][i]);
              refArray8 Uk = reinterpret_cast < refArray8 > (U[0][i]);
              refArray9 Vk = reinterpret_cast < refArray9 > (V[0][i]);
              refArray10 dPdFk = reinterpret_cast < refArray10 > (dPdF[0][i]);
              refArray11 dfk = reinterpret_cast < refArray11 > (df[0][0][i]);
              refArray12 dqk = reinterpret_cast < refArray12 > (dq[i]);
              Add_Force_Differential < __m128, float[16], int[16] > (duk, dpk,
                                                                     alpha_squared_over_kappak,
                                                                     alphak,
                                                                     one_over_hk,
                                                                     cell_volumek,
                                                                     Q_hatk, Uk,
                                                                     Vk, dPdFk,
                                                                     dfk, dqk);
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
      for( int k= 0; k < 20; k++){
      start_timer ();
#pragma omp parallel for copyin(df,dq)
      for (int n = 0; n < NUM_TRIALS; n++)
        {
          typedef T (&refArray1)[3][8][16];
          typedef T (&refArray2)[16];
          typedef T (&refArray3)[16];
          typedef T (&refArray4)[16];
          typedef T (&refArray5)[16];
          typedef T (&refArray6)[16];
          typedef T (&refArray7)[3][16];
          typedef T (&refArray8)[9][16];
          typedef T (&refArray9)[9][16];
          typedef T (&refArray10)[12][16];
          typedef T (&refArray11)[3][8][16];
          typedef T (&refArray12)[16];
          for (int i = 0; i < 16; i += 8)
            {
              refArray1 duk = reinterpret_cast < refArray1 > (du[0][0][i]);
              refArray2 dpk = reinterpret_cast < refArray2 > (dp[i]);
              refArray3 alpha_squared_over_kappak =
                reinterpret_cast < refArray3 > (alpha_squared_over_kappa[i]);
              refArray4 alphak = reinterpret_cast < refArray4 > (alpha[i]);
              refArray5 one_over_hk =
                reinterpret_cast < refArray5 > (one_over_h[i]);
              refArray6 cell_volumek =
                reinterpret_cast < refArray6 > (cell_volume[i]);
              refArray7 Q_hatk = reinterpret_cast < refArray7 > (Q_hat[0][i]);
              refArray8 Uk = reinterpret_cast < refArray8 > (U[0][i]);
              refArray9 Vk = reinterpret_cast < refArray9 > (V[0][i]);
              refArray10 dPdFk = reinterpret_cast < refArray10 > (dPdF[0][i]);
              refArray11 dfk = reinterpret_cast < refArray11 > (df[0][0][i]);
              refArray12 dqk = reinterpret_cast < refArray12 > (dq[i]);
              Add_Force_Differential < __m256, float[16], int[16] > (duk, dpk,
                                                                     alpha_squared_over_kappak,
                                                                     alphak,
                                                                     one_over_hk,
                                                                     cell_volumek,
                                                                     Q_hatk, Uk,
                                                                     Vk, dPdFk,
                                                                     dfk, dqk);
#if 0
              Add_Force_Differential < __m256, float[16], int[16] > (duk, dpk,
                                                                     alpha_squared_over_kappak,
                                                                     alphak,
                                                                     one_over_hk,
                                                                     cell_volumek,
                                                                     Q_hatk, Uk,
                                                                     Vk, dPdFk,
                                                                     dfk, dqk);
#endif
            }

        }
      stop_timer ();
      std::cout << get_time () << "s" << std::endl;
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
      std::cout << "	Running " << NUM_TRIALS << " of NEON :  ";
      start_timer ();
#pragma omp parallel for copyin(df,dq)
      for (int n = 0; n < NUM_TRIALS; n++)
        {
          typedef T (&refArray1)[3][8][16];
          typedef T (&refArray2)[16];
          typedef T (&refArray3)[16];
          typedef T (&refArray4)[16];
          typedef T (&refArray5)[16];
          typedef T (&refArray6)[16];
          typedef T (&refArray7)[3][16];
          typedef T (&refArray8)[9][16];
          typedef T (&refArray9)[9][16];
          typedef T (&refArray10)[12][16];
          typedef T (&refArray11)[3][8][16];
          typedef T (&refArray12)[16];
          for (int i = 0; i < 16; i += 4)
            {
              refArray1 duk = reinterpret_cast < refArray1 > (du[0][0][i]);
              refArray2 dpk = reinterpret_cast < refArray2 > (dp[i]);
              refArray3 alpha_squared_over_kappak =
                reinterpret_cast < refArray3 > (alpha_squared_over_kappa[i]);
              refArray4 alphak = reinterpret_cast < refArray4 > (alpha[i]);
              refArray5 one_over_hk =
                reinterpret_cast < refArray5 > (one_over_h[i]);
              refArray6 cell_volumek =
                reinterpret_cast < refArray6 > (cell_volume[i]);
              refArray7 Q_hatk = reinterpret_cast < refArray7 > (Q_hat[0][i]);
              refArray8 Uk = reinterpret_cast < refArray8 > (U[0][i]);
              refArray9 Vk = reinterpret_cast < refArray9 > (V[0][i]);
              refArray10 dPdFk = reinterpret_cast < refArray10 > (dPdF[0][i]);
              refArray11 dfk = reinterpret_cast < refArray11 > (df[0][0][i]);
              refArray12 dqk = reinterpret_cast < refArray12 > (dq[i]);
              Add_Force_Differential < float32x4_t, float[16], int[16] > (duk,
                                                                          dpk,
                                                                          alpha_squared_over_kappak,
                                                                          alphak,
                                                                          one_over_hk,
                                                                          cell_volumek,
                                                                          Q_hatk,
                                                                          Uk,
                                                                          Vk,
                                                                          dPdFk,
                                                                          dfk,
                                                                          dqk);
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
#pragma omp parallel for copyin(df,dq)
      for (int n = 0; n < NUM_TRIALS; n++)
        {
          typedef T (&refArray1)[3][8][16];
          typedef T (&refArray2)[16];
          typedef T (&refArray3)[16];
          typedef T (&refArray4)[16];
          typedef T (&refArray5)[16];
          typedef T (&refArray6)[16];
          typedef T (&refArray7)[3][16];
          typedef T (&refArray8)[9][16];
          typedef T (&refArray9)[9][16];
          typedef T (&refArray10)[12][16];
          typedef T (&refArray11)[3][8][16];
          typedef T (&refArray12)[16];
          typedef T (&refArray13)[16];
          for (int i = 0; i < 16; i += 16)
            {
              refArray1 duk = reinterpret_cast < refArray1 > (du[0][0][i]);
              refArray2 dpk = reinterpret_cast < refArray2 > (dp[i]);
              refArray3 alpha_squared_over_kappak =
                reinterpret_cast < refArray3 > (alpha_squared_over_kappa[i]);
              refArray4 alphak = reinterpret_cast < refArray4 > (alpha[i]);
              refArray5 one_over_hk =
                reinterpret_cast < refArray5 > (one_over_h[i]);
              refArray6 cell_volumek =
                reinterpret_cast < refArray6 > (cell_volume[i]);
              refArray7 Q_hatk = reinterpret_cast < refArray7 > (Q_hat[0][i]);
              refArray8 Uk = reinterpret_cast < refArray8 > (U[0][i]);
              refArray9 Vk = reinterpret_cast < refArray9 > (V[0][i]);
              refArray10 dPdFk = reinterpret_cast < refArray10 > (dPdF[0][i]);
              refArray11 dfk = reinterpret_cast < refArray11 > (df[0][0][i]);
              refArray12 dqk = reinterpret_cast < refArray12 > (dq[i]);
              refArray13 stab_constantk = reinterpret_cast < refArray13 > (stab_constant[i]);
#if 1
              Add_Force_Differential < __m512, float[16], int[16] > (duk, dpk,
                                                                     alpha_squared_over_kappak,
                                                                     alphak,
                                                                     one_over_hk,
                                                                     cell_volumek,
                                                                     Q_hatk, Uk,
                                                                     Vk, dPdFk,
                                                                     dfk, dqk);
#endif
#if 0
              Add_Force_Differential < __m512, float[16], int[16] > (duk, dpk,
                                                                     alpha_squared_over_kappak,
                                                                     alphak,
                                                                     one_over_hk,
                                                                     cell_volumek,
                                                                     Q_hatk, Uk,
                                                                     Vk, dPdFk,
                                                                     dfk, dqk);
#endif
              //Force_Stabilization < __m512, float[16], int[16] > (duk, stab_constantk, dfk);
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

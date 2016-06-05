
#include <cstdlib>
#include <iostream>
#include "KernelCommon.h"

struct NEOHOOKEAN_TAG;
struct COROTATED_TAG;

#include "Add_Force_Differential.h"
#include "Add_Force_Differential_Reference.h"

template < class T > T Get_Random (const T a = (T) - 1., const T b = (T) 1.)
{
  return ((b - a) * (T) rand ()) / (T) RAND_MAX + a;
}

int
main (int argc, char *argv[])
{
  typedef float   T;

  int             seed = 1;
  if (argc == 2)
    seed = atoi (argv[1]);
  srand (seed);

  std::cout.precision (10);
  std::cout.setf (std::ios::fixed, std::ios::floatfield);



  {
    std::cout << "Running SIMD Test for Add_Force_Differential " << std::endl;


//=======================================================
//
//        DEFINE ALL VARIABLES USED BY KERNEL
//
//=======================================================

    T               du[3][8][16] __attribute__ ((aligned (64)));
    T               dp[16] __attribute__ ((aligned (64)));
    T               alpha_squared_over_kappa[16] __attribute__ ((aligned (64)));
    T               alpha[16] __attribute__ ((aligned (64)));
    T               one_over_h[16] __attribute__ ((aligned (64)));
    T               cell_volume[16] __attribute__ ((aligned (64)));
    T               Q_hat[3][16] __attribute__ ((aligned (64)));
    T               U[9][16] __attribute__ ((aligned (64)));
    T               V[9][16] __attribute__ ((aligned (64)));
    T               dPdF[12][16] __attribute__ ((aligned (64)));
    T               df[3][8][16] __attribute__ ((aligned (64)));
    T               df_reference[3][8][16] __attribute__ ((aligned (64)));
    T               df_original[3][8][16] __attribute__ ((aligned (64)));
    T               dq[16] __attribute__ ((aligned (64)));
    T               dq_reference[16] __attribute__ ((aligned (64)));
    T               dq_original[16] __attribute__ ((aligned (64)));


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


//=======================================================
//
//             COMPUTE REFERENCE RESULTS
//
//=======================================================

    T               __mdu[3][8] __attribute__ ((aligned (4)));
    T __mdp __attribute__ ((aligned (4)));
    T __malpha_squared_over_kappa __attribute__ ((aligned (4)));
    T __malpha __attribute__ ((aligned (4)));
    T __mone_over_h __attribute__ ((aligned (4)));
    T __mcell_volume __attribute__ ((aligned (4)));
    T               __mQ_hat[3] __attribute__ ((aligned (4)));
    T               __mU[9] __attribute__ ((aligned (4)));
    T               __mV[9] __attribute__ ((aligned (4)));
    T               __mdPdF[12] __attribute__ ((aligned (4)));
    T               __mdf[3][8] __attribute__ ((aligned (4)));
    T               __mdf_reference[3][8] __attribute__ ((aligned (4)));
    T               __mdf_original[3][8] __attribute__ ((aligned (4)));
    T __mdq __attribute__ ((aligned (4)));
    T __mdq_reference __attribute__ ((aligned (4)));
    T __mdq_original __attribute__ ((aligned (4)));
    for (int k = 0; k < 16; k++)
    {
      for (int __a = 0; __a < 3; __a++)
        for (int __b = 0; __b < 8; __b++)
          __mdu[__a][__b] = du[__a][__b][k];
      __mdp = dp[k];
      __malpha_squared_over_kappa = alpha_squared_over_kappa[k];
      __malpha = alpha[k];
      __mone_over_h = one_over_h[k];
      __mcell_volume = cell_volume[k];
      for (int __a = 0; __a < 3; __a++)
        __mQ_hat[__a] = Q_hat[__a][k];
      for (int __a = 0; __a < 9; __a++)
        __mU[__a] = U[__a][k];
      for (int __a = 0; __a < 9; __a++)
        __mV[__a] = V[__a][k];
      for (int __a = 0; __a < 12; __a++)
        __mdPdF[__a] = dPdF[__a][k];
      for (int __a = 0; __a < 3; __a++)
        for (int __b = 0; __b < 8; __b++)
          __mdf_reference[__a][__b] = df_reference[__a][__b][k];

      __mdq_reference = dq_reference[k];
      Add_Force_Differential < float, float, int >(__mdu, __mdp,
                                                   __malpha_squared_over_kappa,
                                                   __malpha, __mone_over_h,
                                                   __mcell_volume, __mQ_hat,
                                                   __mU, __mV, __mdPdF,
                                                   __mdf_reference,
                                                   __mdq_reference);
      for (int __a = 0; __a < 3; __a++)
        for (int __b = 0; __b < 8; __b++)
          df_reference[__a][__b][k] = __mdf_reference[__a][__b];
      dq_reference[k] = __mdq_reference;
    }

//=======================================================
//
//               COMPUTE SCALAR RESULTS
//
//=======================================================

    {
      typedef         T (&refArray1)[3][8][16];
      typedef         T (&refArray2)[16];
      typedef         T (&refArray3)[16];
      typedef         T (&refArray4)[16];
      typedef         T (&refArray5)[16];
      typedef         T (&refArray6)[16];
      typedef         T (&refArray7)[3][16];
      typedef         T (&refArray8)[9][16];
      typedef         T (&refArray9)[9][16];
      typedef         T (&refArray10)[12][16];
      typedef         T (&refArray11)[3][8][16];
      typedef         T (&refArray12)[16];
      for (int __a = 0; __a < 3; __a++)
        for (int __b = 0; __b < 8; __b++)
          for (int __c = 0; __c < 16; __c++)
            df[__a][__b][__c] = df_original[__a][__b][__c];
      for (int __a = 0; __a < 16; __a++)
        dq[__a] = dq_original[__a];
      for (int i = 0; i < 16; i += 1)
      {
        refArray1       duk = reinterpret_cast < refArray1 > (du[0][0][i]);
        refArray2       dpk = reinterpret_cast < refArray2 > (dp[i]);
        refArray3       alpha_squared_over_kappak =
          reinterpret_cast < refArray3 > (alpha_squared_over_kappa[i]);
        refArray4       alphak = reinterpret_cast < refArray4 > (alpha[i]);
        refArray5       one_over_hk =
          reinterpret_cast < refArray5 > (one_over_h[i]);
        refArray6       cell_volumek =
          reinterpret_cast < refArray6 > (cell_volume[i]);
        refArray7       Q_hatk = reinterpret_cast < refArray7 > (Q_hat[0][i]);
        refArray8       Uk = reinterpret_cast < refArray8 > (U[0][i]);
        refArray9       Vk = reinterpret_cast < refArray9 > (V[0][i]);
        refArray10      dPdFk = reinterpret_cast < refArray10 > (dPdF[0][i]);
        refArray11      dfk = reinterpret_cast < refArray11 > (df[0][0][i]);
        refArray12      dqk = reinterpret_cast < refArray12 > (dq[i]);
        Add_Force_Differential < float, float[16], int[16] > (duk, dpk,
                                                              alpha_squared_over_kappak,
                                                              alphak,
                                                              one_over_hk,
                                                              cell_volumek,
                                                              Q_hatk, Uk, Vk,
                                                              dPdFk, dfk, dqk);
      }
      for (int __a = 0; __a < 3; __a++)
        for (int __b = 0; __b < 8; __b++)
          for (int __c = 0; __c < 16; __c++)
            if (std::
                abs ((df[__a][__b][__c] -
                      df_reference[__a][__b][__c]) /
                     (df_reference[__a][__b][__c])) > 1)
            {
              std::cerr << "Mismatch detected in SCALAR implementation" << std::
                endl;
              std::cerr << "Variable df:" << std::endl;
              std::
                cerr << "seed=" << seed << ", __a=" << __a << ", __b=" << __b <<
                ", __c=" << __c << std::endl;
              std::cerr << "df SCALAR=  " << df[__a][__b][__c] << std::endl;
              std::
                cerr << "df Reference=  " << df_reference[__a][__b][__c] <<
                std::endl;
              std::cerr << "df Rel Difference=  " << std::
                abs ((df[__a][__b][__c] -
                      df_reference[__a][__b][__c]) /
                     (df_reference[__a][__b][__c])) << std::endl;
              std::cerr << "df Abs Difference=  " << std::
                abs (df[__a][__b][__c] -
                     df_reference[__a][__b][__c]) << std::endl;
              return 1;
            }
      for (int __a = 0; __a < 16; __a++)
        if (std::abs ((dq[__a] - dq_reference[__a]) / (dq_reference[__a])) > 1)
        {
          std::cerr << "Mismatch detected in SCALAR implementation" << std::
            endl;
          std::cerr << "Variable dq:" << std::endl;
          std::cerr << "seed=" << seed << ", __a=" << __a << std::endl;
          std::cerr << "dq SCALAR=  " << dq[__a] << std::endl;
          std::cerr << "dq Reference=  " << dq_reference[__a] << std::endl;
          std::cerr << "dq Rel Difference=  " << std::
            abs ((dq[__a] -
                  dq_reference[__a]) / (dq_reference[__a])) << std::endl;
          std::cerr << "dq Abs Difference=  " << std::abs (dq[__a] -
                                                           dq_reference[__a]) <<
            std::endl;
          return 1;
        }

    }

//=======================================================
//
//               COMPUTE SSE RESULTS
//
//=======================================================

#ifdef ENABLE_SSE_INSTRUCTION_SET
    {
      typedef         T (&refArray1)[3][8][16];
      typedef         T (&refArray2)[16];
      typedef         T (&refArray3)[16];
      typedef         T (&refArray4)[16];
      typedef         T (&refArray5)[16];
      typedef         T (&refArray6)[16];
      typedef         T (&refArray7)[3][16];
      typedef         T (&refArray8)[9][16];
      typedef         T (&refArray9)[9][16];
      typedef         T (&refArray10)[12][16];
      typedef         T (&refArray11)[3][8][16];
      typedef         T (&refArray12)[16];
      for (int __a = 0; __a < 3; __a++)
        for (int __b = 0; __b < 8; __b++)
          for (int __c = 0; __c < 16; __c++)
            df[__a][__b][__c] = df_original[__a][__b][__c];
      for (int __a = 0; __a < 16; __a++)
        dq[__a] = dq_original[__a];
      for (int i = 0; i < 16; i += 4)
      {
        refArray1       duk = reinterpret_cast < refArray1 > (du[0][0][i]);
        refArray2       dpk = reinterpret_cast < refArray2 > (dp[i]);
        refArray3       alpha_squared_over_kappak =
          reinterpret_cast < refArray3 > (alpha_squared_over_kappa[i]);
        refArray4       alphak = reinterpret_cast < refArray4 > (alpha[i]);
        refArray5       one_over_hk =
          reinterpret_cast < refArray5 > (one_over_h[i]);
        refArray6       cell_volumek =
          reinterpret_cast < refArray6 > (cell_volume[i]);
        refArray7       Q_hatk = reinterpret_cast < refArray7 > (Q_hat[0][i]);
        refArray8       Uk = reinterpret_cast < refArray8 > (U[0][i]);
        refArray9       Vk = reinterpret_cast < refArray9 > (V[0][i]);
        refArray10      dPdFk = reinterpret_cast < refArray10 > (dPdF[0][i]);
        refArray11      dfk = reinterpret_cast < refArray11 > (df[0][0][i]);
        refArray12      dqk = reinterpret_cast < refArray12 > (dq[i]);
        Add_Force_Differential < __m128, float[16], int[16] > (duk, dpk,
                                                               alpha_squared_over_kappak,
                                                               alphak,
                                                               one_over_hk,
                                                               cell_volumek,
                                                               Q_hatk, Uk, Vk,
                                                               dPdFk, dfk, dqk);
      }
      for (int __a = 0; __a < 3; __a++)
        for (int __b = 0; __b < 8; __b++)
          for (int __c = 0; __c < 16; __c++)
            if (std::
                abs ((df[__a][__b][__c] -
                      df_reference[__a][__b][__c]) /
                     (df_reference[__a][__b][__c])) > 1)
            {
              std::cerr << "Mismatch detected in SSE implementation" << std::
                endl;
              std::cerr << "Variable df:" << std::endl;
              std::
                cerr << "seed=" << seed << ", __a=" << __a << ", __b=" << __b <<
                ", __c=" << __c << std::endl;
              std::cerr << "df SSE=  " << df[__a][__b][__c] << std::endl;
              std::
                cerr << "df Reference=  " << df_reference[__a][__b][__c] <<
                std::endl;
              std::cerr << "df Rel Difference=  " << std::
                abs ((df[__a][__b][__c] -
                      df_reference[__a][__b][__c]) /
                     (df_reference[__a][__b][__c])) << std::endl;
              std::cerr << "df Abs Difference=  " << std::
                abs (df[__a][__b][__c] -
                     df_reference[__a][__b][__c]) << std::endl;
              return 1;
            }
      for (int __a = 0; __a < 16; __a++)
        if (std::abs ((dq[__a] - dq_reference[__a]) / (dq_reference[__a])) > 1)
        {
          std::cerr << "Mismatch detected in SSE implementation" << std::endl;
          std::cerr << "Variable dq:" << std::endl;
          std::cerr << "seed=" << seed << ", __a=" << __a << std::endl;
          std::cerr << "dq SSE=  " << dq[__a] << std::endl;
          std::cerr << "dq Reference=  " << dq_reference[__a] << std::endl;
          std::cerr << "dq Rel Difference=  " << std::
            abs ((dq[__a] -
                  dq_reference[__a]) / (dq_reference[__a])) << std::endl;
          std::cerr << "dq Abs Difference=  " << std::abs (dq[__a] -
                                                           dq_reference[__a]) <<
            std::endl;
          return 1;
        }

    }
#endif

//=======================================================
//
//               COMPUTE AVX RESULTS
//
//=======================================================

#ifdef ENABLE_AVX_INSTRUCTION_SET
    {
      typedef         T (&refArray1)[3][8][16];
      typedef         T (&refArray2)[16];
      typedef         T (&refArray3)[16];
      typedef         T (&refArray4)[16];
      typedef         T (&refArray5)[16];
      typedef         T (&refArray6)[16];
      typedef         T (&refArray7)[3][16];
      typedef         T (&refArray8)[9][16];
      typedef         T (&refArray9)[9][16];
      typedef         T (&refArray10)[12][16];
      typedef         T (&refArray11)[3][8][16];
      typedef         T (&refArray12)[16];
      for (int __a = 0; __a < 3; __a++)
        for (int __b = 0; __b < 8; __b++)
          for (int __c = 0; __c < 16; __c++)
            df[__a][__b][__c] = df_original[__a][__b][__c];
      for (int __a = 0; __a < 16; __a++)
        dq[__a] = dq_original[__a];
      for (int i = 0; i < 16; i += 8)
      {
        refArray1       duk = reinterpret_cast < refArray1 > (du[0][0][i]);
        refArray2       dpk = reinterpret_cast < refArray2 > (dp[i]);
        refArray3       alpha_squared_over_kappak =
          reinterpret_cast < refArray3 > (alpha_squared_over_kappa[i]);
        refArray4       alphak = reinterpret_cast < refArray4 > (alpha[i]);
        refArray5       one_over_hk =
          reinterpret_cast < refArray5 > (one_over_h[i]);
        refArray6       cell_volumek =
          reinterpret_cast < refArray6 > (cell_volume[i]);
        refArray7       Q_hatk = reinterpret_cast < refArray7 > (Q_hat[0][i]);
        refArray8       Uk = reinterpret_cast < refArray8 > (U[0][i]);
        refArray9       Vk = reinterpret_cast < refArray9 > (V[0][i]);
        refArray10      dPdFk = reinterpret_cast < refArray10 > (dPdF[0][i]);
        refArray11      dfk = reinterpret_cast < refArray11 > (df[0][0][i]);
        refArray12      dqk = reinterpret_cast < refArray12 > (dq[i]);
        Add_Force_Differential < __m256, float[16], int[16] > (duk, dpk,
                                                               alpha_squared_over_kappak,
                                                               alphak,
                                                               one_over_hk,
                                                               cell_volumek,
                                                               Q_hatk, Uk, Vk,
                                                               dPdFk, dfk, dqk);
      }
      for (int __a = 0; __a < 3; __a++)
        for (int __b = 0; __b < 8; __b++)
          for (int __c = 0; __c < 16; __c++)
            if (std::
                abs ((df[__a][__b][__c] -
                      df_reference[__a][__b][__c]) /
                     (df_reference[__a][__b][__c])) > 1)
            {
              std::cerr << "Mismatch detected in AVX implementation" << std::
                endl;
              std::cerr << "Variable df:" << std::endl;
              std::
                cerr << "seed=" << seed << ", __a=" << __a << ", __b=" << __b <<
                ", __c=" << __c << std::endl;
              std::cerr << "df AVX=  " << df[__a][__b][__c] << std::endl;
              std::
                cerr << "df Reference=  " << df_reference[__a][__b][__c] <<
                std::endl;
              std::cerr << "df Rel Difference=  " << std::
                abs ((df[__a][__b][__c] -
                      df_reference[__a][__b][__c]) /
                     (df_reference[__a][__b][__c])) << std::endl;
              std::cerr << "df Abs Difference=  " << std::
                abs (df[__a][__b][__c] -
                     df_reference[__a][__b][__c]) << std::endl;
              return 1;
            }
      for (int __a = 0; __a < 16; __a++)
        if (std::abs ((dq[__a] - dq_reference[__a]) / (dq_reference[__a])) > 1)
        {
          std::cerr << "Mismatch detected in AVX implementation" << std::endl;
          std::cerr << "Variable dq:" << std::endl;
          std::cerr << "seed=" << seed << ", __a=" << __a << std::endl;
          std::cerr << "dq AVX=  " << dq[__a] << std::endl;
          std::cerr << "dq Reference=  " << dq_reference[__a] << std::endl;
          std::cerr << "dq Rel Difference=  " << std::
            abs ((dq[__a] -
                  dq_reference[__a]) / (dq_reference[__a])) << std::endl;
          std::cerr << "dq Abs Difference=  " << std::abs (dq[__a] -
                                                           dq_reference[__a]) <<
            std::endl;
          return 1;
        }

    }
#endif

//=======================================================
//
//               COMPUTE NEON RESULTS
//
//=======================================================

#ifdef ENABLE_NEON_INSTRUCTION_SET
    {
      typedef         T (&refArray1)[3][8][16];
      typedef         T (&refArray2)[16];
      typedef         T (&refArray3)[16];
      typedef         T (&refArray4)[16];
      typedef         T (&refArray5)[16];
      typedef         T (&refArray6)[16];
      typedef         T (&refArray7)[3][16];
      typedef         T (&refArray8)[9][16];
      typedef         T (&refArray9)[9][16];
      typedef         T (&refArray10)[12][16];
      typedef         T (&refArray11)[3][8][16];
      typedef         T (&refArray12)[16];
      for (int __a = 0; __a < 3; __a++)
        for (int __b = 0; __b < 8; __b++)
          for (int __c = 0; __c < 16; __c++)
            df[__a][__b][__c] = df_original[__a][__b][__c];
      for (int __a = 0; __a < 16; __a++)
        dq[__a] = dq_original[__a];
      for (int i = 0; i < 16; i += 4)
      {
        refArray1       duk = reinterpret_cast < refArray1 > (du[0][0][i]);
        refArray2       dpk = reinterpret_cast < refArray2 > (dp[i]);
        refArray3       alpha_squared_over_kappak =
          reinterpret_cast < refArray3 > (alpha_squared_over_kappa[i]);
        refArray4       alphak = reinterpret_cast < refArray4 > (alpha[i]);
        refArray5       one_over_hk =
          reinterpret_cast < refArray5 > (one_over_h[i]);
        refArray6       cell_volumek =
          reinterpret_cast < refArray6 > (cell_volume[i]);
        refArray7       Q_hatk = reinterpret_cast < refArray7 > (Q_hat[0][i]);
        refArray8       Uk = reinterpret_cast < refArray8 > (U[0][i]);
        refArray9       Vk = reinterpret_cast < refArray9 > (V[0][i]);
        refArray10      dPdFk = reinterpret_cast < refArray10 > (dPdF[0][i]);
        refArray11      dfk = reinterpret_cast < refArray11 > (df[0][0][i]);
        refArray12      dqk = reinterpret_cast < refArray12 > (dq[i]);
        Add_Force_Differential < float32x4_t, float[16], int[16] > (duk, dpk,
                                                                    alpha_squared_over_kappak,
                                                                    alphak,
                                                                    one_over_hk,
                                                                    cell_volumek,
                                                                    Q_hatk, Uk,
                                                                    Vk, dPdFk,
                                                                    dfk, dqk);
      }
      for (int __a = 0; __a < 3; __a++)
        for (int __b = 0; __b < 8; __b++)
          for (int __c = 0; __c < 16; __c++)
            if (std::
                abs ((df[__a][__b][__c] -
                      df_reference[__a][__b][__c]) /
                     (df_reference[__a][__b][__c])) > 1)
            {
              std::cerr << "Mismatch detected in NEON implementation" << std::
                endl;
              std::cerr << "Variable df:" << std::endl;
              std::
                cerr << "seed=" << seed << ", __a=" << __a << ", __b=" << __b <<
                ", __c=" << __c << std::endl;
              std::cerr << "df NEON=  " << df[__a][__b][__c] << std::endl;
              std::
                cerr << "df Reference=  " << df_reference[__a][__b][__c] <<
                std::endl;
              std::cerr << "df Rel Difference=  " << std::
                abs ((df[__a][__b][__c] -
                      df_reference[__a][__b][__c]) /
                     (df_reference[__a][__b][__c])) << std::endl;
              std::cerr << "df Abs Difference=  " << std::
                abs (df[__a][__b][__c] -
                     df_reference[__a][__b][__c]) << std::endl;
              return 1;
            }
      for (int __a = 0; __a < 16; __a++)
        if (std::abs ((dq[__a] - dq_reference[__a]) / (dq_reference[__a])) > 1)
        {
          std::cerr << "Mismatch detected in NEON implementation" << std::endl;
          std::cerr << "Variable dq:" << std::endl;
          std::cerr << "seed=" << seed << ", __a=" << __a << std::endl;
          std::cerr << "dq NEON=  " << dq[__a] << std::endl;
          std::cerr << "dq Reference=  " << dq_reference[__a] << std::endl;
          std::cerr << "dq Rel Difference=  " << std::
            abs ((dq[__a] -
                  dq_reference[__a]) / (dq_reference[__a])) << std::endl;
          std::cerr << "dq Abs Difference=  " << std::abs (dq[__a] -
                                                           dq_reference[__a]) <<
            std::endl;
          return 1;
        }

    }
#endif

//=======================================================
//
//               COMPUTE MIC RESULTS
//
//=======================================================

#ifdef ENABLE_MIC_INSTRUCTION_SET
    {
      typedef         T (&refArray1)[3][8][16];
      typedef         T (&refArray2)[16];
      typedef         T (&refArray3)[16];
      typedef         T (&refArray4)[16];
      typedef         T (&refArray5)[16];
      typedef         T (&refArray6)[16];
      typedef         T (&refArray7)[3][16];
      typedef         T (&refArray8)[9][16];
      typedef         T (&refArray9)[9][16];
      typedef         T (&refArray10)[12][16];
      typedef         T (&refArray11)[3][8][16];
      typedef         T (&refArray12)[16];
      for (int __a = 0; __a < 3; __a++)
        for (int __b = 0; __b < 8; __b++)
          for (int __c = 0; __c < 16; __c++)
            df[__a][__b][__c] = df_original[__a][__b][__c];
      for (int __a = 0; __a < 16; __a++)
        dq[__a] = dq_original[__a];
      for (int i = 0; i < 16; i += 16)
      {
        refArray1       duk = reinterpret_cast < refArray1 > (du[0][0][i]);
        refArray2       dpk = reinterpret_cast < refArray2 > (dp[i]);
        refArray3       alpha_squared_over_kappak =
          reinterpret_cast < refArray3 > (alpha_squared_over_kappa[i]);
        refArray4       alphak = reinterpret_cast < refArray4 > (alpha[i]);
        refArray5       one_over_hk =
          reinterpret_cast < refArray5 > (one_over_h[i]);
        refArray6       cell_volumek =
          reinterpret_cast < refArray6 > (cell_volume[i]);
        refArray7       Q_hatk = reinterpret_cast < refArray7 > (Q_hat[0][i]);
        refArray8       Uk = reinterpret_cast < refArray8 > (U[0][i]);
        refArray9       Vk = reinterpret_cast < refArray9 > (V[0][i]);
        refArray10      dPdFk = reinterpret_cast < refArray10 > (dPdF[0][i]);
        refArray11      dfk = reinterpret_cast < refArray11 > (df[0][0][i]);
        refArray12      dqk = reinterpret_cast < refArray12 > (dq[i]);
        Add_Force_Differential < __m512, float[16], int[16] > (duk, dpk,
                                                               alpha_squared_over_kappak,
                                                               alphak,
                                                               one_over_hk,
                                                               cell_volumek,
                                                               Q_hatk, Uk, Vk,
                                                               dPdFk, dfk, dqk);
      }
      for (int __a = 0; __a < 3; __a++)
        for (int __b = 0; __b < 8; __b++)
          for (int __c = 0; __c < 16; __c++)
            if (std::
                abs ((df[__a][__b][__c] -
                      df_reference[__a][__b][__c]) /
                     (df_reference[__a][__b][__c])) > 1)
            {
              std::cerr << "Mismatch detected in MIC implementation" << std::
                endl;
              std::cerr << "Variable df:" << std::endl;
              std::
                cerr << "seed=" << seed << ", __a=" << __a << ", __b=" << __b <<
                ", __c=" << __c << std::endl;
              std::cerr << "df MIC=  " << df[__a][__b][__c] << std::endl;
              std::
                cerr << "df Reference=  " << df_reference[__a][__b][__c] <<
                std::endl;
              std::cerr << "df Rel Difference=  " << std::
                abs ((df[__a][__b][__c] -
                      df_reference[__a][__b][__c]) /
                     (df_reference[__a][__b][__c])) << std::endl;
              std::cerr << "df Abs Difference=  " << std::
                abs (df[__a][__b][__c] -
                     df_reference[__a][__b][__c]) << std::endl;
              return 1;
            }
      for (int __a = 0; __a < 16; __a++)
        if (std::abs ((dq[__a] - dq_reference[__a]) / (dq_reference[__a])) > 1)
        {
          std::cerr << "Mismatch detected in MIC implementation" << std::endl;
          std::cerr << "Variable dq:" << std::endl;
          std::cerr << "seed=" << seed << ", __a=" << __a << std::endl;
          std::cerr << "dq MIC=  " << dq[__a] << std::endl;
          std::cerr << "dq Reference=  " << dq_reference[__a] << std::endl;
          std::cerr << "dq Rel Difference=  " << std::
            abs ((dq[__a] -
                  dq_reference[__a]) / (dq_reference[__a])) << std::endl;
          std::cerr << "dq Abs Difference=  " << std::abs (dq[__a] -
                                                           dq_reference[__a]) <<
            std::endl;
          return 1;
        }

    }
#endif

  }



  std::cout << "SIMD check successful!" << std::endl;

  return 0;

}

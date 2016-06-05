
#include <cstdlib>
#include <iostream>

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
  typedef float T;

  int seed = 1;
  if (argc == 2)
    seed = atoi (argv[1]);
  srand (seed);



  {
    T du[3][8] __attribute__ ((aligned (4)));
    T dp __attribute__ ((aligned (4)));
    T alpha_squared_over_kappa __attribute__ ((aligned (4)));
    T alpha __attribute__ ((aligned (4)));
    T one_over_h __attribute__ ((aligned (4)));
    T cell_volume __attribute__ ((aligned (4)));
    T Q_hat[3] __attribute__ ((aligned (4)));
    T U[9] __attribute__ ((aligned (4)));
    T V[9] __attribute__ ((aligned (4)));
    T dPdF[12] __attribute__ ((aligned (4)));
    T df[3][8] __attribute__ ((aligned (4)));
    T df_reference[3][8] __attribute__ ((aligned (4)));
    T df_original[3][8] __attribute__ ((aligned (4)));
    T dq __attribute__ ((aligned (4)));
    T dq_reference __attribute__ ((aligned (4)));
    T dq_original __attribute__ ((aligned (4)));


    for (int __a = 0; __a < 3; __a++)
      for (int __b = 0; __b < 8; __b++)
        du[__a][__b] = Get_Random < float >();
    dp = Get_Random < float >();
    alpha_squared_over_kappa = Get_Random < float >();
    alpha = Get_Random < float >();
    one_over_h = Get_Random < float >();
    cell_volume = Get_Random < float >();
    for (int __a = 0; __a < 3; __a++)
      Q_hat[__a] = Get_Random < float >();
    for (int __a = 0; __a < 9; __a++)
      U[__a] = Get_Random < float >();
    for (int __a = 0; __a < 9; __a++)
      V[__a] = Get_Random < float >();
    for (int __a = 0; __a < 12; __a++)
      dPdF[__a] = Get_Random < float >();
    for (int __a = 0; __a < 3; __a++)
      for (int __b = 0; __b < 8; __b++)
        {
          df_original[__a][__b] = Get_Random < float >();
          df[__a][__b] = df_original[__a][__b];
          df_reference[__a][__b] = df_original[__a][__b];
        }
    {
      dq_original = Get_Random < float >();
      dq = dq_original;
      dq_reference = dq_original;
    }


    for (int __a = 0; __a < 3; __a++)
      for (int __b = 0; __b < 8; __b++)
        df[__a][__b] = df_original[__a][__b];
    dq = dq_original;
    for (int i = 0; i < 1; i += 1)
      {
        Add_Force_Differential < float, float, int >(du, dp,
                                                     alpha_squared_over_kappa,
                                                     alpha, one_over_h,
                                                     cell_volume, Q_hat, U, V,
                                                     dPdF, df, dq);
      }

    Add_Force_Differential_Reference < float >(du, dp, alpha_squared_over_kappa,
                                               alpha, one_over_h, cell_volume,
                                               Q_hat, U, V, dPdF, df_reference,
                                               dq_reference);
    if (!
        (Add_Force_Differential_Compare <
         float >(df, dq, df_reference, dq_reference)))
      {
        std::
          cout << "Failed to confirm unit test for Add_Force_Differential " <<
          std::endl;
        return 1;
      }

  }



  return 0;

}

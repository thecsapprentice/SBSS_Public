
#include <cstdlib>
#include <iostream>

struct NEOHOOKEAN_TAG;
struct COROTATED_TAG;

#include "Stress_Tensor_Differential.h"
#include "Stress_Tensor_Differential_Reference.h"

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
    T dP_hat[9] __attribute__ ((aligned (4)));
    T dP_hat_reference[9] __attribute__ ((aligned (4)));
    T dP_hat_original[9] __attribute__ ((aligned (4)));
    T dPdF[12] __attribute__ ((aligned (4)));
    T dF_hat[9] __attribute__ ((aligned (4)));
    T Q_hat[3] __attribute__ ((aligned (4)));
    T dp __attribute__ ((aligned (4)));
    T alpha __attribute__ ((aligned (4)));


    for (int __a = 0; __a < 9; __a++)
      {
        dP_hat_original[__a] = Get_Random < float >();
        dP_hat[__a] = dP_hat_original[__a];
        dP_hat_reference[__a] = dP_hat_original[__a];
      }
    for (int __a = 0; __a < 12; __a++)
      dPdF[__a] = Get_Random < float >();
    for (int __a = 0; __a < 9; __a++)
      dF_hat[__a] = Get_Random < float >();
    for (int __a = 0; __a < 3; __a++)
      Q_hat[__a] = Get_Random < float >();
    dp = Get_Random < float >();
    alpha = Get_Random < float >();

    for (int __a = 0; __a < 9; __a++)
      dP_hat[__a] = dP_hat_original[__a];
    for (int i = 0; i < 1; i += 1)
      {
        Stress_Tensor_Differential < float, float, int >(dP_hat, dPdF, dF_hat,
                                                         Q_hat, dp, alpha);
      }

    Stress_Tensor_Differential_Reference < float >(dP_hat_reference, dPdF,
                                                   dF_hat, Q_hat, dp, alpha);
    if (!
        (Stress_Tensor_Differential_Compare <
         float >(dP_hat, dP_hat_reference)))
      {
        std::
          cout << "Failed to confirm unit test for Stress_Tensor_Differential "
          << std::endl;
        return 1;
      }

  }



  return 0;

}

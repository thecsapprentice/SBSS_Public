
#include <cstdlib>
#include <iostream>

struct NEOHOOKEAN_TAG;
struct COROTATED_TAG;

#include "Pressure_Force_Differential.h"
#include "Pressure_Force_Differential_Reference.h"

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
    T dq __attribute__ ((aligned (4)));
    T dq_reference __attribute__ ((aligned (4)));
    T dq_original __attribute__ ((aligned (4)));
    T Q_hat[3] __attribute__ ((aligned (4)));
    T dF_hat[9] __attribute__ ((aligned (4)));
    T dp __attribute__ ((aligned (4)));
    T alpha __attribute__ ((aligned (4)));
    T alpha_squared_over_kappa __attribute__ ((aligned (4)));


    {
      dq_original = Get_Random < float >();
      dq = dq_original;
      dq_reference = dq_original;
    }
    for (int __a = 0; __a < 3; __a++)
      Q_hat[__a] = Get_Random < float >();
    for (int __a = 0; __a < 9; __a++)
      dF_hat[__a] = Get_Random < float >();
    dp = Get_Random < float >();
    alpha = Get_Random < float >();
    alpha_squared_over_kappa = Get_Random < float >();

    dq = dq_original;
    for (int i = 0; i < 1; i += 1)
      {
        Pressure_Force_Differential < float, float, int >(dq, Q_hat, dF_hat, dp,
                                                          alpha,
                                                          alpha_squared_over_kappa);
      }

    Pressure_Force_Differential_Reference < float >(dq_reference, Q_hat, dF_hat,
                                                    dp, alpha,
                                                    alpha_squared_over_kappa);
    if (!(Pressure_Force_Differential_Compare < float >(dq, dq_reference)))
      {
        std::
          cout << "Failed to confirm unit test for Pressure_Force_Differential "
          << std::endl;
        return 1;
      }

  }



  return 0;

}

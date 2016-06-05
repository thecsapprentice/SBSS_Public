
#include <cstdlib>
#include <iostream>

struct NEOHOOKEAN_TAG;
struct COROTATED_TAG;

#include "Penalty_Measure_Gradient.h"
#include "Penalty_Measure_Gradient_Reference.h"

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
    T Sigma[3] __attribute__ ((aligned (4)));
    T Q_hat[3] __attribute__ ((aligned (4)));
    T Q_hat_reference[3] __attribute__ ((aligned (4)));
    T Q_hat_original[3] __attribute__ ((aligned (4)));


    for (int __a = 0; __a < 3; __a++)
      Sigma[__a] = Get_Random < float >();
    for (int __a = 0; __a < 3; __a++)
      {
        Q_hat_original[__a] = Get_Random < float >();
        Q_hat[__a] = Q_hat_original[__a];
        Q_hat_reference[__a] = Q_hat_original[__a];
      }


    for (int __a = 0; __a < 3; __a++)
      Q_hat[__a] = Q_hat_original[__a];
    for (int i = 0; i < 1; i += 1)
      {
        Penalty_Measure_Gradient < NEOHOOKEAN_TAG, float, float,
          int >::Run (Sigma, Q_hat);
      }

    Penalty_Measure_Gradient_Neohookean_Reference < float >(Sigma,
                                                            Q_hat_reference);
    if (!
        (Penalty_Measure_Gradient_Neohookean_Compare <
         float >(Q_hat, Q_hat_reference)))
      {
        std::
          cout <<
          "Failed to confirm unit test for Penalty_Measure_Gradient with material NEOHOOKEAN"
          << std::endl;
        return 1;
      }

  }



  {
    T Sigma[3] __attribute__ ((aligned (4)));
    T Q_hat[3] __attribute__ ((aligned (4)));
    T Q_hat_reference[3] __attribute__ ((aligned (4)));
    T Q_hat_original[3] __attribute__ ((aligned (4)));


    for (int __a = 0; __a < 3; __a++)
      Sigma[__a] = Get_Random < float >();
    for (int __a = 0; __a < 3; __a++)
      {
        Q_hat_original[__a] = Get_Random < float >();
        Q_hat[__a] = Q_hat_original[__a];
        Q_hat_reference[__a] = Q_hat_original[__a];
      }


    for (int __a = 0; __a < 3; __a++)
      Q_hat[__a] = Q_hat_original[__a];
    for (int i = 0; i < 1; i += 1)
      {
        Penalty_Measure_Gradient < COROTATED_TAG, float, float,
          int >::Run (Sigma, Q_hat);
      }

    Penalty_Measure_Gradient_Corotated_Reference < float >(Sigma,
                                                           Q_hat_reference);
    if (!
        (Penalty_Measure_Gradient_Corotated_Compare <
         float >(Q_hat, Q_hat_reference)))
      {
        std::
          cout <<
          "Failed to confirm unit test for Penalty_Measure_Gradient with material COROTATED"
          << std::endl;
        return 1;
      }

  }



  return 0;

}

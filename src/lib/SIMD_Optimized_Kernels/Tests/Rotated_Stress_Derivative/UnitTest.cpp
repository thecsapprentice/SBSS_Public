
#include <cstdlib>
#include <iostream>

struct NEOHOOKEAN_TAG;
struct COROTATED_TAG;

#include "Rotated_Stress_Derivative.h"
#include "Rotated_Stress_Derivative_Reference.h"

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
    T dPdF[12] __attribute__ ((aligned (4)));
    T dPdF_reference[12] __attribute__ ((aligned (4)));
    T dPdF_original[12] __attribute__ ((aligned (4)));
    T Sigma[3] __attribute__ ((aligned (4)));
    T mu __attribute__ ((aligned (4)));
    T kappa __attribute__ ((aligned (4)));


    for (int __a = 0; __a < 12; __a++)
      {
        dPdF_original[__a] = Get_Random < float >();
        dPdF[__a] = dPdF_original[__a];
        dPdF_reference[__a] = dPdF_original[__a];
      }
    for (int __a = 0; __a < 3; __a++)
        Sigma[__a] = Get_Random < float >(0,1);
    mu = Get_Random < float >(0,1);
    kappa = Get_Random < float >(0,1);

    for (int __a = 0; __a < 12; __a++)
      dPdF[__a] = dPdF_original[__a];
    for (int i = 0; i < 1; i += 1)
      {
        Rotated_Stress_Derivative < NEOHOOKEAN_TAG, float, float,
          int >::Run (dPdF, Sigma, mu, kappa);
      }

    Rotated_Stress_Derivative_Neohookean_Reference < float >(dPdF_reference,
                                                             Sigma, mu, kappa);
    if (!
        (Rotated_Stress_Derivative_Neohookean_Compare <
         float >(dPdF, dPdF_reference)))
      {
        std::
          cout <<
          "Failed to confirm unit test for Rotated_Stress_Derivative with material NEOHOOKEAN"
          << std::endl;
        return 1;
      }

  }



  {
    T dPdF[12] __attribute__ ((aligned (4)));
    T dPdF_reference[12] __attribute__ ((aligned (4)));
    T dPdF_original[12] __attribute__ ((aligned (4)));
    T Sigma[3] __attribute__ ((aligned (4)));
    T mu __attribute__ ((aligned (4)));
    T kappa __attribute__ ((aligned (4)));


    for (int __a = 0; __a < 12; __a++)
      {
        dPdF_original[__a] = Get_Random < float >();
        dPdF[__a] = dPdF_original[__a];
        dPdF_reference[__a] = dPdF_original[__a];
      }
    for (int __a = 0; __a < 3; __a++)
      Sigma[__a] = Get_Random < float >();
    mu = Get_Random < float >();
    kappa = Get_Random < float >();

    for (int __a = 0; __a < 12; __a++)
      dPdF[__a] = dPdF_original[__a];
    for (int i = 0; i < 1; i += 1)
      {
        Rotated_Stress_Derivative < COROTATED_TAG, float, float,
          int >::Run (dPdF, Sigma, mu, kappa);
      }

    Rotated_Stress_Derivative_Corotated_Reference < float >(dPdF_reference,
                                                            Sigma, mu, kappa);
    if (!
        (Rotated_Stress_Derivative_Corotated_Compare <
         float >(dPdF, dPdF_reference)))
      {
        std::
          cout <<
          "Failed to confirm unit test for Rotated_Stress_Derivative with material COROTATED"
          << std::endl;
        return 1;
      }

  }



  return 0;

}


#include <cstdlib>
#include <iostream>

struct NEOHOOKEAN_TAG;
struct COROTATED_TAG;

#include "Rotated_Stress_Derivative_Contraction.h"
#include "Rotated_Stress_Derivative_Contraction_Reference.h"

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
    T dF_Hat[9] __attribute__ ((aligned (4)));
    T dP_Hat[9] __attribute__ ((aligned (4)));
    T dP_Hat_reference[9] __attribute__ ((aligned (4)));
    T dP_Hat_original[9] __attribute__ ((aligned (4)));


    for (int __a = 0; __a < 12; __a++)
      dPdF[__a] = Get_Random < float >();
    for (int __a = 0; __a < 9; __a++)
      dF_Hat[__a] = Get_Random < float >();
    for (int __a = 0; __a < 9; __a++)
      {
        dP_Hat_original[__a] = Get_Random < float >();
        dP_Hat[__a] = dP_Hat_original[__a];
        dP_Hat_reference[__a] = dP_Hat_original[__a];
      }


    for (int __a = 0; __a < 9; __a++)
      dP_Hat[__a] = dP_Hat_original[__a];
    for (int i = 0; i < 1; i += 1)
      {
        Rotated_Stress_Derivative_Contraction < float, float, int >(dPdF,
                                                                    dF_Hat,
                                                                    dP_Hat);
      }

    Rotated_Stress_Derivative_Contraction_Reference < float >(dPdF, dF_Hat,
                                                              dP_Hat_reference);
    if (!
        (Rotated_Stress_Derivative_Contraction_Compare <
         float >(dP_Hat, dP_Hat_reference)))
      {
        std::
          cout <<
          "Failed to confirm unit test for Rotated_Stress_Derivative_Contraction "
          << std::endl;
        return 1;
      }

  }



  return 0;

}

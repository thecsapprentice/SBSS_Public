
#include <cstdlib>
#include <iostream>

struct NEOHOOKEAN_TAG;
struct COROTATED_TAG;

#include "Unweighted_Gradient.h"
#include "Unweighted_Gradient_Reference.h"

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
    T u[3][8] __attribute__ ((aligned (4)));
    T F[9] __attribute__ ((aligned (4)));
    T F_reference[9] __attribute__ ((aligned (4)));
    T F_original[9] __attribute__ ((aligned (4)));
    T one_over_h __attribute__ ((aligned (4)));


    for (int __a = 0; __a < 3; __a++)
      for (int __b = 0; __b < 8; __b++)
        u[__a][__b] = Get_Random < float >();
    for (int __a = 0; __a < 9; __a++)
      {
        F_original[__a] = Get_Random < float >();
        F[__a] = F_original[__a];
        F_reference[__a] = F_original[__a];
      }

    one_over_h = Get_Random < float >();

    for (int __a = 0; __a < 9; __a++)
      F[__a] = F_original[__a];
    for (int i = 0; i < 1; i += 1)
      {
        Unweighted_Gradient < float, float, int >(u, F, one_over_h);
      }

    Unweighted_Gradient_Reference < float >(u, F_reference, one_over_h);
    if (!(Unweighted_Gradient_Compare < float >(F, F_reference)))
      {
        std::
          cout << "Failed to confirm unit test for Unweighted_Gradient " <<
          std::endl;
        return 1;
      }

  }



  return 0;

}

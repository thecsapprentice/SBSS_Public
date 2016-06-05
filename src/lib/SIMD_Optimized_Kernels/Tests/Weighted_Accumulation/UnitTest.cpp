
#include <cstdlib>
#include <iostream>

struct NEOHOOKEAN_TAG;
struct COROTATED_TAG;

#include "Weighted_Accumulation.h"
#include "Weighted_Accumulation_Reference.h"

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
    T u_reference[3][8] __attribute__ ((aligned (4)));
    T u_original[3][8] __attribute__ ((aligned (4)));
    T F[9] __attribute__ ((aligned (4)));
    T W[3] __attribute__ ((aligned (4)));
    T one_over_h __attribute__ ((aligned (4)));
    T scale __attribute__ ((aligned (4)));


    for (int __a = 0; __a < 3; __a++)
      for (int __b = 0; __b < 8; __b++)
        {
          u_original[__a][__b] = Get_Random < float >();
          u[__a][__b] = u_original[__a][__b];
          u_reference[__a][__b] = u_original[__a][__b];
        }
    for (int __a = 0; __a < 9; __a++)
      F[__a] = Get_Random < float >();
    for (int __a = 0; __a < 3; __a++)
      W[__a] = Get_Random < float >();
    one_over_h = Get_Random < float >();
    scale = Get_Random < float >();

    for (int __a = 0; __a < 3; __a++)
      for (int __b = 0; __b < 8; __b++)
        u[__a][__b] = u_original[__a][__b];
    for (int i = 0; i < 1; i += 1)
      {
        Weighted_Accumulation < float, float, int >(u, F, W, one_over_h, scale);
      }

    Weighted_Accumulation_Reference < float >(u_reference, F, W, one_over_h,
                                              scale);
    if (!(Weighted_Accumulation_Compare < float >(u, u_reference)))
      {
        std::
          cout << "Failed to confirm unit test for Weighted_Accumulation " <<
          std::endl;
        return 1;
      }

  }



  return 0;

}

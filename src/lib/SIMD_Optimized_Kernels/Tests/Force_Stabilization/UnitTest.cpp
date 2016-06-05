
#include <cstdlib>
#include <iostream>

struct NEOHOOKEAN_TAG;
struct COROTATED_TAG;

#include "Force_Stabilization.h"
#include "Force_Stabilization_Reference.h"

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
    T Du[3][8] __attribute__ ((aligned (4)));
    T constant __attribute__ ((aligned (4)));
    T dH[3][8] __attribute__ ((aligned (4)));
    T dH_reference[3][8] __attribute__ ((aligned (4)));
    T dH_original[3][8] __attribute__ ((aligned (4)));


    for (int __a = 0; __a < 3; __a++)
      for (int __b = 0; __b < 8; __b++)
        Du[__a][__b] = Get_Random < float >();
    constant = Get_Random < float >();
    for (int __a = 0; __a < 3; __a++)
      for (int __b = 0; __b < 8; __b++)
        {
          dH_original[__a][__b] = Get_Random < float >();
          dH[__a][__b] = dH_original[__a][__b];
          dH_reference[__a][__b] = dH_original[__a][__b];
        }


    for (int __a = 0; __a < 3; __a++)
      for (int __b = 0; __b < 8; __b++)
        dH[__a][__b] = dH_original[__a][__b];
    for (int i = 0; i < 1; i += 1)
      {
        Force_Stabilization < float, float, int >(Du, constant, dH);
      }

    Force_Stabilization_Reference < float >(Du, constant, dH_reference);
    if (!(Force_Stabilization_Compare < float >(dH, dH_reference)))
      {
        std::
          cout << "Failed to confirm unit test for Force_Stabilization " <<
          std::endl;
        return 1;
      }

  }



  return 0;

}

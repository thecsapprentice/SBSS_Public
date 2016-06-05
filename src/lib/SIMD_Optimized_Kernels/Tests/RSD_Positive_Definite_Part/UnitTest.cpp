
#include <cstdlib>
#include <iostream>

struct NEOHOOKEAN_TAG;
struct COROTATED_TAG;

#include "RSD_Positive_Definite_Part.h"
#include "RSD_Positive_Definite_Part_Reference.h"

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
    T RSD[12] __attribute__ ((aligned (4)));
    T RSDpd[12] __attribute__ ((aligned (4)));
    T RSDpd_reference[12] __attribute__ ((aligned (4)));
    T RSDpd_original[12] __attribute__ ((aligned (4)));


    for (int __a = 0; __a < 12; __a++)
      RSD[__a] = Get_Random < float >();
    for (int __a = 0; __a < 12; __a++)
      {
        RSDpd_original[__a] = Get_Random < float >();
        RSDpd[__a] = RSDpd_original[__a];
        RSDpd_reference[__a] = RSDpd_original[__a];
      }


    for (int __a = 0; __a < 12; __a++)
      RSDpd[__a] = RSDpd_original[__a];
    for (int i = 0; i < 1; i += 1)
      {
        RSD_Positive_Definite_Part < float, float, int >(RSD, RSDpd);
      }

    RSD_Positive_Definite_Part_Reference < float >(RSD, RSDpd_reference);
    if (!(RSD_Positive_Definite_Part_Compare < float >(RSDpd, RSDpd_reference)))
      {
        std::
          cout << "Failed to confirm unit test for RSD_Positive_Definite_Part "
          << std::endl;
        return 1;
      }

  }



  return 0;

}

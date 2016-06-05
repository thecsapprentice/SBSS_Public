
#include <cstdlib>
#include <iostream>

struct NEOHOOKEAN_TAG;
struct COROTATED_TAG;

#include "Symmetric_Definite_Projection.h"
#include "Symmetric_Definite_Projection_Reference.h"

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
    T A[6] __attribute__ ((aligned (4)));
    T Apd[6] __attribute__ ((aligned (4)));
    T Apd_reference[6] __attribute__ ((aligned (4)));
    T Apd_original[6] __attribute__ ((aligned (4)));


    for (int __a = 0; __a < 6; __a++)
      A[__a] = Get_Random < float >();
    for (int __a = 0; __a < 6; __a++)
      {
        Apd_original[__a] = Get_Random < float >();
        Apd[__a] = Apd_original[__a];
        Apd_reference[__a] = Apd_original[__a];
      }


    for (int __a = 0; __a < 6; __a++)
      Apd[__a] = Apd_original[__a];
    for (int i = 0; i < 1; i += 1)
      {
        Symmetric_Definite_Projection < float, float, int >(A, Apd);
      }

    Symmetric_Definite_Projection_Reference < float >(A, Apd_reference);
    if (!(Symmetric_Definite_Projection_Compare < float >(Apd, Apd_reference)))
      {
        std::
          cout <<
          "Failed to confirm unit test for Symmetric_Definite_Projection " <<
          std::endl;
        return 1;
      }

  }



  return 0;

}

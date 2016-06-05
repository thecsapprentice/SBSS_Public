
#include <cstdlib>
#include <iostream>

struct NEOHOOKEAN_TAG;
struct COROTATED_TAG;

#include "Muscle_Differential.h"
#include "Muscle_Differential_Reference.h"

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
    T dP_fiber[9] __attribute__ ((aligned (4)));
    T dP_fiber_reference[9] __attribute__ ((aligned (4)));
    T dP_fiber_original[9] __attribute__ ((aligned (4)));
    T dF[9] __attribute__ ((aligned (4)));
    T fiber[3] __attribute__ ((aligned (4)));
    T Ffiber[3] __attribute__ ((aligned (4)));
    T c1 __attribute__ ((aligned (4)));
    T c2 __attribute__ ((aligned (4)));


    for (int __a = 0; __a < 9; __a++)
      {
        dP_fiber_original[__a] = Get_Random < float >();
        dP_fiber[__a] = dP_fiber_original[__a];
        dP_fiber_reference[__a] = dP_fiber_original[__a];
      }
    for (int __a = 0; __a < 9; __a++)
      dF[__a] = Get_Random < float >();
    for (int __a = 0; __a < 3; __a++)
      fiber[__a] = Get_Random < float >();
    for (int __a = 0; __a < 3; __a++)
      Ffiber[__a] = Get_Random < float >();
    c1 = Get_Random < float >();
    c2 = Get_Random < float >();

    for (int __a = 0; __a < 9; __a++)
      dP_fiber[__a] = dP_fiber_original[__a];
    for (int i = 0; i < 1; i += 1)
      {
        Muscle_Differential < float, float, int >(dP_fiber, dF, fiber, Ffiber,
                                                  c1, c2);
      }

    Muscle_Differential_Reference < float >(dP_fiber_reference, dF, fiber,
                                            Ffiber, c1, c2);
    if (!(Muscle_Differential_Compare < float >(dP_fiber, dP_fiber_reference)))
      {
        std::
          cout << "Failed to confirm unit test for Muscle_Differential " <<
          std::endl;
        return 1;
      }

  }



  return 0;

}

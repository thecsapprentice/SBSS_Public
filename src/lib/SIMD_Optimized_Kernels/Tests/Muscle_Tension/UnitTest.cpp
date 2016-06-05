
#include <cstdlib>
#include <iostream>

struct NEOHOOKEAN_TAG;
struct COROTATED_TAG;

#include "Muscle_Tension.h"
#include "Muscle_Tension_Reference.h"

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
    T tension __attribute__ ((aligned (4)));
    T tension_reference __attribute__ ((aligned (4)));
    T tension_original __attribute__ ((aligned (4)));
    T stretch __attribute__ ((aligned (4)));
    T activation __attribute__ ((aligned (4)));
    T density __attribute__ ((aligned (4)));
    T fiber_max_stress __attribute__ ((aligned (4)));


    {
      tension_original = Get_Random < float >();
      tension = tension_original;
      tension_reference = tension_original;
    }

    stretch = Get_Random < float >();
    activation = Get_Random < float >();
    density = Get_Random < float >();
    fiber_max_stress = Get_Random < float >();

    tension = tension_original;
    for (int i = 0; i < 1; i += 1)
      {
        Tension < float, float, int >(tension, stretch, activation, density,
                                      fiber_max_stress);
      }

    Tension_Reference < float >(tension_reference, stretch, activation, density,
                                fiber_max_stress);
    if (!(Tension_Compare < float >(tension, tension_reference)))
      {
        std::cout << "Failed to confirm unit test for Tension " << std::endl;
        return 1;
      }

  }



  {
    T tension_derivative __attribute__ ((aligned (4)));
    T tension_derivative_reference __attribute__ ((aligned (4)));
    T tension_derivative_original __attribute__ ((aligned (4)));
    T stretch __attribute__ ((aligned (4)));
    T activation __attribute__ ((aligned (4)));
    T density __attribute__ ((aligned (4)));
    T fiber_max_stress __attribute__ ((aligned (4)));


    {
      tension_derivative_original = Get_Random < float >();
      tension_derivative = tension_derivative_original;
      tension_derivative_reference = tension_derivative_original;
    }

    stretch = Get_Random < float >();
    activation = Get_Random < float >();
    density = Get_Random < float >();
    fiber_max_stress = Get_Random < float >();

    tension_derivative = tension_derivative_original;
    for (int i = 0; i < 1; i += 1)
      {
        Tension_Derivative < float, float, int >(tension_derivative, stretch,
                                                 activation, density,
                                                 fiber_max_stress);
      }

    Tension_Derivative_Reference < float >(tension_derivative_reference,
                                           stretch, activation, density,
                                           fiber_max_stress);
    if (!
        (Tension_Derivative_Compare <
         float >(tension_derivative, tension_derivative_reference)))
      {
        std::
          cout << "Failed to confirm unit test for Tension_Derivative " << std::
          endl;
        return 1;
      }

  }



  return 0;

}

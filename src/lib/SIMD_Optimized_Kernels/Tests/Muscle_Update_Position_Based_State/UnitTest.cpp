
#include <cstdlib>
#include <iostream>

struct NEOHOOKEAN_TAG;
struct COROTATED_TAG;

#include "Muscle_Update_Position_Based_State.h"
#include "Muscle_Update_Position_Based_State_Reference.h"

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
    int muscle_id __attribute__ ((aligned (4)));
    T fiber[3] __attribute__ ((aligned (4)));
    T density __attribute__ ((aligned (4)));
    T one_over_h __attribute__ ((aligned (4)));
    T c1 __attribute__ ((aligned (4)));
    T c1_reference __attribute__ ((aligned (4)));
    T c1_original __attribute__ ((aligned (4)));
    T c2 __attribute__ ((aligned (4)));
    T c2_reference __attribute__ ((aligned (4)));
    T c2_original __attribute__ ((aligned (4)));
    T F_fiber[3] __attribute__ ((aligned (4)));
    T F_fiber_reference[3] __attribute__ ((aligned (4)));
    T F_fiber_original[3] __attribute__ ((aligned (4)));
    float *activations __attribute__ ((aligned (4)));
    float *fiber_max_stresses __attribute__ ((aligned (4)));


    for (int __a = 0; __a < 3; __a++)
      for (int __b = 0; __b < 8; __b++)
        u[__a][__b] = Get_Random < float >();
    muscle_id = Get_Random < int >(1, 99);
    for (int __a = 0; __a < 3; __a++)
      fiber[__a] = Get_Random < float >();
    density = Get_Random < float >();
    one_over_h = Get_Random < float >();
    {
      c1_original = Get_Random < float >();
      c1 = c1_original;
      c1_reference = c1_original;
    }
    {
      c2_original = Get_Random < float >();
      c2 = c2_original;
      c2_reference = c2_original;
    }
    for (int __a = 0; __a < 3; __a++)
      {
        F_fiber_original[__a] = Get_Random < float >();
        F_fiber[__a] = F_fiber_original[__a];
        F_fiber_reference[__a] = F_fiber_original[__a];
      }

    activations = new float[100];
    for (int __x__ = 0; __x__ < 100; __x__++)
      activations[__x__] = Get_Random < float >();;
    fiber_max_stresses = new float[100];
    for (int __x__ = 0; __x__ < 100; __x__++)
      fiber_max_stresses[__x__] = Get_Random < float >();;

    c1 = c1_original;
    c2 = c2_original;
    for (int __a = 0; __a < 3; __a++)
      F_fiber[__a] = F_fiber_original[__a];
    for (int i = 0; i < 1; i += 1)
      {
        Muscle_Update_Position_Based_State < float, float, int >(u, muscle_id,
                                                                 fiber, density,
                                                                 one_over_h, c1,
                                                                 c2, F_fiber,
                                                                 activations,
                                                                 fiber_max_stresses);
      }

    Muscle_Update_Position_Based_State_Reference < float >(u, muscle_id, fiber,
                                                           density, one_over_h,
                                                           c1_reference,
                                                           c2_reference,
                                                           F_fiber_reference,
                                                           activations,
                                                           fiber_max_stresses);
    if (!
        (Muscle_Update_Position_Based_State_Compare <
         float >(c1, c2, F_fiber, c1_reference, c2_reference,
                 F_fiber_reference)))
      {
        std::
          cout <<
          "Failed to confirm unit test for Muscle_Update_Position_Based_State "
          << std::endl;
        return 1;
      }

  }



  return 0;

}

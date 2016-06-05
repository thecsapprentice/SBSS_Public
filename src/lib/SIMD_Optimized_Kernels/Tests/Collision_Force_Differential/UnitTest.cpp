
#include <cstdlib>
#include <iostream>

struct NEOHOOKEAN_TAG;
struct COROTATED_TAG;

#include "Collision_Force_Differential.h"
#include "Collision_Force_Differential_Reference.h"

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
    T f[3][8] __attribute__ ((aligned (4)));
    T f_reference[3][8] __attribute__ ((aligned (4)));
    T f_original[3][8] __attribute__ ((aligned (4)));
    T u[3][8] __attribute__ ((aligned (4)));
    T W[3] __attribute__ ((aligned (4)));
    int spring_id __attribute__ ((aligned (4)));
    float *extern_collision_stiffness __attribute__ ((aligned (4)));


    for (int __a = 0; __a < 3; __a++)
      for (int __b = 0; __b < 8; __b++)
        {
          f_original[__a][__b] = Get_Random < float >();
          f[__a][__b] = f_original[__a][__b];
          f_reference[__a][__b] = f_original[__a][__b];
        }
    for (int __a = 0; __a < 3; __a++)
      for (int __b = 0; __b < 8; __b++)
        u[__a][__b] = Get_Random < float >();
    for (int __a = 0; __a < 3; __a++)
      W[__a] = Get_Random < float >();
    spring_id = Get_Random < int >(1, 99);
    extern_collision_stiffness = new float[100];
    for (int __x__ = 0; __x__ < 100; __x__++)
      extern_collision_stiffness[__x__] = Get_Random < float >();;

    for (int __a = 0; __a < 3; __a++)
      for (int __b = 0; __b < 8; __b++)
        f[__a][__b] = f_original[__a][__b];
    for (int i = 0; i < 1; i += 1)
      {
        Collision_Force_Differential < float, float, int >(f, u, W, spring_id,
                                                           extern_collision_stiffness);
      }

    Collision_Force_Differential_Reference < float >(f_reference, u, W,
                                                     spring_id,
                                                     extern_collision_stiffness);
    if (!(Collision_Force_Differential_Compare < float >(f, f_reference)))
      {
        std::
          cout <<
          "Failed to confirm unit test for Collision_Force_Differential " <<
          std::endl;
        return 1;
      }

  }



  return 0;

}

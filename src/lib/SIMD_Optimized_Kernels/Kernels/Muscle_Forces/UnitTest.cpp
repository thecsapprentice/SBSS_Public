
#include <cstdlib>

#include "Muscle_Forces.h"
#include "Muscle_Forces_Reference.h"

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
    T f[3][8];
    T f_reference[3][8];
    T f_original[3][8];
    T fiber[3];
    T Ffiber[3];
    T c1;
    T one_over_h;
    T cell_volume;


    for (int __a = 0; __a < 3; __a++)
      for (int __b = 0; __b < 8; __b++)
        {

          f_original[__a][__b] = Get_Random < float >();

          f[__a][__b] = f_original[__a][__b];

          f_reference[__a][__b] = f_original[__a][__b];
        }
    for (int __a = 0; __a < 3; __a++)
      fiber[__a] = Get_Random < float >();
    for (int __a = 0; __a < 3; __a++)
      Ffiber[__a] = Get_Random < float >();

    c1 = Get_Random < float >();

    one_over_h = Get_Random < float >();

    cell_volume = Get_Random < float >();


    Muscle_Forces < T, T, 1 > (f, fiber, Ffiber, c1, one_over_h, cell_volume);
    Muscle_Forces_Reference < T > (f_reference, fiber, Ffiber, c1, one_over_h,
                                   cell_volume);
    return !(Muscle_Forces_Compare < T > (f, f_reference));
  }



  return 0;

}

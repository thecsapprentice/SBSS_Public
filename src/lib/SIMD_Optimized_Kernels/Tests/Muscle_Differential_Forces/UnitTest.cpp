
#include <cstdlib>
#include <iostream>

struct NEOHOOKEAN_TAG;
struct COROTATED_TAG;

#include "Muscle_Differential_Forces.h"
#include "Muscle_Differential_Forces_Reference.h"

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
    T df[3][8] __attribute__ ((aligned (4)));
    T df_reference[3][8] __attribute__ ((aligned (4)));
    T df_original[3][8] __attribute__ ((aligned (4)));
    T du[3][8] __attribute__ ((aligned (4)));
    T fiber[3] __attribute__ ((aligned (4)));
    T Ffiber[3] __attribute__ ((aligned (4)));
    T c1 __attribute__ ((aligned (4)));
    T c2 __attribute__ ((aligned (4)));
    T one_over_h __attribute__ ((aligned (4)));
    T cell_volume __attribute__ ((aligned (4)));


    for (int __a = 0; __a < 3; __a++)
      for (int __b = 0; __b < 8; __b++)
        {
          df_original[__a][__b] = Get_Random < float >();
          df[__a][__b] = df_original[__a][__b];
          df_reference[__a][__b] = df_original[__a][__b];
        }
    for (int __a = 0; __a < 3; __a++)
      for (int __b = 0; __b < 8; __b++)
        du[__a][__b] = Get_Random < float >();
    for (int __a = 0; __a < 3; __a++)
      fiber[__a] = Get_Random < float >();
    for (int __a = 0; __a < 3; __a++)
      Ffiber[__a] = Get_Random < float >();
    c1 = Get_Random < float >();
    c2 = Get_Random < float >();
    one_over_h = Get_Random < float >();
    cell_volume = Get_Random < float >();

    for (int __a = 0; __a < 3; __a++)
      for (int __b = 0; __b < 8; __b++)
        df[__a][__b] = df_original[__a][__b];
    for (int i = 0; i < 1; i += 1)
      {
        Muscle_Differential_Forces < float, float, int >(df, du, fiber, Ffiber,
                                                         c1, c2, one_over_h,
                                                         cell_volume);
      }

    Muscle_Differential_Forces_Reference < float >(df_reference, du, fiber,
                                                   Ffiber, c1, c2, one_over_h,
                                                   cell_volume);
    if (!(Muscle_Differential_Forces_Compare < float >(df, df_reference)))
      {
        std::
          cout << "Failed to confirm unit test for Muscle_Differential_Forces "
          << std::endl;
        return 1;
      }

  }



  return 0;

}

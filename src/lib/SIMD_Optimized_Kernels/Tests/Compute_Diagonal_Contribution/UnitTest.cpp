
#include <cstdlib>
#include <iostream>

struct NEOHOOKEAN_TAG;
struct COROTATED_TAG;

#include "Compute_Diagonal_Contribution.h"
#include "Compute_Diagonal_Contribution_Reference.h"

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
    T one_over_h __attribute__ ((aligned (4)));
    T mu_stab __attribute__ ((aligned (4)));
    T cell_volume __attribute__ ((aligned (4)));
    T U[9] __attribute__ ((aligned (4)));
    T V[9] __attribute__ ((aligned (4)));
    T dPdF[12] __attribute__ ((aligned (4)));
    T d[3][8] __attribute__ ((aligned (4)));
    T d_reference[3][8] __attribute__ ((aligned (4)));
    T d_original[3][8] __attribute__ ((aligned (4)));



    one_over_h = Get_Random < float >();
    mu_stab = Get_Random < float >();
    cell_volume = Get_Random < float >();
    for (int __a = 0; __a < 9; __a++)
      U[__a] = Get_Random < float >();
    for (int __a = 0; __a < 9; __a++)
      V[__a] = Get_Random < float >();
    for (int __a = 0; __a < 12; __a++)
      dPdF[__a] = Get_Random < float >();
    for (int __a = 0; __a < 3; __a++)
      for (int __b = 0; __b < 8; __b++)
        {
          d_original[__a][__b] = Get_Random < float >();
          d[__a][__b] = d_original[__a][__b];
          d_reference[__a][__b] = d_original[__a][__b];
        }


    for (int __a = 0; __a < 3; __a++)
      for (int __b = 0; __b < 8; __b++)
        d[__a][__b] = d_original[__a][__b];
    for (int i = 0; i < 1; i += 1)
      {
        Compute_Diagonal_Contribution < float, float, int >(one_over_h, mu_stab,
                                                            cell_volume, U, V,
                                                            dPdF, d);
      }

    Compute_Diagonal_Contribution_Reference < float >(one_over_h, mu_stab,
                                                      cell_volume, U, V, dPdF,
                                                      d_reference);
    if (!(Compute_Diagonal_Contribution_Compare < float >(d, d_reference)))
      {
        std::
          cout <<
          "Failed to confirm unit test for Compute_Diagonal_Contribution " <<
          std::endl;
        return 1;
      }

  }



  return 0;

}

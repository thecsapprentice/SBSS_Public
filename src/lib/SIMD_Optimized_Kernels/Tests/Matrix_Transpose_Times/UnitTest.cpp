
#include <cstdlib>
#include <iostream>

struct NEOHOOKEAN_TAG;
struct COROTATED_TAG;

#include "Matrix_Transpose_Times.h"
#include "Matrix_Transpose_Times_Reference.h"

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
    T A[9] __attribute__ ((aligned (4)));
    T B[9] __attribute__ ((aligned (4)));
    T C[9] __attribute__ ((aligned (4)));
    T C_reference[9] __attribute__ ((aligned (4)));
    T C_original[9] __attribute__ ((aligned (4)));


    for (int __a = 0; __a < 9; __a++)
      A[__a] = Get_Random < float >();
    for (int __a = 0; __a < 9; __a++)
      B[__a] = Get_Random < float >();
    for (int __a = 0; __a < 9; __a++)
      {
        C_original[__a] = Get_Random < float >();
        C[__a] = C_original[__a];
        C_reference[__a] = C_original[__a];
      }


    for (int __a = 0; __a < 9; __a++)
      C[__a] = C_original[__a];
    for (int i = 0; i < 1; i += 1)
      {
        Matrix_Transpose_Times < float, float, int >(A, B, C);
      }

    Matrix_Transpose_Times_Reference < float >(A, B, C_reference);
    if (!(Matrix_Transpose_Times_Compare < float >(C, C_reference)))
      {
        std::
          cout << "Failed to confirm unit test for Matrix_Transpose_Times " <<
          std::endl;
        return 1;
      }

  }



  return 0;

}

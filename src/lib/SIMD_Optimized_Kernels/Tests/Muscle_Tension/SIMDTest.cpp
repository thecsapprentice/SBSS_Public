
#include <cstdlib>
#include <iostream>
#include "KernelCommon.h"

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
  typedef float   T;

  int             seed = 1;
  if (argc == 2)
    seed = atoi (argv[1]);
  srand (seed);

  std::cout.precision (10);
  std::cout.setf (std::ios::fixed, std::ios::floatfield);



  {
    std::cout << "Running SIMD Test for Tension " << std::endl;


//=======================================================
//
//        DEFINE ALL VARIABLES USED BY KERNEL
//
//=======================================================

    T               tension[16] __attribute__ ((aligned (64)));
    T               tension_reference[16] __attribute__ ((aligned (64)));
    T               tension_original[16] __attribute__ ((aligned (64)));
    T               stretch[16] __attribute__ ((aligned (64)));
    T               activation[16] __attribute__ ((aligned (64)));
    T               density[16] __attribute__ ((aligned (64)));
    T               fiber_max_stress[16] __attribute__ ((aligned (64)));


    for (int __a = 0; __a < 16; __a++)
    {
      tension_original[__a] = Get_Random < float >();
      tension[__a] = tension_original[__a];
      tension_reference[__a] = tension_original[__a];
    }
    for (int __a = 0; __a < 16; __a++)
      stretch[__a] = Get_Random < float >();
    for (int __a = 0; __a < 16; __a++)
      activation[__a] = Get_Random < float >();
    for (int __a = 0; __a < 16; __a++)
      density[__a] = Get_Random < float >();
    for (int __a = 0; __a < 16; __a++)
      fiber_max_stress[__a] = Get_Random < float >();

//=======================================================
//
//             COMPUTE REFERENCE RESULTS
//
//=======================================================

    T __mtension __attribute__ ((aligned (4)));
    T __mtension_reference __attribute__ ((aligned (4)));
    T __mtension_original __attribute__ ((aligned (4)));
    T __mstretch __attribute__ ((aligned (4)));
    T __mactivation __attribute__ ((aligned (4)));
    T __mdensity __attribute__ ((aligned (4)));
    T __mfiber_max_stress __attribute__ ((aligned (4)));
    for (int k = 0; k < 16; k++)
    {
      __mtension_reference = tension_reference[k];

      __mstretch = stretch[k];
      __mactivation = activation[k];
      __mdensity = density[k];
      __mfiber_max_stress = fiber_max_stress[k];
      Tension < float, float, int >(__mtension_reference, __mstretch,
                                    __mactivation, __mdensity,
                                    __mfiber_max_stress);
      tension_reference[k] = __mtension_reference;
    }

//=======================================================
//
//               COMPUTE SCALAR RESULTS
//
//=======================================================

    {
      typedef         T (&refArray1)[16];
      typedef         T (&refArray2)[16];
      typedef         T (&refArray3)[16];
      typedef         T (&refArray4)[16];
      typedef         T (&refArray5)[16];
      for (int __a = 0; __a < 16; __a++)
        tension[__a] = tension_original[__a];
      for (int i = 0; i < 16; i += 1)
      {
        refArray1       tensionk = reinterpret_cast < refArray1 > (tension[i]);
        refArray2       stretchk = reinterpret_cast < refArray2 > (stretch[i]);
        refArray3       activationk =
          reinterpret_cast < refArray3 > (activation[i]);
        refArray4       densityk = reinterpret_cast < refArray4 > (density[i]);
        refArray5       fiber_max_stressk =
          reinterpret_cast < refArray5 > (fiber_max_stress[i]);
        Tension < float, float[16], int[16] > (tensionk, stretchk, activationk,
                                               densityk, fiber_max_stressk);
      }
      for (int __a = 0; __a < 16; __a++)
        if (std::
            abs ((tension[__a] -
                  tension_reference[__a]) / (tension_reference[__a])) > 1)
        {
          std::cerr << "Mismatch detected in SCALAR implementation" << std::
            endl;
          std::cerr << "Variable tension:" << std::endl;
          std::cerr << "seed=" << seed << ", __a=" << __a << std::endl;
          std::cerr << "tension SCALAR=  " << tension[__a] << std::endl;
          std::cerr << "tension Reference=  " << tension_reference[__a] << std::
            endl;
          std::cerr << "tension Rel Difference=  " << std::
            abs ((tension[__a] -
                  tension_reference[__a]) /
                 (tension_reference[__a])) << std::endl;
          std::cerr << "tension Abs Difference=  " << std::abs (tension[__a] -
                                                                tension_reference
                                                                [__a]) << std::
            endl;
          return 1;
        }

    }

//=======================================================
//
//               COMPUTE SSE RESULTS
//
//=======================================================

#ifdef ENABLE_SSE_INSTRUCTION_SET
    {
      typedef         T (&refArray1)[16];
      typedef         T (&refArray2)[16];
      typedef         T (&refArray3)[16];
      typedef         T (&refArray4)[16];
      typedef         T (&refArray5)[16];
      for (int __a = 0; __a < 16; __a++)
        tension[__a] = tension_original[__a];
      for (int i = 0; i < 16; i += 4)
      {
        refArray1       tensionk = reinterpret_cast < refArray1 > (tension[i]);
        refArray2       stretchk = reinterpret_cast < refArray2 > (stretch[i]);
        refArray3       activationk =
          reinterpret_cast < refArray3 > (activation[i]);
        refArray4       densityk = reinterpret_cast < refArray4 > (density[i]);
        refArray5       fiber_max_stressk =
          reinterpret_cast < refArray5 > (fiber_max_stress[i]);
        Tension < __m128, float[16], int[16] > (tensionk, stretchk, activationk,
                                                densityk, fiber_max_stressk);
      }
      for (int __a = 0; __a < 16; __a++)
        if (std::
            abs ((tension[__a] -
                  tension_reference[__a]) / (tension_reference[__a])) > 1)
        {
          std::cerr << "Mismatch detected in SSE implementation" << std::endl;
          std::cerr << "Variable tension:" << std::endl;
          std::cerr << "seed=" << seed << ", __a=" << __a << std::endl;
          std::cerr << "tension SSE=  " << tension[__a] << std::endl;
          std::cerr << "tension Reference=  " << tension_reference[__a] << std::
            endl;
          std::cerr << "tension Rel Difference=  " << std::
            abs ((tension[__a] -
                  tension_reference[__a]) /
                 (tension_reference[__a])) << std::endl;
          std::cerr << "tension Abs Difference=  " << std::abs (tension[__a] -
                                                                tension_reference
                                                                [__a]) << std::
            endl;
          return 1;
        }

    }
#endif

//=======================================================
//
//               COMPUTE AVX RESULTS
//
//=======================================================

#ifdef ENABLE_AVX_INSTRUCTION_SET
    {
      typedef         T (&refArray1)[16];
      typedef         T (&refArray2)[16];
      typedef         T (&refArray3)[16];
      typedef         T (&refArray4)[16];
      typedef         T (&refArray5)[16];
      for (int __a = 0; __a < 16; __a++)
        tension[__a] = tension_original[__a];
      for (int i = 0; i < 16; i += 8)
      {
        refArray1       tensionk = reinterpret_cast < refArray1 > (tension[i]);
        refArray2       stretchk = reinterpret_cast < refArray2 > (stretch[i]);
        refArray3       activationk =
          reinterpret_cast < refArray3 > (activation[i]);
        refArray4       densityk = reinterpret_cast < refArray4 > (density[i]);
        refArray5       fiber_max_stressk =
          reinterpret_cast < refArray5 > (fiber_max_stress[i]);
        Tension < __m256, float[16], int[16] > (tensionk, stretchk, activationk,
                                                densityk, fiber_max_stressk);
      }
      for (int __a = 0; __a < 16; __a++)
        if (std::
            abs ((tension[__a] -
                  tension_reference[__a]) / (tension_reference[__a])) > 1)
        {
          std::cerr << "Mismatch detected in AVX implementation" << std::endl;
          std::cerr << "Variable tension:" << std::endl;
          std::cerr << "seed=" << seed << ", __a=" << __a << std::endl;
          std::cerr << "tension AVX=  " << tension[__a] << std::endl;
          std::cerr << "tension Reference=  " << tension_reference[__a] << std::
            endl;
          std::cerr << "tension Rel Difference=  " << std::
            abs ((tension[__a] -
                  tension_reference[__a]) /
                 (tension_reference[__a])) << std::endl;
          std::cerr << "tension Abs Difference=  " << std::abs (tension[__a] -
                                                                tension_reference
                                                                [__a]) << std::
            endl;
          return 1;
        }

    }
#endif

//=======================================================
//
//               COMPUTE NEON RESULTS
//
//=======================================================

#ifdef ENABLE_NEON_INSTRUCTION_SET
    {
      typedef         T (&refArray1)[16];
      typedef         T (&refArray2)[16];
      typedef         T (&refArray3)[16];
      typedef         T (&refArray4)[16];
      typedef         T (&refArray5)[16];
      for (int __a = 0; __a < 16; __a++)
        tension[__a] = tension_original[__a];
      for (int i = 0; i < 16; i += 4)
      {
        refArray1       tensionk = reinterpret_cast < refArray1 > (tension[i]);
        refArray2       stretchk = reinterpret_cast < refArray2 > (stretch[i]);
        refArray3       activationk =
          reinterpret_cast < refArray3 > (activation[i]);
        refArray4       densityk = reinterpret_cast < refArray4 > (density[i]);
        refArray5       fiber_max_stressk =
          reinterpret_cast < refArray5 > (fiber_max_stress[i]);
        Tension < float32x4_t, float[16], int[16] > (tensionk, stretchk,
                                                     activationk, densityk,
                                                     fiber_max_stressk);
      }
      for (int __a = 0; __a < 16; __a++)
        if (std::
            abs ((tension[__a] -
                  tension_reference[__a]) / (tension_reference[__a])) > 1)
        {
          std::cerr << "Mismatch detected in NEON implementation" << std::endl;
          std::cerr << "Variable tension:" << std::endl;
          std::cerr << "seed=" << seed << ", __a=" << __a << std::endl;
          std::cerr << "tension NEON=  " << tension[__a] << std::endl;
          std::cerr << "tension Reference=  " << tension_reference[__a] << std::
            endl;
          std::cerr << "tension Rel Difference=  " << std::
            abs ((tension[__a] -
                  tension_reference[__a]) /
                 (tension_reference[__a])) << std::endl;
          std::cerr << "tension Abs Difference=  " << std::abs (tension[__a] -
                                                                tension_reference
                                                                [__a]) << std::
            endl;
          return 1;
        }

    }
#endif

//=======================================================
//
//               COMPUTE MIC RESULTS
//
//=======================================================

#ifdef ENABLE_MIC_INSTRUCTION_SET
    {
      typedef         T (&refArray1)[16];
      typedef         T (&refArray2)[16];
      typedef         T (&refArray3)[16];
      typedef         T (&refArray4)[16];
      typedef         T (&refArray5)[16];
      for (int __a = 0; __a < 16; __a++)
        tension[__a] = tension_original[__a];
      for (int i = 0; i < 16; i += 16)
      {
        refArray1       tensionk = reinterpret_cast < refArray1 > (tension[i]);
        refArray2       stretchk = reinterpret_cast < refArray2 > (stretch[i]);
        refArray3       activationk =
          reinterpret_cast < refArray3 > (activation[i]);
        refArray4       densityk = reinterpret_cast < refArray4 > (density[i]);
        refArray5       fiber_max_stressk =
          reinterpret_cast < refArray5 > (fiber_max_stress[i]);
        Tension < __m512, float[16], int[16] > (tensionk, stretchk, activationk,
                                                densityk, fiber_max_stressk);
      }
      for (int __a = 0; __a < 16; __a++)
        if (std::
            abs ((tension[__a] -
                  tension_reference[__a]) / (tension_reference[__a])) > 1)
        {
          std::cerr << "Mismatch detected in MIC implementation" << std::endl;
          std::cerr << "Variable tension:" << std::endl;
          std::cerr << "seed=" << seed << ", __a=" << __a << std::endl;
          std::cerr << "tension MIC=  " << tension[__a] << std::endl;
          std::cerr << "tension Reference=  " << tension_reference[__a] << std::
            endl;
          std::cerr << "tension Rel Difference=  " << std::
            abs ((tension[__a] -
                  tension_reference[__a]) /
                 (tension_reference[__a])) << std::endl;
          std::cerr << "tension Abs Difference=  " << std::abs (tension[__a] -
                                                                tension_reference
                                                                [__a]) << std::
            endl;
          return 1;
        }

    }
#endif

  }



  {
    std::cout << "Running SIMD Test for Tension_Derivative " << std::endl;


//=======================================================
//
//        DEFINE ALL VARIABLES USED BY KERNEL
//
//=======================================================

    T               tension_derivative[16] __attribute__ ((aligned (64)));
    T               tension_derivative_reference[16]
      __attribute__ ((aligned (64)));
    T               tension_derivative_original[16]
      __attribute__ ((aligned (64)));
    T               stretch[16] __attribute__ ((aligned (64)));
    T               activation[16] __attribute__ ((aligned (64)));
    T               density[16] __attribute__ ((aligned (64)));
    T               fiber_max_stress[16] __attribute__ ((aligned (64)));


    for (int __a = 0; __a < 16; __a++)
    {
      tension_derivative_original[__a] = Get_Random < float >();
      tension_derivative[__a] = tension_derivative_original[__a];
      tension_derivative_reference[__a] = tension_derivative_original[__a];
    }
    for (int __a = 0; __a < 16; __a++)
      stretch[__a] = Get_Random < float >();
    for (int __a = 0; __a < 16; __a++)
      activation[__a] = Get_Random < float >();
    for (int __a = 0; __a < 16; __a++)
      density[__a] = Get_Random < float >();
    for (int __a = 0; __a < 16; __a++)
      fiber_max_stress[__a] = Get_Random < float >();

//=======================================================
//
//             COMPUTE REFERENCE RESULTS
//
//=======================================================

    T __mtension_derivative __attribute__ ((aligned (4)));
    T __mtension_derivative_reference __attribute__ ((aligned (4)));
    T __mtension_derivative_original __attribute__ ((aligned (4)));
    T __mstretch __attribute__ ((aligned (4)));
    T __mactivation __attribute__ ((aligned (4)));
    T __mdensity __attribute__ ((aligned (4)));
    T __mfiber_max_stress __attribute__ ((aligned (4)));
    for (int k = 0; k < 16; k++)
    {
      __mtension_derivative_reference = tension_derivative_reference[k];

      __mstretch = stretch[k];
      __mactivation = activation[k];
      __mdensity = density[k];
      __mfiber_max_stress = fiber_max_stress[k];
      Tension_Derivative < float, float, int >(__mtension_derivative_reference,
                                               __mstretch, __mactivation,
                                               __mdensity, __mfiber_max_stress);
      tension_derivative_reference[k] = __mtension_derivative_reference;
    }

//=======================================================
//
//               COMPUTE SCALAR RESULTS
//
//=======================================================

    {
      typedef         T (&refArray1)[16];
      typedef         T (&refArray2)[16];
      typedef         T (&refArray3)[16];
      typedef         T (&refArray4)[16];
      typedef         T (&refArray5)[16];
      for (int __a = 0; __a < 16; __a++)
        tension_derivative[__a] = tension_derivative_original[__a];
      for (int i = 0; i < 16; i += 1)
      {
        refArray1       tension_derivativek =
          reinterpret_cast < refArray1 > (tension_derivative[i]);
        refArray2       stretchk = reinterpret_cast < refArray2 > (stretch[i]);
        refArray3       activationk =
          reinterpret_cast < refArray3 > (activation[i]);
        refArray4       densityk = reinterpret_cast < refArray4 > (density[i]);
        refArray5       fiber_max_stressk =
          reinterpret_cast < refArray5 > (fiber_max_stress[i]);
        Tension_Derivative < float, float[16], int[16] > (tension_derivativek,
                                                          stretchk, activationk,
                                                          densityk,
                                                          fiber_max_stressk);
      }
      for (int __a = 0; __a < 16; __a++)
        if (std::
            abs ((tension_derivative[__a] -
                  tension_derivative_reference[__a]) /
                 (tension_derivative_reference[__a])) > 1)
        {
          std::cerr << "Mismatch detected in SCALAR implementation" << std::
            endl;
          std::cerr << "Variable tension_derivative:" << std::endl;
          std::cerr << "seed=" << seed << ", __a=" << __a << std::endl;
          std::
            cerr << "tension_derivative SCALAR=  " << tension_derivative[__a] <<
            std::endl;
          std::
            cerr << "tension_derivative Reference=  " <<
            tension_derivative_reference[__a] << std::endl;
          std::cerr << "tension_derivative Rel Difference=  " << std::
            abs ((tension_derivative[__a] -
                  tension_derivative_reference[__a]) /
                 (tension_derivative_reference[__a])) << std::endl;
          std::cerr << "tension_derivative Abs Difference=  " << std::
            abs (tension_derivative[__a] -
                 tension_derivative_reference[__a]) << std::endl;
          return 1;
        }

    }

//=======================================================
//
//               COMPUTE SSE RESULTS
//
//=======================================================

#ifdef ENABLE_SSE_INSTRUCTION_SET
    {
      typedef         T (&refArray1)[16];
      typedef         T (&refArray2)[16];
      typedef         T (&refArray3)[16];
      typedef         T (&refArray4)[16];
      typedef         T (&refArray5)[16];
      for (int __a = 0; __a < 16; __a++)
        tension_derivative[__a] = tension_derivative_original[__a];
      for (int i = 0; i < 16; i += 4)
      {
        refArray1       tension_derivativek =
          reinterpret_cast < refArray1 > (tension_derivative[i]);
        refArray2       stretchk = reinterpret_cast < refArray2 > (stretch[i]);
        refArray3       activationk =
          reinterpret_cast < refArray3 > (activation[i]);
        refArray4       densityk = reinterpret_cast < refArray4 > (density[i]);
        refArray5       fiber_max_stressk =
          reinterpret_cast < refArray5 > (fiber_max_stress[i]);
        Tension_Derivative < __m128, float[16], int[16] > (tension_derivativek,
                                                           stretchk,
                                                           activationk,
                                                           densityk,
                                                           fiber_max_stressk);
      }
      for (int __a = 0; __a < 16; __a++)
        if (std::
            abs ((tension_derivative[__a] -
                  tension_derivative_reference[__a]) /
                 (tension_derivative_reference[__a])) > 1)
        {
          std::cerr << "Mismatch detected in SSE implementation" << std::endl;
          std::cerr << "Variable tension_derivative:" << std::endl;
          std::cerr << "seed=" << seed << ", __a=" << __a << std::endl;
          std::
            cerr << "tension_derivative SSE=  " << tension_derivative[__a] <<
            std::endl;
          std::
            cerr << "tension_derivative Reference=  " <<
            tension_derivative_reference[__a] << std::endl;
          std::cerr << "tension_derivative Rel Difference=  " << std::
            abs ((tension_derivative[__a] -
                  tension_derivative_reference[__a]) /
                 (tension_derivative_reference[__a])) << std::endl;
          std::cerr << "tension_derivative Abs Difference=  " << std::
            abs (tension_derivative[__a] -
                 tension_derivative_reference[__a]) << std::endl;
          return 1;
        }

    }
#endif

//=======================================================
//
//               COMPUTE AVX RESULTS
//
//=======================================================

#ifdef ENABLE_AVX_INSTRUCTION_SET
    {
      typedef         T (&refArray1)[16];
      typedef         T (&refArray2)[16];
      typedef         T (&refArray3)[16];
      typedef         T (&refArray4)[16];
      typedef         T (&refArray5)[16];
      for (int __a = 0; __a < 16; __a++)
        tension_derivative[__a] = tension_derivative_original[__a];
      for (int i = 0; i < 16; i += 8)
      {
        refArray1       tension_derivativek =
          reinterpret_cast < refArray1 > (tension_derivative[i]);
        refArray2       stretchk = reinterpret_cast < refArray2 > (stretch[i]);
        refArray3       activationk =
          reinterpret_cast < refArray3 > (activation[i]);
        refArray4       densityk = reinterpret_cast < refArray4 > (density[i]);
        refArray5       fiber_max_stressk =
          reinterpret_cast < refArray5 > (fiber_max_stress[i]);
        Tension_Derivative < __m256, float[16], int[16] > (tension_derivativek,
                                                           stretchk,
                                                           activationk,
                                                           densityk,
                                                           fiber_max_stressk);
      }
      for (int __a = 0; __a < 16; __a++)
        if (std::
            abs ((tension_derivative[__a] -
                  tension_derivative_reference[__a]) /
                 (tension_derivative_reference[__a])) > 1)
        {
          std::cerr << "Mismatch detected in AVX implementation" << std::endl;
          std::cerr << "Variable tension_derivative:" << std::endl;
          std::cerr << "seed=" << seed << ", __a=" << __a << std::endl;
          std::
            cerr << "tension_derivative AVX=  " << tension_derivative[__a] <<
            std::endl;
          std::
            cerr << "tension_derivative Reference=  " <<
            tension_derivative_reference[__a] << std::endl;
          std::cerr << "tension_derivative Rel Difference=  " << std::
            abs ((tension_derivative[__a] -
                  tension_derivative_reference[__a]) /
                 (tension_derivative_reference[__a])) << std::endl;
          std::cerr << "tension_derivative Abs Difference=  " << std::
            abs (tension_derivative[__a] -
                 tension_derivative_reference[__a]) << std::endl;
          return 1;
        }

    }
#endif

//=======================================================
//
//               COMPUTE NEON RESULTS
//
//=======================================================

#ifdef ENABLE_NEON_INSTRUCTION_SET
    {
      typedef         T (&refArray1)[16];
      typedef         T (&refArray2)[16];
      typedef         T (&refArray3)[16];
      typedef         T (&refArray4)[16];
      typedef         T (&refArray5)[16];
      for (int __a = 0; __a < 16; __a++)
        tension_derivative[__a] = tension_derivative_original[__a];
      for (int i = 0; i < 16; i += 4)
      {
        refArray1       tension_derivativek =
          reinterpret_cast < refArray1 > (tension_derivative[i]);
        refArray2       stretchk = reinterpret_cast < refArray2 > (stretch[i]);
        refArray3       activationk =
          reinterpret_cast < refArray3 > (activation[i]);
        refArray4       densityk = reinterpret_cast < refArray4 > (density[i]);
        refArray5       fiber_max_stressk =
          reinterpret_cast < refArray5 > (fiber_max_stress[i]);
        Tension_Derivative < float32x4_t, float[16],
          int[16] > (tension_derivativek, stretchk, activationk, densityk,
                     fiber_max_stressk);
      }
      for (int __a = 0; __a < 16; __a++)
        if (std::
            abs ((tension_derivative[__a] -
                  tension_derivative_reference[__a]) /
                 (tension_derivative_reference[__a])) > 1)
        {
          std::cerr << "Mismatch detected in NEON implementation" << std::endl;
          std::cerr << "Variable tension_derivative:" << std::endl;
          std::cerr << "seed=" << seed << ", __a=" << __a << std::endl;
          std::
            cerr << "tension_derivative NEON=  " << tension_derivative[__a] <<
            std::endl;
          std::
            cerr << "tension_derivative Reference=  " <<
            tension_derivative_reference[__a] << std::endl;
          std::cerr << "tension_derivative Rel Difference=  " << std::
            abs ((tension_derivative[__a] -
                  tension_derivative_reference[__a]) /
                 (tension_derivative_reference[__a])) << std::endl;
          std::cerr << "tension_derivative Abs Difference=  " << std::
            abs (tension_derivative[__a] -
                 tension_derivative_reference[__a]) << std::endl;
          return 1;
        }

    }
#endif

//=======================================================
//
//               COMPUTE MIC RESULTS
//
//=======================================================

#ifdef ENABLE_MIC_INSTRUCTION_SET
    {
      typedef         T (&refArray1)[16];
      typedef         T (&refArray2)[16];
      typedef         T (&refArray3)[16];
      typedef         T (&refArray4)[16];
      typedef         T (&refArray5)[16];
      for (int __a = 0; __a < 16; __a++)
        tension_derivative[__a] = tension_derivative_original[__a];
      for (int i = 0; i < 16; i += 16)
      {
        refArray1       tension_derivativek =
          reinterpret_cast < refArray1 > (tension_derivative[i]);
        refArray2       stretchk = reinterpret_cast < refArray2 > (stretch[i]);
        refArray3       activationk =
          reinterpret_cast < refArray3 > (activation[i]);
        refArray4       densityk = reinterpret_cast < refArray4 > (density[i]);
        refArray5       fiber_max_stressk =
          reinterpret_cast < refArray5 > (fiber_max_stress[i]);
        Tension_Derivative < __m512, float[16], int[16] > (tension_derivativek,
                                                           stretchk,
                                                           activationk,
                                                           densityk,
                                                           fiber_max_stressk);
      }
      for (int __a = 0; __a < 16; __a++)
        if (std::
            abs ((tension_derivative[__a] -
                  tension_derivative_reference[__a]) /
                 (tension_derivative_reference[__a])) > 1)
        {
          std::cerr << "Mismatch detected in MIC implementation" << std::endl;
          std::cerr << "Variable tension_derivative:" << std::endl;
          std::cerr << "seed=" << seed << ", __a=" << __a << std::endl;
          std::
            cerr << "tension_derivative MIC=  " << tension_derivative[__a] <<
            std::endl;
          std::
            cerr << "tension_derivative Reference=  " <<
            tension_derivative_reference[__a] << std::endl;
          std::cerr << "tension_derivative Rel Difference=  " << std::
            abs ((tension_derivative[__a] -
                  tension_derivative_reference[__a]) /
                 (tension_derivative_reference[__a])) << std::endl;
          std::cerr << "tension_derivative Abs Difference=  " << std::
            abs (tension_derivative[__a] -
                 tension_derivative_reference[__a]) << std::endl;
          return 1;
        }

    }
#endif

  }



  std::cout << "SIMD check successful!" << std::endl;

  return 0;

}

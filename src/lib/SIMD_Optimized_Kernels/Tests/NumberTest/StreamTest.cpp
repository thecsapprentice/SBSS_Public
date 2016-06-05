
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <sys/time.h>
#include "KernelCommon.h"
struct NEOHOOKEAN_TAG;
struct COROTATED_TAG;

#include "NumberTest.h"

#define NUM_TRIALS 1000000

template < class T > T Get_Random (const T a = (T) - 1., const T b = (T) 1.)
{
  return ((b - a) * (T) rand ()) / (T) RAND_MAX + a;
}

struct timeval starttime, stoptime;
void
start_timer ()
{
  gettimeofday (&starttime, NULL);
}

void
stop_timer ()
{
  gettimeofday (&stoptime, NULL);
}

double
get_time ()
{
  return (double) stoptime.tv_sec - (double) starttime.tv_sec +
    (double) 1e-6 *(double) stoptime.tv_usec -
    (double) 1e-6 *(double) starttime.tv_usec;
}

int
main (int argc, char *argv[])
{
  typedef float T;

  std::cout << "Preparing to Run " << NUM_TRIALS << " of all kernels." << std::
    endl;

  int seed = 1;
  if (argc == 2)
    seed = atoi (argv[1]);
  srand (seed);



  {
    std::cout << "Running Stream Test for NumberTest " << std::endl;

//=======================================================
//
//        DEFINE ALL VARIABLES USED BY KERNEL
//
//=======================================================

    T A[16] __attribute__ ((aligned (64)));
    T B[16] __attribute__ ((aligned (64)));
    T C[16] __attribute__ ((aligned (64)));
    T D[16] __attribute__ ((aligned (64)));
    T Add[16] __attribute__ ((aligned (64)));
    T Add_reference[16] __attribute__ ((aligned (64)));
    T Add_original[16] __attribute__ ((aligned (64)));
    T Multiply[16] __attribute__ ((aligned (64)));
    T Multiply_reference[16] __attribute__ ((aligned (64)));
    T Multiply_original[16] __attribute__ ((aligned (64)));
    T Subtract[16] __attribute__ ((aligned (64)));
    T Subtract_reference[16] __attribute__ ((aligned (64)));
    T Subtract_original[16] __attribute__ ((aligned (64)));
    T Divide[16] __attribute__ ((aligned (64)));
    T Divide_reference[16] __attribute__ ((aligned (64)));
    T Divide_original[16] __attribute__ ((aligned (64)));
    T LessThan[16] __attribute__ ((aligned (64)));
    T LessThan_reference[16] __attribute__ ((aligned (64)));
    T LessThan_original[16] __attribute__ ((aligned (64)));
    T GreaterThan[16] __attribute__ ((aligned (64)));
    T GreaterThan_reference[16] __attribute__ ((aligned (64)));
    T GreaterThan_original[16] __attribute__ ((aligned (64)));
    T LessEquals[16] __attribute__ ((aligned (64)));
    T LessEquals_reference[16] __attribute__ ((aligned (64)));
    T LessEquals_original[16] __attribute__ ((aligned (64)));
    T GreaterEquals[16] __attribute__ ((aligned (64)));
    T GreaterEquals_reference[16] __attribute__ ((aligned (64)));
    T GreaterEquals_original[16] __attribute__ ((aligned (64)));
    T Equals[16] __attribute__ ((aligned (64)));
    T Equals_reference[16] __attribute__ ((aligned (64)));
    T Equals_original[16] __attribute__ ((aligned (64)));
    T And[16] __attribute__ ((aligned (64)));
    T And_reference[16] __attribute__ ((aligned (64)));
    T And_original[16] __attribute__ ((aligned (64)));
    T Or[16] __attribute__ ((aligned (64)));
    T Or_reference[16] __attribute__ ((aligned (64)));
    T Or_original[16] __attribute__ ((aligned (64)));
    T Xor[16] __attribute__ ((aligned (64)));
    T Xor_reference[16] __attribute__ ((aligned (64)));
    T Xor_original[16] __attribute__ ((aligned (64)));
    T AndNot[16] __attribute__ ((aligned (64)));
    T AndNot_reference[16] __attribute__ ((aligned (64)));
    T AndNot_original[16] __attribute__ ((aligned (64)));
    T Not[16] __attribute__ ((aligned (64)));
    T Not_reference[16] __attribute__ ((aligned (64)));
    T Not_original[16] __attribute__ ((aligned (64)));
    T Sqrt[16] __attribute__ ((aligned (64)));
    T Sqrt_reference[16] __attribute__ ((aligned (64)));
    T Sqrt_original[16] __attribute__ ((aligned (64)));
    T RSqrt[16] __attribute__ ((aligned (64)));
    T RSqrt_reference[16] __attribute__ ((aligned (64)));
    T RSqrt_original[16] __attribute__ ((aligned (64)));
    T Log[16] __attribute__ ((aligned (64)));
    T Log_reference[16] __attribute__ ((aligned (64)));
    T Log_original[16] __attribute__ ((aligned (64)));
    T Exp[16] __attribute__ ((aligned (64)));
    T Exp_reference[16] __attribute__ ((aligned (64)));
    T Exp_original[16] __attribute__ ((aligned (64)));
    T Inverse[16] __attribute__ ((aligned (64)));
    T Inverse_reference[16] __attribute__ ((aligned (64)));
    T Inverse_original[16] __attribute__ ((aligned (64)));
    T AbsoluteValue[16] __attribute__ ((aligned (64)));
    T AbsoluteValue_reference[16] __attribute__ ((aligned (64)));
    T AbsoluteValue_original[16] __attribute__ ((aligned (64)));
    T Sign[16] __attribute__ ((aligned (64)));
    T Sign_reference[16] __attribute__ ((aligned (64)));
    T Sign_original[16] __attribute__ ((aligned (64)));
    T Minimum[16] __attribute__ ((aligned (64)));
    T Minimum_reference[16] __attribute__ ((aligned (64)));
    T Minimum_original[16] __attribute__ ((aligned (64)));
    T Maximum[16] __attribute__ ((aligned (64)));
    T Maximum_reference[16] __attribute__ ((aligned (64)));
    T Maximum_original[16] __attribute__ ((aligned (64)));
    T Blend[16] __attribute__ ((aligned (64)));
    T Blend_reference[16] __attribute__ ((aligned (64)));
    T Blend_original[16] __attribute__ ((aligned (64)));


    for (int __a = 0; __a < 16; __a++)
      A[__a] = Get_Random < float >();
    for (int __a = 0; __a < 16; __a++)
      B[__a] = Get_Random < float >();
    for (int __a = 0; __a < 16; __a++)
      C[__a] = Get_Random < float >();
    for (int __a = 0; __a < 16; __a++)
      D[__a] = Get_Random < float >();
    for (int __a = 0; __a < 16; __a++)
      {
        Add_original[__a] = Get_Random < float >();
        Add[__a] = Add_original[__a];
        Add_reference[__a] = Add_original[__a];
      }
    for (int __a = 0; __a < 16; __a++)
      {
        Multiply_original[__a] = Get_Random < float >();
        Multiply[__a] = Multiply_original[__a];
        Multiply_reference[__a] = Multiply_original[__a];
      }
    for (int __a = 0; __a < 16; __a++)
      {
        Subtract_original[__a] = Get_Random < float >();
        Subtract[__a] = Subtract_original[__a];
        Subtract_reference[__a] = Subtract_original[__a];
      }
    for (int __a = 0; __a < 16; __a++)
      {
        Divide_original[__a] = Get_Random < float >();
        Divide[__a] = Divide_original[__a];
        Divide_reference[__a] = Divide_original[__a];
      }
    for (int __a = 0; __a < 16; __a++)
      {
        LessThan_original[__a] = Get_Random < float >();
        LessThan[__a] = LessThan_original[__a];
        LessThan_reference[__a] = LessThan_original[__a];
      }
    for (int __a = 0; __a < 16; __a++)
      {
        GreaterThan_original[__a] = Get_Random < float >();
        GreaterThan[__a] = GreaterThan_original[__a];
        GreaterThan_reference[__a] = GreaterThan_original[__a];
      }
    for (int __a = 0; __a < 16; __a++)
      {
        LessEquals_original[__a] = Get_Random < float >();
        LessEquals[__a] = LessEquals_original[__a];
        LessEquals_reference[__a] = LessEquals_original[__a];
      }
    for (int __a = 0; __a < 16; __a++)
      {
        GreaterEquals_original[__a] = Get_Random < float >();
        GreaterEquals[__a] = GreaterEquals_original[__a];
        GreaterEquals_reference[__a] = GreaterEquals_original[__a];
      }
    for (int __a = 0; __a < 16; __a++)
      {
        Equals_original[__a] = Get_Random < float >();
        Equals[__a] = Equals_original[__a];
        Equals_reference[__a] = Equals_original[__a];
      }
    for (int __a = 0; __a < 16; __a++)
      {
        And_original[__a] = Get_Random < float >();
        And[__a] = And_original[__a];
        And_reference[__a] = And_original[__a];
      }
    for (int __a = 0; __a < 16; __a++)
      {
        Or_original[__a] = Get_Random < float >();
        Or[__a] = Or_original[__a];
        Or_reference[__a] = Or_original[__a];
      }
    for (int __a = 0; __a < 16; __a++)
      {
        Xor_original[__a] = Get_Random < float >();
        Xor[__a] = Xor_original[__a];
        Xor_reference[__a] = Xor_original[__a];
      }
    for (int __a = 0; __a < 16; __a++)
      {
        AndNot_original[__a] = Get_Random < float >();
        AndNot[__a] = AndNot_original[__a];
        AndNot_reference[__a] = AndNot_original[__a];
      }
    for (int __a = 0; __a < 16; __a++)
      {
        Not_original[__a] = Get_Random < float >();
        Not[__a] = Not_original[__a];
        Not_reference[__a] = Not_original[__a];
      }
    for (int __a = 0; __a < 16; __a++)
      {
        Sqrt_original[__a] = Get_Random < float >();
        Sqrt[__a] = Sqrt_original[__a];
        Sqrt_reference[__a] = Sqrt_original[__a];
      }
    for (int __a = 0; __a < 16; __a++)
      {
        RSqrt_original[__a] = Get_Random < float >();
        RSqrt[__a] = RSqrt_original[__a];
        RSqrt_reference[__a] = RSqrt_original[__a];
      }
    for (int __a = 0; __a < 16; __a++)
      {
        Log_original[__a] = Get_Random < float >();
        Log[__a] = Log_original[__a];
        Log_reference[__a] = Log_original[__a];
      }
    for (int __a = 0; __a < 16; __a++)
      {
        Exp_original[__a] = Get_Random < float >();
        Exp[__a] = Exp_original[__a];
        Exp_reference[__a] = Exp_original[__a];
      }
    for (int __a = 0; __a < 16; __a++)
      {
        Inverse_original[__a] = Get_Random < float >();
        Inverse[__a] = Inverse_original[__a];
        Inverse_reference[__a] = Inverse_original[__a];
      }
    for (int __a = 0; __a < 16; __a++)
      {
        AbsoluteValue_original[__a] = Get_Random < float >();
        AbsoluteValue[__a] = AbsoluteValue_original[__a];
        AbsoluteValue_reference[__a] = AbsoluteValue_original[__a];
      }
    for (int __a = 0; __a < 16; __a++)
      {
        Sign_original[__a] = Get_Random < float >();
        Sign[__a] = Sign_original[__a];
        Sign_reference[__a] = Sign_original[__a];
      }
    for (int __a = 0; __a < 16; __a++)
      {
        Minimum_original[__a] = Get_Random < float >();
        Minimum[__a] = Minimum_original[__a];
        Minimum_reference[__a] = Minimum_original[__a];
      }
    for (int __a = 0; __a < 16; __a++)
      {
        Maximum_original[__a] = Get_Random < float >();
        Maximum[__a] = Maximum_original[__a];
        Maximum_reference[__a] = Maximum_original[__a];
      }
    for (int __a = 0; __a < 16; __a++)
      {
        Blend_original[__a] = Get_Random < float >();
        Blend[__a] = Blend_original[__a];
        Blend_reference[__a] = Blend_original[__a];
      }


//=======================================================
//
//             COMPUTE SCALAR RESULTS
//
//=======================================================

    {
      std::cout << "	Running " << NUM_TRIALS << " of SCALAR :  ";
      start_timer ();
      for (int n = 0; n < NUM_TRIALS; n++)
        {
          typedef T (&refArray1)[16];
          typedef T (&refArray2)[16];
          typedef T (&refArray3)[16];
          typedef T (&refArray4)[16];
          typedef T (&refArray5)[16];
          typedef T (&refArray6)[16];
          typedef T (&refArray7)[16];
          typedef T (&refArray8)[16];
          typedef T (&refArray9)[16];
          typedef T (&refArray10)[16];
          typedef T (&refArray11)[16];
          typedef T (&refArray12)[16];
          typedef T (&refArray13)[16];
          typedef T (&refArray14)[16];
          typedef T (&refArray15)[16];
          typedef T (&refArray16)[16];
          typedef T (&refArray17)[16];
          typedef T (&refArray18)[16];
          typedef T (&refArray19)[16];
          typedef T (&refArray20)[16];
          typedef T (&refArray21)[16];
          typedef T (&refArray22)[16];
          typedef T (&refArray23)[16];
          typedef T (&refArray24)[16];
          typedef T (&refArray25)[16];
          typedef T (&refArray26)[16];
          typedef T (&refArray27)[16];
          typedef T (&refArray28)[16];
          for (int i = 0; i < 16; i += 1)
            {
              refArray1 Ak = reinterpret_cast < refArray1 > (A[i]);
              refArray2 Bk = reinterpret_cast < refArray2 > (B[i]);
              refArray3 Ck = reinterpret_cast < refArray3 > (C[i]);
              refArray4 Dk = reinterpret_cast < refArray4 > (D[i]);
              refArray5 Addk = reinterpret_cast < refArray5 > (Add[i]);
              refArray6 Multiplyk =
                reinterpret_cast < refArray6 > (Multiply[i]);
              refArray7 Subtractk =
                reinterpret_cast < refArray7 > (Subtract[i]);
              refArray8 Dividek = reinterpret_cast < refArray8 > (Divide[i]);
              refArray9 LessThank =
                reinterpret_cast < refArray9 > (LessThan[i]);
              refArray10 GreaterThank =
                reinterpret_cast < refArray10 > (GreaterThan[i]);
              refArray11 LessEqualsk =
                reinterpret_cast < refArray11 > (LessEquals[i]);
              refArray12 GreaterEqualsk =
                reinterpret_cast < refArray12 > (GreaterEquals[i]);
              refArray13 Equalsk = reinterpret_cast < refArray13 > (Equals[i]);
              refArray14 Andk = reinterpret_cast < refArray14 > (And[i]);
              refArray15 Ork = reinterpret_cast < refArray15 > (Or[i]);
              refArray16 Xork = reinterpret_cast < refArray16 > (Xor[i]);
              refArray17 AndNotk = reinterpret_cast < refArray17 > (AndNot[i]);
              refArray18 Notk = reinterpret_cast < refArray18 > (Not[i]);
              refArray19 Sqrtk = reinterpret_cast < refArray19 > (Sqrt[i]);
              refArray20 RSqrtk = reinterpret_cast < refArray20 > (RSqrt[i]);
              refArray21 Logk = reinterpret_cast < refArray21 > (Log[i]);
              refArray22 Expk = reinterpret_cast < refArray22 > (Exp[i]);
              refArray23 Inversek =
                reinterpret_cast < refArray23 > (Inverse[i]);
              refArray24 AbsoluteValuek =
                reinterpret_cast < refArray24 > (AbsoluteValue[i]);
              refArray25 Signk = reinterpret_cast < refArray25 > (Sign[i]);
              refArray26 Minimumk =
                reinterpret_cast < refArray26 > (Minimum[i]);
              refArray27 Maximumk =
                reinterpret_cast < refArray27 > (Maximum[i]);
              refArray28 Blendk = reinterpret_cast < refArray28 > (Blend[i]);
              NumberTest < float, float[16], int[16] > (Ak, Bk, Ck, Dk, Addk,
                                                        Multiplyk, Subtractk,
                                                        Dividek, LessThank,
                                                        GreaterThank,
                                                        LessEqualsk,
                                                        GreaterEqualsk, Equalsk,
                                                        Andk, Ork, Xork,
                                                        AndNotk, Notk, Sqrtk,
                                                        RSqrtk, Logk, Expk,
                                                        Inversek,
                                                        AbsoluteValuek, Signk,
                                                        Minimumk, Maximumk,
                                                        Blendk);
            }

        }
      stop_timer ();
      std::cout << get_time () << "s" << std::endl;
    }

//=======================================================
//
//             COMPUTE SSE RESULTS
//
//=======================================================

#ifdef ENABLE_SSE_INSTRUCTION_SET
    {
      std::cout << "	Running " << NUM_TRIALS << " of SSE :  ";
      start_timer ();
      for (int n = 0; n < NUM_TRIALS; n++)
        {
          typedef T (&refArray1)[16];
          typedef T (&refArray2)[16];
          typedef T (&refArray3)[16];
          typedef T (&refArray4)[16];
          typedef T (&refArray5)[16];
          typedef T (&refArray6)[16];
          typedef T (&refArray7)[16];
          typedef T (&refArray8)[16];
          typedef T (&refArray9)[16];
          typedef T (&refArray10)[16];
          typedef T (&refArray11)[16];
          typedef T (&refArray12)[16];
          typedef T (&refArray13)[16];
          typedef T (&refArray14)[16];
          typedef T (&refArray15)[16];
          typedef T (&refArray16)[16];
          typedef T (&refArray17)[16];
          typedef T (&refArray18)[16];
          typedef T (&refArray19)[16];
          typedef T (&refArray20)[16];
          typedef T (&refArray21)[16];
          typedef T (&refArray22)[16];
          typedef T (&refArray23)[16];
          typedef T (&refArray24)[16];
          typedef T (&refArray25)[16];
          typedef T (&refArray26)[16];
          typedef T (&refArray27)[16];
          typedef T (&refArray28)[16];
          for (int i = 0; i < 16; i += 4)
            {
              refArray1 Ak = reinterpret_cast < refArray1 > (A[i]);
              refArray2 Bk = reinterpret_cast < refArray2 > (B[i]);
              refArray3 Ck = reinterpret_cast < refArray3 > (C[i]);
              refArray4 Dk = reinterpret_cast < refArray4 > (D[i]);
              refArray5 Addk = reinterpret_cast < refArray5 > (Add[i]);
              refArray6 Multiplyk =
                reinterpret_cast < refArray6 > (Multiply[i]);
              refArray7 Subtractk =
                reinterpret_cast < refArray7 > (Subtract[i]);
              refArray8 Dividek = reinterpret_cast < refArray8 > (Divide[i]);
              refArray9 LessThank =
                reinterpret_cast < refArray9 > (LessThan[i]);
              refArray10 GreaterThank =
                reinterpret_cast < refArray10 > (GreaterThan[i]);
              refArray11 LessEqualsk =
                reinterpret_cast < refArray11 > (LessEquals[i]);
              refArray12 GreaterEqualsk =
                reinterpret_cast < refArray12 > (GreaterEquals[i]);
              refArray13 Equalsk = reinterpret_cast < refArray13 > (Equals[i]);
              refArray14 Andk = reinterpret_cast < refArray14 > (And[i]);
              refArray15 Ork = reinterpret_cast < refArray15 > (Or[i]);
              refArray16 Xork = reinterpret_cast < refArray16 > (Xor[i]);
              refArray17 AndNotk = reinterpret_cast < refArray17 > (AndNot[i]);
              refArray18 Notk = reinterpret_cast < refArray18 > (Not[i]);
              refArray19 Sqrtk = reinterpret_cast < refArray19 > (Sqrt[i]);
              refArray20 RSqrtk = reinterpret_cast < refArray20 > (RSqrt[i]);
              refArray21 Logk = reinterpret_cast < refArray21 > (Log[i]);
              refArray22 Expk = reinterpret_cast < refArray22 > (Exp[i]);
              refArray23 Inversek =
                reinterpret_cast < refArray23 > (Inverse[i]);
              refArray24 AbsoluteValuek =
                reinterpret_cast < refArray24 > (AbsoluteValue[i]);
              refArray25 Signk = reinterpret_cast < refArray25 > (Sign[i]);
              refArray26 Minimumk =
                reinterpret_cast < refArray26 > (Minimum[i]);
              refArray27 Maximumk =
                reinterpret_cast < refArray27 > (Maximum[i]);
              refArray28 Blendk = reinterpret_cast < refArray28 > (Blend[i]);
              NumberTest < __m128, float[16], int[16] > (Ak, Bk, Ck, Dk, Addk,
                                                         Multiplyk, Subtractk,
                                                         Dividek, LessThank,
                                                         GreaterThank,
                                                         LessEqualsk,
                                                         GreaterEqualsk,
                                                         Equalsk, Andk, Ork,
                                                         Xork, AndNotk, Notk,
                                                         Sqrtk, RSqrtk, Logk,
                                                         Expk, Inversek,
                                                         AbsoluteValuek, Signk,
                                                         Minimumk, Maximumk,
                                                         Blendk);
            }

        }
      stop_timer ();
      std::cout << get_time () << "s" << std::endl;
    }
#endif

//=======================================================
//
//             COMPUTE AVX RESULTS
//
//=======================================================

#ifdef ENABLE_AVX_INSTRUCTION_SET
    {
      std::cout << "	Running " << NUM_TRIALS << " of AVX :  ";
      start_timer ();
      for (int n = 0; n < NUM_TRIALS; n++)
        {
          typedef T (&refArray1)[16];
          typedef T (&refArray2)[16];
          typedef T (&refArray3)[16];
          typedef T (&refArray4)[16];
          typedef T (&refArray5)[16];
          typedef T (&refArray6)[16];
          typedef T (&refArray7)[16];
          typedef T (&refArray8)[16];
          typedef T (&refArray9)[16];
          typedef T (&refArray10)[16];
          typedef T (&refArray11)[16];
          typedef T (&refArray12)[16];
          typedef T (&refArray13)[16];
          typedef T (&refArray14)[16];
          typedef T (&refArray15)[16];
          typedef T (&refArray16)[16];
          typedef T (&refArray17)[16];
          typedef T (&refArray18)[16];
          typedef T (&refArray19)[16];
          typedef T (&refArray20)[16];
          typedef T (&refArray21)[16];
          typedef T (&refArray22)[16];
          typedef T (&refArray23)[16];
          typedef T (&refArray24)[16];
          typedef T (&refArray25)[16];
          typedef T (&refArray26)[16];
          typedef T (&refArray27)[16];
          typedef T (&refArray28)[16];
          for (int i = 0; i < 16; i += 8)
            {
              refArray1 Ak = reinterpret_cast < refArray1 > (A[i]);
              refArray2 Bk = reinterpret_cast < refArray2 > (B[i]);
              refArray3 Ck = reinterpret_cast < refArray3 > (C[i]);
              refArray4 Dk = reinterpret_cast < refArray4 > (D[i]);
              refArray5 Addk = reinterpret_cast < refArray5 > (Add[i]);
              refArray6 Multiplyk =
                reinterpret_cast < refArray6 > (Multiply[i]);
              refArray7 Subtractk =
                reinterpret_cast < refArray7 > (Subtract[i]);
              refArray8 Dividek = reinterpret_cast < refArray8 > (Divide[i]);
              refArray9 LessThank =
                reinterpret_cast < refArray9 > (LessThan[i]);
              refArray10 GreaterThank =
                reinterpret_cast < refArray10 > (GreaterThan[i]);
              refArray11 LessEqualsk =
                reinterpret_cast < refArray11 > (LessEquals[i]);
              refArray12 GreaterEqualsk =
                reinterpret_cast < refArray12 > (GreaterEquals[i]);
              refArray13 Equalsk = reinterpret_cast < refArray13 > (Equals[i]);
              refArray14 Andk = reinterpret_cast < refArray14 > (And[i]);
              refArray15 Ork = reinterpret_cast < refArray15 > (Or[i]);
              refArray16 Xork = reinterpret_cast < refArray16 > (Xor[i]);
              refArray17 AndNotk = reinterpret_cast < refArray17 > (AndNot[i]);
              refArray18 Notk = reinterpret_cast < refArray18 > (Not[i]);
              refArray19 Sqrtk = reinterpret_cast < refArray19 > (Sqrt[i]);
              refArray20 RSqrtk = reinterpret_cast < refArray20 > (RSqrt[i]);
              refArray21 Logk = reinterpret_cast < refArray21 > (Log[i]);
              refArray22 Expk = reinterpret_cast < refArray22 > (Exp[i]);
              refArray23 Inversek =
                reinterpret_cast < refArray23 > (Inverse[i]);
              refArray24 AbsoluteValuek =
                reinterpret_cast < refArray24 > (AbsoluteValue[i]);
              refArray25 Signk = reinterpret_cast < refArray25 > (Sign[i]);
              refArray26 Minimumk =
                reinterpret_cast < refArray26 > (Minimum[i]);
              refArray27 Maximumk =
                reinterpret_cast < refArray27 > (Maximum[i]);
              refArray28 Blendk = reinterpret_cast < refArray28 > (Blend[i]);
              NumberTest < __m256, float[16], int[16] > (Ak, Bk, Ck, Dk, Addk,
                                                         Multiplyk, Subtractk,
                                                         Dividek, LessThank,
                                                         GreaterThank,
                                                         LessEqualsk,
                                                         GreaterEqualsk,
                                                         Equalsk, Andk, Ork,
                                                         Xork, AndNotk, Notk,
                                                         Sqrtk, RSqrtk, Logk,
                                                         Expk, Inversek,
                                                         AbsoluteValuek, Signk,
                                                         Minimumk, Maximumk,
                                                         Blendk);
            }

        }
      stop_timer ();
      std::cout << get_time () << "s" << std::endl;
    }
#endif

//=======================================================
//
//             COMPUTE NEON RESULTS
//
//=======================================================

#ifdef ENABLE_NEON_INSTRUCTION_SET
    {
      std::cout << "	Running " << NUM_TRIALS << " of NEON :  ";
      start_timer ();
      for (int n = 0; n < NUM_TRIALS; n++)
        {
          typedef T (&refArray1)[16];
          typedef T (&refArray2)[16];
          typedef T (&refArray3)[16];
          typedef T (&refArray4)[16];
          typedef T (&refArray5)[16];
          typedef T (&refArray6)[16];
          typedef T (&refArray7)[16];
          typedef T (&refArray8)[16];
          typedef T (&refArray9)[16];
          typedef T (&refArray10)[16];
          typedef T (&refArray11)[16];
          typedef T (&refArray12)[16];
          typedef T (&refArray13)[16];
          typedef T (&refArray14)[16];
          typedef T (&refArray15)[16];
          typedef T (&refArray16)[16];
          typedef T (&refArray17)[16];
          typedef T (&refArray18)[16];
          typedef T (&refArray19)[16];
          typedef T (&refArray20)[16];
          typedef T (&refArray21)[16];
          typedef T (&refArray22)[16];
          typedef T (&refArray23)[16];
          typedef T (&refArray24)[16];
          typedef T (&refArray25)[16];
          typedef T (&refArray26)[16];
          typedef T (&refArray27)[16];
          typedef T (&refArray28)[16];
          for (int i = 0; i < 16; i += 4)
            {
              refArray1 Ak = reinterpret_cast < refArray1 > (A[i]);
              refArray2 Bk = reinterpret_cast < refArray2 > (B[i]);
              refArray3 Ck = reinterpret_cast < refArray3 > (C[i]);
              refArray4 Dk = reinterpret_cast < refArray4 > (D[i]);
              refArray5 Addk = reinterpret_cast < refArray5 > (Add[i]);
              refArray6 Multiplyk =
                reinterpret_cast < refArray6 > (Multiply[i]);
              refArray7 Subtractk =
                reinterpret_cast < refArray7 > (Subtract[i]);
              refArray8 Dividek = reinterpret_cast < refArray8 > (Divide[i]);
              refArray9 LessThank =
                reinterpret_cast < refArray9 > (LessThan[i]);
              refArray10 GreaterThank =
                reinterpret_cast < refArray10 > (GreaterThan[i]);
              refArray11 LessEqualsk =
                reinterpret_cast < refArray11 > (LessEquals[i]);
              refArray12 GreaterEqualsk =
                reinterpret_cast < refArray12 > (GreaterEquals[i]);
              refArray13 Equalsk = reinterpret_cast < refArray13 > (Equals[i]);
              refArray14 Andk = reinterpret_cast < refArray14 > (And[i]);
              refArray15 Ork = reinterpret_cast < refArray15 > (Or[i]);
              refArray16 Xork = reinterpret_cast < refArray16 > (Xor[i]);
              refArray17 AndNotk = reinterpret_cast < refArray17 > (AndNot[i]);
              refArray18 Notk = reinterpret_cast < refArray18 > (Not[i]);
              refArray19 Sqrtk = reinterpret_cast < refArray19 > (Sqrt[i]);
              refArray20 RSqrtk = reinterpret_cast < refArray20 > (RSqrt[i]);
              refArray21 Logk = reinterpret_cast < refArray21 > (Log[i]);
              refArray22 Expk = reinterpret_cast < refArray22 > (Exp[i]);
              refArray23 Inversek =
                reinterpret_cast < refArray23 > (Inverse[i]);
              refArray24 AbsoluteValuek =
                reinterpret_cast < refArray24 > (AbsoluteValue[i]);
              refArray25 Signk = reinterpret_cast < refArray25 > (Sign[i]);
              refArray26 Minimumk =
                reinterpret_cast < refArray26 > (Minimum[i]);
              refArray27 Maximumk =
                reinterpret_cast < refArray27 > (Maximum[i]);
              refArray28 Blendk = reinterpret_cast < refArray28 > (Blend[i]);
              NumberTest < float32x4_t, float[16], int[16] > (Ak, Bk, Ck, Dk,
                                                              Addk, Multiplyk,
                                                              Subtractk,
                                                              Dividek,
                                                              LessThank,
                                                              GreaterThank,
                                                              LessEqualsk,
                                                              GreaterEqualsk,
                                                              Equalsk, Andk,
                                                              Ork, Xork,
                                                              AndNotk, Notk,
                                                              Sqrtk, RSqrtk,
                                                              Logk, Expk,
                                                              Inversek,
                                                              AbsoluteValuek,
                                                              Signk, Minimumk,
                                                              Maximumk, Blendk);
            }

        }
      stop_timer ();
      std::cout << get_time () << "s" << std::endl;
    }
#endif

//=======================================================
//
//             COMPUTE MIC RESULTS
//
//=======================================================

#ifdef ENABLE_MIC_INSTRUCTION_SET
    {
      std::cout << "	Running " << NUM_TRIALS << " of MIC :  ";
      start_timer ();
      for (int n = 0; n < NUM_TRIALS; n++)
        {
          typedef T (&refArray1)[16];
          typedef T (&refArray2)[16];
          typedef T (&refArray3)[16];
          typedef T (&refArray4)[16];
          typedef T (&refArray5)[16];
          typedef T (&refArray6)[16];
          typedef T (&refArray7)[16];
          typedef T (&refArray8)[16];
          typedef T (&refArray9)[16];
          typedef T (&refArray10)[16];
          typedef T (&refArray11)[16];
          typedef T (&refArray12)[16];
          typedef T (&refArray13)[16];
          typedef T (&refArray14)[16];
          typedef T (&refArray15)[16];
          typedef T (&refArray16)[16];
          typedef T (&refArray17)[16];
          typedef T (&refArray18)[16];
          typedef T (&refArray19)[16];
          typedef T (&refArray20)[16];
          typedef T (&refArray21)[16];
          typedef T (&refArray22)[16];
          typedef T (&refArray23)[16];
          typedef T (&refArray24)[16];
          typedef T (&refArray25)[16];
          typedef T (&refArray26)[16];
          typedef T (&refArray27)[16];
          typedef T (&refArray28)[16];
          for (int i = 0; i < 16; i += 16)
            {
              refArray1 Ak = reinterpret_cast < refArray1 > (A[i]);
              refArray2 Bk = reinterpret_cast < refArray2 > (B[i]);
              refArray3 Ck = reinterpret_cast < refArray3 > (C[i]);
              refArray4 Dk = reinterpret_cast < refArray4 > (D[i]);
              refArray5 Addk = reinterpret_cast < refArray5 > (Add[i]);
              refArray6 Multiplyk =
                reinterpret_cast < refArray6 > (Multiply[i]);
              refArray7 Subtractk =
                reinterpret_cast < refArray7 > (Subtract[i]);
              refArray8 Dividek = reinterpret_cast < refArray8 > (Divide[i]);
              refArray9 LessThank =
                reinterpret_cast < refArray9 > (LessThan[i]);
              refArray10 GreaterThank =
                reinterpret_cast < refArray10 > (GreaterThan[i]);
              refArray11 LessEqualsk =
                reinterpret_cast < refArray11 > (LessEquals[i]);
              refArray12 GreaterEqualsk =
                reinterpret_cast < refArray12 > (GreaterEquals[i]);
              refArray13 Equalsk = reinterpret_cast < refArray13 > (Equals[i]);
              refArray14 Andk = reinterpret_cast < refArray14 > (And[i]);
              refArray15 Ork = reinterpret_cast < refArray15 > (Or[i]);
              refArray16 Xork = reinterpret_cast < refArray16 > (Xor[i]);
              refArray17 AndNotk = reinterpret_cast < refArray17 > (AndNot[i]);
              refArray18 Notk = reinterpret_cast < refArray18 > (Not[i]);
              refArray19 Sqrtk = reinterpret_cast < refArray19 > (Sqrt[i]);
              refArray20 RSqrtk = reinterpret_cast < refArray20 > (RSqrt[i]);
              refArray21 Logk = reinterpret_cast < refArray21 > (Log[i]);
              refArray22 Expk = reinterpret_cast < refArray22 > (Exp[i]);
              refArray23 Inversek =
                reinterpret_cast < refArray23 > (Inverse[i]);
              refArray24 AbsoluteValuek =
                reinterpret_cast < refArray24 > (AbsoluteValue[i]);
              refArray25 Signk = reinterpret_cast < refArray25 > (Sign[i]);
              refArray26 Minimumk =
                reinterpret_cast < refArray26 > (Minimum[i]);
              refArray27 Maximumk =
                reinterpret_cast < refArray27 > (Maximum[i]);
              refArray28 Blendk = reinterpret_cast < refArray28 > (Blend[i]);
              NumberTest < __m512, float[16], int[16] > (Ak, Bk, Ck, Dk, Addk,
                                                         Multiplyk, Subtractk,
                                                         Dividek, LessThank,
                                                         GreaterThank,
                                                         LessEqualsk,
                                                         GreaterEqualsk,
                                                         Equalsk, Andk, Ork,
                                                         Xork, AndNotk, Notk,
                                                         Sqrtk, RSqrtk, Logk,
                                                         Expk, Inversek,
                                                         AbsoluteValuek, Signk,
                                                         Minimumk, Maximumk,
                                                         Blendk);
            }

        }
      stop_timer ();
      std::cout << get_time () << "s" << std::endl;
    }
#endif

  }



  return 0;

}

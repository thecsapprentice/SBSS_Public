
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <sys/time.h>
#include "KernelCommon.h"
struct NEOHOOKEAN_TAG;
struct COROTATED_TAG;

#include "NumberTest.h"

#include <Thread_Queueing/PTHREAD_QUEUE.h>
#include <Kernel_Serial_Base_Helper.h>


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


template < int SIZE > class NumberTest_SCALAR_None
{
private:
  // Generate Variables Here
  float *_local_A;
  float *_local_B;
  float *_local_C;
  float *_local_D;
  float *_local_Add;
  float *_local_Multiply;
  float *_local_Subtract;
  float *_local_Divide;
  float *_local_LessThan;
  float *_local_GreaterThan;
  float *_local_LessEquals;
  float *_local_GreaterEquals;
  float *_local_Equals;
  float *_local_And;
  float *_local_Or;
  float *_local_Xor;
  float *_local_AndNot;
  float *_local_Not;
  float *_local_Sqrt;
  float *_local_RSqrt;
  float *_local_Log;
  float *_local_Exp;
  float *_local_Inverse;
  float *_local_AbsoluteValue;
  float *_local_Sign;
  float *_local_Minimum;
  float *_local_Maximum;
  float *_local_Blend;

public:
    explicit NumberTest_SCALAR_None (float *A_in, float *B_in, float *C_in,
                                     float *D_in, float *Add_in,
                                     float *Multiply_in, float *Subtract_in,
                                     float *Divide_in, float *LessThan_in,
                                     float *GreaterThan_in,
                                     float *LessEquals_in,
                                     float *GreaterEquals_in, float *Equals_in,
                                     float *And_in, float *Or_in, float *Xor_in,
                                     float *AndNot_in, float *Not_in,
                                     float *Sqrt_in, float *RSqrt_in,
                                     float *Log_in, float *Exp_in,
                                     float *Inverse_in, float *AbsoluteValue_in,
                                     float *Sign_in, float *Minimum_in,
                                     float *Maximum_in,
                                     float *Blend_in):_local_A (A_in),
    _local_B (B_in), _local_C (C_in), _local_D (D_in), _local_Add (Add_in),
    _local_Multiply (Multiply_in), _local_Subtract (Subtract_in),
    _local_Divide (Divide_in), _local_LessThan (LessThan_in),
    _local_GreaterThan (GreaterThan_in), _local_LessEquals (LessEquals_in),
    _local_GreaterEquals (GreaterEquals_in), _local_Equals (Equals_in),
    _local_And (And_in), _local_Or (Or_in), _local_Xor (Xor_in),
    _local_AndNot (AndNot_in), _local_Not (Not_in), _local_Sqrt (Sqrt_in),
    _local_RSqrt (RSqrt_in), _local_Log (Log_in), _local_Exp (Exp_in),
    _local_Inverse (Inverse_in), _local_AbsoluteValue (AbsoluteValue_in),
    _local_Sign (Sign_in), _local_Minimum (Minimum_in),
    _local_Maximum (Maximum_in), _local_Blend (Blend_in)
  {
  }
  void Execute (int index)
  {
    // full array typedefs
    //typedef int (&refArray)[SIZE][3][8];
    typedef float (&fullArray1)[SIZE][16];
    typedef float (&fullArray2)[SIZE][16];
    typedef float (&fullArray3)[SIZE][16];
    typedef float (&fullArray4)[SIZE][16];
    typedef float (&fullArray5)[SIZE][16];
    typedef float (&fullArray6)[SIZE][16];
    typedef float (&fullArray7)[SIZE][16];
    typedef float (&fullArray8)[SIZE][16];
    typedef float (&fullArray9)[SIZE][16];
    typedef float (&fullArray10)[SIZE][16];
    typedef float (&fullArray11)[SIZE][16];
    typedef float (&fullArray12)[SIZE][16];
    typedef float (&fullArray13)[SIZE][16];
    typedef float (&fullArray14)[SIZE][16];
    typedef float (&fullArray15)[SIZE][16];
    typedef float (&fullArray16)[SIZE][16];
    typedef float (&fullArray17)[SIZE][16];
    typedef float (&fullArray18)[SIZE][16];
    typedef float (&fullArray19)[SIZE][16];
    typedef float (&fullArray20)[SIZE][16];
    typedef float (&fullArray21)[SIZE][16];
    typedef float (&fullArray22)[SIZE][16];
    typedef float (&fullArray23)[SIZE][16];
    typedef float (&fullArray24)[SIZE][16];
    typedef float (&fullArray25)[SIZE][16];
    typedef float (&fullArray26)[SIZE][16];
    typedef float (&fullArray27)[SIZE][16];
    typedef float (&fullArray28)[SIZE][16];

    // full array extractions
    //refArray _rA = reinterpret_cast < refArray >(*_A);
    fullArray1 _rA = reinterpret_cast < fullArray1 > (*_local_A);
    fullArray2 _rB = reinterpret_cast < fullArray2 > (*_local_B);
    fullArray3 _rC = reinterpret_cast < fullArray3 > (*_local_C);
    fullArray4 _rD = reinterpret_cast < fullArray4 > (*_local_D);
    fullArray5 _rAdd = reinterpret_cast < fullArray5 > (*_local_Add);
    fullArray6 _rMultiply = reinterpret_cast < fullArray6 > (*_local_Multiply);
    fullArray7 _rSubtract = reinterpret_cast < fullArray7 > (*_local_Subtract);
    fullArray8 _rDivide = reinterpret_cast < fullArray8 > (*_local_Divide);
    fullArray9 _rLessThan = reinterpret_cast < fullArray9 > (*_local_LessThan);
    fullArray10 _rGreaterThan =
      reinterpret_cast < fullArray10 > (*_local_GreaterThan);
    fullArray11 _rLessEquals =
      reinterpret_cast < fullArray11 > (*_local_LessEquals);
    fullArray12 _rGreaterEquals =
      reinterpret_cast < fullArray12 > (*_local_GreaterEquals);
    fullArray13 _rEquals = reinterpret_cast < fullArray13 > (*_local_Equals);
    fullArray14 _rAnd = reinterpret_cast < fullArray14 > (*_local_And);
    fullArray15 _rOr = reinterpret_cast < fullArray15 > (*_local_Or);
    fullArray16 _rXor = reinterpret_cast < fullArray16 > (*_local_Xor);
    fullArray17 _rAndNot = reinterpret_cast < fullArray17 > (*_local_AndNot);
    fullArray18 _rNot = reinterpret_cast < fullArray18 > (*_local_Not);
    fullArray19 _rSqrt = reinterpret_cast < fullArray19 > (*_local_Sqrt);
    fullArray20 _rRSqrt = reinterpret_cast < fullArray20 > (*_local_RSqrt);
    fullArray21 _rLog = reinterpret_cast < fullArray21 > (*_local_Log);
    fullArray22 _rExp = reinterpret_cast < fullArray22 > (*_local_Exp);
    fullArray23 _rInverse = reinterpret_cast < fullArray23 > (*_local_Inverse);
    fullArray24 _rAbsoluteValue =
      reinterpret_cast < fullArray24 > (*_local_AbsoluteValue);
    fullArray25 _rSign = reinterpret_cast < fullArray25 > (*_local_Sign);
    fullArray26 _rMinimum = reinterpret_cast < fullArray26 > (*_local_Minimum);
    fullArray27 _rMaximum = reinterpret_cast < fullArray27 > (*_local_Maximum);
    fullArray28 _rBlend = reinterpret_cast < fullArray28 > (*_local_Blend);

    const int ChunkSize = 1;
    // chunk typedef
    //typedef int (&refArray1)[3][8];
    typedef float (&refArray1)[16];
    typedef float (&refArray2)[16];
    typedef float (&refArray3)[16];
    typedef float (&refArray4)[16];
    typedef float (&refArray5)[16];
    typedef float (&refArray6)[16];
    typedef float (&refArray7)[16];
    typedef float (&refArray8)[16];
    typedef float (&refArray9)[16];
    typedef float (&refArray10)[16];
    typedef float (&refArray11)[16];
    typedef float (&refArray12)[16];
    typedef float (&refArray13)[16];
    typedef float (&refArray14)[16];
    typedef float (&refArray15)[16];
    typedef float (&refArray16)[16];
    typedef float (&refArray17)[16];
    typedef float (&refArray18)[16];
    typedef float (&refArray19)[16];
    typedef float (&refArray20)[16];
    typedef float (&refArray21)[16];
    typedef float (&refArray22)[16];
    typedef float (&refArray23)[16];
    typedef float (&refArray24)[16];
    typedef float (&refArray25)[16];
    typedef float (&refArray26)[16];
    typedef float (&refArray27)[16];
    typedef float (&refArray28)[16];

    for (int chunk_offset = 0; chunk_offset < 16; chunk_offset += ChunkSize)
      {
        //refArray1 Ak = reinterpret_cast < refArray1 > (_rA[entry][0][chunk_offset]);
        refArray1 Ak =
          reinterpret_cast < refArray1 > (_rA[index][chunk_offset]);
        refArray2 Bk =
          reinterpret_cast < refArray2 > (_rB[index][chunk_offset]);
        refArray3 Ck =
          reinterpret_cast < refArray3 > (_rC[index][chunk_offset]);
        refArray4 Dk =
          reinterpret_cast < refArray4 > (_rD[index][chunk_offset]);
        refArray5 Addk =
          reinterpret_cast < refArray5 > (_rAdd[index][chunk_offset]);
        refArray6 Multiplyk =
          reinterpret_cast < refArray6 > (_rMultiply[index][chunk_offset]);
        refArray7 Subtractk =
          reinterpret_cast < refArray7 > (_rSubtract[index][chunk_offset]);
        refArray8 Dividek =
          reinterpret_cast < refArray8 > (_rDivide[index][chunk_offset]);
        refArray9 LessThank =
          reinterpret_cast < refArray9 > (_rLessThan[index][chunk_offset]);
        refArray10 GreaterThank =
          reinterpret_cast < refArray10 > (_rGreaterThan[index][chunk_offset]);
        refArray11 LessEqualsk =
          reinterpret_cast < refArray11 > (_rLessEquals[index][chunk_offset]);
        refArray12 GreaterEqualsk =
          reinterpret_cast < refArray12 >
          (_rGreaterEquals[index][chunk_offset]);
        refArray13 Equalsk =
          reinterpret_cast < refArray13 > (_rEquals[index][chunk_offset]);
        refArray14 Andk =
          reinterpret_cast < refArray14 > (_rAnd[index][chunk_offset]);
        refArray15 Ork =
          reinterpret_cast < refArray15 > (_rOr[index][chunk_offset]);
        refArray16 Xork =
          reinterpret_cast < refArray16 > (_rXor[index][chunk_offset]);
        refArray17 AndNotk =
          reinterpret_cast < refArray17 > (_rAndNot[index][chunk_offset]);
        refArray18 Notk =
          reinterpret_cast < refArray18 > (_rNot[index][chunk_offset]);
        refArray19 Sqrtk =
          reinterpret_cast < refArray19 > (_rSqrt[index][chunk_offset]);
        refArray20 RSqrtk =
          reinterpret_cast < refArray20 > (_rRSqrt[index][chunk_offset]);
        refArray21 Logk =
          reinterpret_cast < refArray21 > (_rLog[index][chunk_offset]);
        refArray22 Expk =
          reinterpret_cast < refArray22 > (_rExp[index][chunk_offset]);
        refArray23 Inversek =
          reinterpret_cast < refArray23 > (_rInverse[index][chunk_offset]);
        refArray24 AbsoluteValuek =
          reinterpret_cast < refArray24 >
          (_rAbsoluteValue[index][chunk_offset]);
        refArray25 Signk =
          reinterpret_cast < refArray25 > (_rSign[index][chunk_offset]);
        refArray26 Minimumk =
          reinterpret_cast < refArray26 > (_rMinimum[index][chunk_offset]);
        refArray27 Maximumk =
          reinterpret_cast < refArray27 > (_rMaximum[index][chunk_offset]);
        refArray28 Blendk =
          reinterpret_cast < refArray28 > (_rBlend[index][chunk_offset]);

        NumberTest < float, float[16], int[16] > (Ak, Bk, Ck, Dk, Addk,
                                                  Multiplyk, Subtractk, Dividek,
                                                  LessThank, GreaterThank,
                                                  LessEqualsk, GreaterEqualsk,
                                                  Equalsk, Andk, Ork, Xork,
                                                  AndNotk, Notk, Sqrtk, RSqrtk,
                                                  Logk, Expk, Inversek,
                                                  AbsoluteValuek, Signk,
                                                  Minimumk, Maximumk, Blendk);
      }

  }
};


#ifdef ENABLE_SSE_INSTRUCTION_SET

template < int SIZE > class NumberTest_SSE_None
{
private:
  // Generate Variables Here
  float *_local_A;
  float *_local_B;
  float *_local_C;
  float *_local_D;
  float *_local_Add;
  float *_local_Multiply;
  float *_local_Subtract;
  float *_local_Divide;
  float *_local_LessThan;
  float *_local_GreaterThan;
  float *_local_LessEquals;
  float *_local_GreaterEquals;
  float *_local_Equals;
  float *_local_And;
  float *_local_Or;
  float *_local_Xor;
  float *_local_AndNot;
  float *_local_Not;
  float *_local_Sqrt;
  float *_local_RSqrt;
  float *_local_Log;
  float *_local_Exp;
  float *_local_Inverse;
  float *_local_AbsoluteValue;
  float *_local_Sign;
  float *_local_Minimum;
  float *_local_Maximum;
  float *_local_Blend;

public:
    explicit NumberTest_SSE_None (float *A_in, float *B_in, float *C_in,
                                  float *D_in, float *Add_in,
                                  float *Multiply_in, float *Subtract_in,
                                  float *Divide_in, float *LessThan_in,
                                  float *GreaterThan_in, float *LessEquals_in,
                                  float *GreaterEquals_in, float *Equals_in,
                                  float *And_in, float *Or_in, float *Xor_in,
                                  float *AndNot_in, float *Not_in,
                                  float *Sqrt_in, float *RSqrt_in,
                                  float *Log_in, float *Exp_in,
                                  float *Inverse_in, float *AbsoluteValue_in,
                                  float *Sign_in, float *Minimum_in,
                                  float *Maximum_in,
                                  float *Blend_in):_local_A (A_in),
    _local_B (B_in), _local_C (C_in), _local_D (D_in), _local_Add (Add_in),
    _local_Multiply (Multiply_in), _local_Subtract (Subtract_in),
    _local_Divide (Divide_in), _local_LessThan (LessThan_in),
    _local_GreaterThan (GreaterThan_in), _local_LessEquals (LessEquals_in),
    _local_GreaterEquals (GreaterEquals_in), _local_Equals (Equals_in),
    _local_And (And_in), _local_Or (Or_in), _local_Xor (Xor_in),
    _local_AndNot (AndNot_in), _local_Not (Not_in), _local_Sqrt (Sqrt_in),
    _local_RSqrt (RSqrt_in), _local_Log (Log_in), _local_Exp (Exp_in),
    _local_Inverse (Inverse_in), _local_AbsoluteValue (AbsoluteValue_in),
    _local_Sign (Sign_in), _local_Minimum (Minimum_in),
    _local_Maximum (Maximum_in), _local_Blend (Blend_in)
  {
  }
  void Execute (int index)
  {
    // full array typedefs
    //typedef int (&refArray)[SIZE][3][8];
    typedef float (&fullArray1)[SIZE][16];
    typedef float (&fullArray2)[SIZE][16];
    typedef float (&fullArray3)[SIZE][16];
    typedef float (&fullArray4)[SIZE][16];
    typedef float (&fullArray5)[SIZE][16];
    typedef float (&fullArray6)[SIZE][16];
    typedef float (&fullArray7)[SIZE][16];
    typedef float (&fullArray8)[SIZE][16];
    typedef float (&fullArray9)[SIZE][16];
    typedef float (&fullArray10)[SIZE][16];
    typedef float (&fullArray11)[SIZE][16];
    typedef float (&fullArray12)[SIZE][16];
    typedef float (&fullArray13)[SIZE][16];
    typedef float (&fullArray14)[SIZE][16];
    typedef float (&fullArray15)[SIZE][16];
    typedef float (&fullArray16)[SIZE][16];
    typedef float (&fullArray17)[SIZE][16];
    typedef float (&fullArray18)[SIZE][16];
    typedef float (&fullArray19)[SIZE][16];
    typedef float (&fullArray20)[SIZE][16];
    typedef float (&fullArray21)[SIZE][16];
    typedef float (&fullArray22)[SIZE][16];
    typedef float (&fullArray23)[SIZE][16];
    typedef float (&fullArray24)[SIZE][16];
    typedef float (&fullArray25)[SIZE][16];
    typedef float (&fullArray26)[SIZE][16];
    typedef float (&fullArray27)[SIZE][16];
    typedef float (&fullArray28)[SIZE][16];

    // full array extractions
    //refArray _rA = reinterpret_cast < refArray >(*_A);
    fullArray1 _rA = reinterpret_cast < fullArray1 > (*_local_A);
    fullArray2 _rB = reinterpret_cast < fullArray2 > (*_local_B);
    fullArray3 _rC = reinterpret_cast < fullArray3 > (*_local_C);
    fullArray4 _rD = reinterpret_cast < fullArray4 > (*_local_D);
    fullArray5 _rAdd = reinterpret_cast < fullArray5 > (*_local_Add);
    fullArray6 _rMultiply = reinterpret_cast < fullArray6 > (*_local_Multiply);
    fullArray7 _rSubtract = reinterpret_cast < fullArray7 > (*_local_Subtract);
    fullArray8 _rDivide = reinterpret_cast < fullArray8 > (*_local_Divide);
    fullArray9 _rLessThan = reinterpret_cast < fullArray9 > (*_local_LessThan);
    fullArray10 _rGreaterThan =
      reinterpret_cast < fullArray10 > (*_local_GreaterThan);
    fullArray11 _rLessEquals =
      reinterpret_cast < fullArray11 > (*_local_LessEquals);
    fullArray12 _rGreaterEquals =
      reinterpret_cast < fullArray12 > (*_local_GreaterEquals);
    fullArray13 _rEquals = reinterpret_cast < fullArray13 > (*_local_Equals);
    fullArray14 _rAnd = reinterpret_cast < fullArray14 > (*_local_And);
    fullArray15 _rOr = reinterpret_cast < fullArray15 > (*_local_Or);
    fullArray16 _rXor = reinterpret_cast < fullArray16 > (*_local_Xor);
    fullArray17 _rAndNot = reinterpret_cast < fullArray17 > (*_local_AndNot);
    fullArray18 _rNot = reinterpret_cast < fullArray18 > (*_local_Not);
    fullArray19 _rSqrt = reinterpret_cast < fullArray19 > (*_local_Sqrt);
    fullArray20 _rRSqrt = reinterpret_cast < fullArray20 > (*_local_RSqrt);
    fullArray21 _rLog = reinterpret_cast < fullArray21 > (*_local_Log);
    fullArray22 _rExp = reinterpret_cast < fullArray22 > (*_local_Exp);
    fullArray23 _rInverse = reinterpret_cast < fullArray23 > (*_local_Inverse);
    fullArray24 _rAbsoluteValue =
      reinterpret_cast < fullArray24 > (*_local_AbsoluteValue);
    fullArray25 _rSign = reinterpret_cast < fullArray25 > (*_local_Sign);
    fullArray26 _rMinimum = reinterpret_cast < fullArray26 > (*_local_Minimum);
    fullArray27 _rMaximum = reinterpret_cast < fullArray27 > (*_local_Maximum);
    fullArray28 _rBlend = reinterpret_cast < fullArray28 > (*_local_Blend);

    const int ChunkSize = 4;
    // chunk typedef
    //typedef int (&refArray1)[3][8];
    typedef float (&refArray1)[16];
    typedef float (&refArray2)[16];
    typedef float (&refArray3)[16];
    typedef float (&refArray4)[16];
    typedef float (&refArray5)[16];
    typedef float (&refArray6)[16];
    typedef float (&refArray7)[16];
    typedef float (&refArray8)[16];
    typedef float (&refArray9)[16];
    typedef float (&refArray10)[16];
    typedef float (&refArray11)[16];
    typedef float (&refArray12)[16];
    typedef float (&refArray13)[16];
    typedef float (&refArray14)[16];
    typedef float (&refArray15)[16];
    typedef float (&refArray16)[16];
    typedef float (&refArray17)[16];
    typedef float (&refArray18)[16];
    typedef float (&refArray19)[16];
    typedef float (&refArray20)[16];
    typedef float (&refArray21)[16];
    typedef float (&refArray22)[16];
    typedef float (&refArray23)[16];
    typedef float (&refArray24)[16];
    typedef float (&refArray25)[16];
    typedef float (&refArray26)[16];
    typedef float (&refArray27)[16];
    typedef float (&refArray28)[16];

    for (int chunk_offset = 0; chunk_offset < 16; chunk_offset += ChunkSize)
      {
        //refArray1 Ak = reinterpret_cast < refArray1 > (_rA[entry][0][chunk_offset]);
        refArray1 Ak =
          reinterpret_cast < refArray1 > (_rA[index][chunk_offset]);
        refArray2 Bk =
          reinterpret_cast < refArray2 > (_rB[index][chunk_offset]);
        refArray3 Ck =
          reinterpret_cast < refArray3 > (_rC[index][chunk_offset]);
        refArray4 Dk =
          reinterpret_cast < refArray4 > (_rD[index][chunk_offset]);
        refArray5 Addk =
          reinterpret_cast < refArray5 > (_rAdd[index][chunk_offset]);
        refArray6 Multiplyk =
          reinterpret_cast < refArray6 > (_rMultiply[index][chunk_offset]);
        refArray7 Subtractk =
          reinterpret_cast < refArray7 > (_rSubtract[index][chunk_offset]);
        refArray8 Dividek =
          reinterpret_cast < refArray8 > (_rDivide[index][chunk_offset]);
        refArray9 LessThank =
          reinterpret_cast < refArray9 > (_rLessThan[index][chunk_offset]);
        refArray10 GreaterThank =
          reinterpret_cast < refArray10 > (_rGreaterThan[index][chunk_offset]);
        refArray11 LessEqualsk =
          reinterpret_cast < refArray11 > (_rLessEquals[index][chunk_offset]);
        refArray12 GreaterEqualsk =
          reinterpret_cast < refArray12 >
          (_rGreaterEquals[index][chunk_offset]);
        refArray13 Equalsk =
          reinterpret_cast < refArray13 > (_rEquals[index][chunk_offset]);
        refArray14 Andk =
          reinterpret_cast < refArray14 > (_rAnd[index][chunk_offset]);
        refArray15 Ork =
          reinterpret_cast < refArray15 > (_rOr[index][chunk_offset]);
        refArray16 Xork =
          reinterpret_cast < refArray16 > (_rXor[index][chunk_offset]);
        refArray17 AndNotk =
          reinterpret_cast < refArray17 > (_rAndNot[index][chunk_offset]);
        refArray18 Notk =
          reinterpret_cast < refArray18 > (_rNot[index][chunk_offset]);
        refArray19 Sqrtk =
          reinterpret_cast < refArray19 > (_rSqrt[index][chunk_offset]);
        refArray20 RSqrtk =
          reinterpret_cast < refArray20 > (_rRSqrt[index][chunk_offset]);
        refArray21 Logk =
          reinterpret_cast < refArray21 > (_rLog[index][chunk_offset]);
        refArray22 Expk =
          reinterpret_cast < refArray22 > (_rExp[index][chunk_offset]);
        refArray23 Inversek =
          reinterpret_cast < refArray23 > (_rInverse[index][chunk_offset]);
        refArray24 AbsoluteValuek =
          reinterpret_cast < refArray24 >
          (_rAbsoluteValue[index][chunk_offset]);
        refArray25 Signk =
          reinterpret_cast < refArray25 > (_rSign[index][chunk_offset]);
        refArray26 Minimumk =
          reinterpret_cast < refArray26 > (_rMinimum[index][chunk_offset]);
        refArray27 Maximumk =
          reinterpret_cast < refArray27 > (_rMaximum[index][chunk_offset]);
        refArray28 Blendk =
          reinterpret_cast < refArray28 > (_rBlend[index][chunk_offset]);

        NumberTest < __m128, float[16], int[16] > (Ak, Bk, Ck, Dk, Addk,
                                                   Multiplyk, Subtractk,
                                                   Dividek, LessThank,
                                                   GreaterThank, LessEqualsk,
                                                   GreaterEqualsk, Equalsk,
                                                   Andk, Ork, Xork, AndNotk,
                                                   Notk, Sqrtk, RSqrtk, Logk,
                                                   Expk, Inversek,
                                                   AbsoluteValuek, Signk,
                                                   Minimumk, Maximumk, Blendk);
      }

  }
};

#endif

#ifdef ENABLE_AVX_INSTRUCTION_SET

template < int SIZE > class NumberTest_AVX_None
{
private:
  // Generate Variables Here
  float *_local_A;
  float *_local_B;
  float *_local_C;
  float *_local_D;
  float *_local_Add;
  float *_local_Multiply;
  float *_local_Subtract;
  float *_local_Divide;
  float *_local_LessThan;
  float *_local_GreaterThan;
  float *_local_LessEquals;
  float *_local_GreaterEquals;
  float *_local_Equals;
  float *_local_And;
  float *_local_Or;
  float *_local_Xor;
  float *_local_AndNot;
  float *_local_Not;
  float *_local_Sqrt;
  float *_local_RSqrt;
  float *_local_Log;
  float *_local_Exp;
  float *_local_Inverse;
  float *_local_AbsoluteValue;
  float *_local_Sign;
  float *_local_Minimum;
  float *_local_Maximum;
  float *_local_Blend;

public:
    explicit NumberTest_AVX_None (float *A_in, float *B_in, float *C_in,
                                  float *D_in, float *Add_in,
                                  float *Multiply_in, float *Subtract_in,
                                  float *Divide_in, float *LessThan_in,
                                  float *GreaterThan_in, float *LessEquals_in,
                                  float *GreaterEquals_in, float *Equals_in,
                                  float *And_in, float *Or_in, float *Xor_in,
                                  float *AndNot_in, float *Not_in,
                                  float *Sqrt_in, float *RSqrt_in,
                                  float *Log_in, float *Exp_in,
                                  float *Inverse_in, float *AbsoluteValue_in,
                                  float *Sign_in, float *Minimum_in,
                                  float *Maximum_in,
                                  float *Blend_in):_local_A (A_in),
    _local_B (B_in), _local_C (C_in), _local_D (D_in), _local_Add (Add_in),
    _local_Multiply (Multiply_in), _local_Subtract (Subtract_in),
    _local_Divide (Divide_in), _local_LessThan (LessThan_in),
    _local_GreaterThan (GreaterThan_in), _local_LessEquals (LessEquals_in),
    _local_GreaterEquals (GreaterEquals_in), _local_Equals (Equals_in),
    _local_And (And_in), _local_Or (Or_in), _local_Xor (Xor_in),
    _local_AndNot (AndNot_in), _local_Not (Not_in), _local_Sqrt (Sqrt_in),
    _local_RSqrt (RSqrt_in), _local_Log (Log_in), _local_Exp (Exp_in),
    _local_Inverse (Inverse_in), _local_AbsoluteValue (AbsoluteValue_in),
    _local_Sign (Sign_in), _local_Minimum (Minimum_in),
    _local_Maximum (Maximum_in), _local_Blend (Blend_in)
  {
  }
  void Execute (int index)
  {
    // full array typedefs
    //typedef int (&refArray)[SIZE][3][8];
    typedef float (&fullArray1)[SIZE][16];
    typedef float (&fullArray2)[SIZE][16];
    typedef float (&fullArray3)[SIZE][16];
    typedef float (&fullArray4)[SIZE][16];
    typedef float (&fullArray5)[SIZE][16];
    typedef float (&fullArray6)[SIZE][16];
    typedef float (&fullArray7)[SIZE][16];
    typedef float (&fullArray8)[SIZE][16];
    typedef float (&fullArray9)[SIZE][16];
    typedef float (&fullArray10)[SIZE][16];
    typedef float (&fullArray11)[SIZE][16];
    typedef float (&fullArray12)[SIZE][16];
    typedef float (&fullArray13)[SIZE][16];
    typedef float (&fullArray14)[SIZE][16];
    typedef float (&fullArray15)[SIZE][16];
    typedef float (&fullArray16)[SIZE][16];
    typedef float (&fullArray17)[SIZE][16];
    typedef float (&fullArray18)[SIZE][16];
    typedef float (&fullArray19)[SIZE][16];
    typedef float (&fullArray20)[SIZE][16];
    typedef float (&fullArray21)[SIZE][16];
    typedef float (&fullArray22)[SIZE][16];
    typedef float (&fullArray23)[SIZE][16];
    typedef float (&fullArray24)[SIZE][16];
    typedef float (&fullArray25)[SIZE][16];
    typedef float (&fullArray26)[SIZE][16];
    typedef float (&fullArray27)[SIZE][16];
    typedef float (&fullArray28)[SIZE][16];

    // full array extractions
    //refArray _rA = reinterpret_cast < refArray >(*_A);
    fullArray1 _rA = reinterpret_cast < fullArray1 > (*_local_A);
    fullArray2 _rB = reinterpret_cast < fullArray2 > (*_local_B);
    fullArray3 _rC = reinterpret_cast < fullArray3 > (*_local_C);
    fullArray4 _rD = reinterpret_cast < fullArray4 > (*_local_D);
    fullArray5 _rAdd = reinterpret_cast < fullArray5 > (*_local_Add);
    fullArray6 _rMultiply = reinterpret_cast < fullArray6 > (*_local_Multiply);
    fullArray7 _rSubtract = reinterpret_cast < fullArray7 > (*_local_Subtract);
    fullArray8 _rDivide = reinterpret_cast < fullArray8 > (*_local_Divide);
    fullArray9 _rLessThan = reinterpret_cast < fullArray9 > (*_local_LessThan);
    fullArray10 _rGreaterThan =
      reinterpret_cast < fullArray10 > (*_local_GreaterThan);
    fullArray11 _rLessEquals =
      reinterpret_cast < fullArray11 > (*_local_LessEquals);
    fullArray12 _rGreaterEquals =
      reinterpret_cast < fullArray12 > (*_local_GreaterEquals);
    fullArray13 _rEquals = reinterpret_cast < fullArray13 > (*_local_Equals);
    fullArray14 _rAnd = reinterpret_cast < fullArray14 > (*_local_And);
    fullArray15 _rOr = reinterpret_cast < fullArray15 > (*_local_Or);
    fullArray16 _rXor = reinterpret_cast < fullArray16 > (*_local_Xor);
    fullArray17 _rAndNot = reinterpret_cast < fullArray17 > (*_local_AndNot);
    fullArray18 _rNot = reinterpret_cast < fullArray18 > (*_local_Not);
    fullArray19 _rSqrt = reinterpret_cast < fullArray19 > (*_local_Sqrt);
    fullArray20 _rRSqrt = reinterpret_cast < fullArray20 > (*_local_RSqrt);
    fullArray21 _rLog = reinterpret_cast < fullArray21 > (*_local_Log);
    fullArray22 _rExp = reinterpret_cast < fullArray22 > (*_local_Exp);
    fullArray23 _rInverse = reinterpret_cast < fullArray23 > (*_local_Inverse);
    fullArray24 _rAbsoluteValue =
      reinterpret_cast < fullArray24 > (*_local_AbsoluteValue);
    fullArray25 _rSign = reinterpret_cast < fullArray25 > (*_local_Sign);
    fullArray26 _rMinimum = reinterpret_cast < fullArray26 > (*_local_Minimum);
    fullArray27 _rMaximum = reinterpret_cast < fullArray27 > (*_local_Maximum);
    fullArray28 _rBlend = reinterpret_cast < fullArray28 > (*_local_Blend);

    const int ChunkSize = 8;
    // chunk typedef
    //typedef int (&refArray1)[3][8];
    typedef float (&refArray1)[16];
    typedef float (&refArray2)[16];
    typedef float (&refArray3)[16];
    typedef float (&refArray4)[16];
    typedef float (&refArray5)[16];
    typedef float (&refArray6)[16];
    typedef float (&refArray7)[16];
    typedef float (&refArray8)[16];
    typedef float (&refArray9)[16];
    typedef float (&refArray10)[16];
    typedef float (&refArray11)[16];
    typedef float (&refArray12)[16];
    typedef float (&refArray13)[16];
    typedef float (&refArray14)[16];
    typedef float (&refArray15)[16];
    typedef float (&refArray16)[16];
    typedef float (&refArray17)[16];
    typedef float (&refArray18)[16];
    typedef float (&refArray19)[16];
    typedef float (&refArray20)[16];
    typedef float (&refArray21)[16];
    typedef float (&refArray22)[16];
    typedef float (&refArray23)[16];
    typedef float (&refArray24)[16];
    typedef float (&refArray25)[16];
    typedef float (&refArray26)[16];
    typedef float (&refArray27)[16];
    typedef float (&refArray28)[16];

    for (int chunk_offset = 0; chunk_offset < 16; chunk_offset += ChunkSize)
      {
        //refArray1 Ak = reinterpret_cast < refArray1 > (_rA[entry][0][chunk_offset]);
        refArray1 Ak =
          reinterpret_cast < refArray1 > (_rA[index][chunk_offset]);
        refArray2 Bk =
          reinterpret_cast < refArray2 > (_rB[index][chunk_offset]);
        refArray3 Ck =
          reinterpret_cast < refArray3 > (_rC[index][chunk_offset]);
        refArray4 Dk =
          reinterpret_cast < refArray4 > (_rD[index][chunk_offset]);
        refArray5 Addk =
          reinterpret_cast < refArray5 > (_rAdd[index][chunk_offset]);
        refArray6 Multiplyk =
          reinterpret_cast < refArray6 > (_rMultiply[index][chunk_offset]);
        refArray7 Subtractk =
          reinterpret_cast < refArray7 > (_rSubtract[index][chunk_offset]);
        refArray8 Dividek =
          reinterpret_cast < refArray8 > (_rDivide[index][chunk_offset]);
        refArray9 LessThank =
          reinterpret_cast < refArray9 > (_rLessThan[index][chunk_offset]);
        refArray10 GreaterThank =
          reinterpret_cast < refArray10 > (_rGreaterThan[index][chunk_offset]);
        refArray11 LessEqualsk =
          reinterpret_cast < refArray11 > (_rLessEquals[index][chunk_offset]);
        refArray12 GreaterEqualsk =
          reinterpret_cast < refArray12 >
          (_rGreaterEquals[index][chunk_offset]);
        refArray13 Equalsk =
          reinterpret_cast < refArray13 > (_rEquals[index][chunk_offset]);
        refArray14 Andk =
          reinterpret_cast < refArray14 > (_rAnd[index][chunk_offset]);
        refArray15 Ork =
          reinterpret_cast < refArray15 > (_rOr[index][chunk_offset]);
        refArray16 Xork =
          reinterpret_cast < refArray16 > (_rXor[index][chunk_offset]);
        refArray17 AndNotk =
          reinterpret_cast < refArray17 > (_rAndNot[index][chunk_offset]);
        refArray18 Notk =
          reinterpret_cast < refArray18 > (_rNot[index][chunk_offset]);
        refArray19 Sqrtk =
          reinterpret_cast < refArray19 > (_rSqrt[index][chunk_offset]);
        refArray20 RSqrtk =
          reinterpret_cast < refArray20 > (_rRSqrt[index][chunk_offset]);
        refArray21 Logk =
          reinterpret_cast < refArray21 > (_rLog[index][chunk_offset]);
        refArray22 Expk =
          reinterpret_cast < refArray22 > (_rExp[index][chunk_offset]);
        refArray23 Inversek =
          reinterpret_cast < refArray23 > (_rInverse[index][chunk_offset]);
        refArray24 AbsoluteValuek =
          reinterpret_cast < refArray24 >
          (_rAbsoluteValue[index][chunk_offset]);
        refArray25 Signk =
          reinterpret_cast < refArray25 > (_rSign[index][chunk_offset]);
        refArray26 Minimumk =
          reinterpret_cast < refArray26 > (_rMinimum[index][chunk_offset]);
        refArray27 Maximumk =
          reinterpret_cast < refArray27 > (_rMaximum[index][chunk_offset]);
        refArray28 Blendk =
          reinterpret_cast < refArray28 > (_rBlend[index][chunk_offset]);

        NumberTest < __m256, float[16], int[16] > (Ak, Bk, Ck, Dk, Addk,
                                                   Multiplyk, Subtractk,
                                                   Dividek, LessThank,
                                                   GreaterThank, LessEqualsk,
                                                   GreaterEqualsk, Equalsk,
                                                   Andk, Ork, Xork, AndNotk,
                                                   Notk, Sqrtk, RSqrtk, Logk,
                                                   Expk, Inversek,
                                                   AbsoluteValuek, Signk,
                                                   Minimumk, Maximumk, Blendk);
      }

  }
};

#endif

#ifdef ENABLE_NEON_INSTRUCTION_SET

template < int SIZE > class NumberTest_NEON_None
{
private:
  // Generate Variables Here
  float *_local_A;
  float *_local_B;
  float *_local_C;
  float *_local_D;
  float *_local_Add;
  float *_local_Multiply;
  float *_local_Subtract;
  float *_local_Divide;
  float *_local_LessThan;
  float *_local_GreaterThan;
  float *_local_LessEquals;
  float *_local_GreaterEquals;
  float *_local_Equals;
  float *_local_And;
  float *_local_Or;
  float *_local_Xor;
  float *_local_AndNot;
  float *_local_Not;
  float *_local_Sqrt;
  float *_local_RSqrt;
  float *_local_Log;
  float *_local_Exp;
  float *_local_Inverse;
  float *_local_AbsoluteValue;
  float *_local_Sign;
  float *_local_Minimum;
  float *_local_Maximum;
  float *_local_Blend;

public:
    explicit NumberTest_NEON_None (float *A_in, float *B_in, float *C_in,
                                   float *D_in, float *Add_in,
                                   float *Multiply_in, float *Subtract_in,
                                   float *Divide_in, float *LessThan_in,
                                   float *GreaterThan_in, float *LessEquals_in,
                                   float *GreaterEquals_in, float *Equals_in,
                                   float *And_in, float *Or_in, float *Xor_in,
                                   float *AndNot_in, float *Not_in,
                                   float *Sqrt_in, float *RSqrt_in,
                                   float *Log_in, float *Exp_in,
                                   float *Inverse_in, float *AbsoluteValue_in,
                                   float *Sign_in, float *Minimum_in,
                                   float *Maximum_in,
                                   float *Blend_in):_local_A (A_in),
    _local_B (B_in), _local_C (C_in), _local_D (D_in), _local_Add (Add_in),
    _local_Multiply (Multiply_in), _local_Subtract (Subtract_in),
    _local_Divide (Divide_in), _local_LessThan (LessThan_in),
    _local_GreaterThan (GreaterThan_in), _local_LessEquals (LessEquals_in),
    _local_GreaterEquals (GreaterEquals_in), _local_Equals (Equals_in),
    _local_And (And_in), _local_Or (Or_in), _local_Xor (Xor_in),
    _local_AndNot (AndNot_in), _local_Not (Not_in), _local_Sqrt (Sqrt_in),
    _local_RSqrt (RSqrt_in), _local_Log (Log_in), _local_Exp (Exp_in),
    _local_Inverse (Inverse_in), _local_AbsoluteValue (AbsoluteValue_in),
    _local_Sign (Sign_in), _local_Minimum (Minimum_in),
    _local_Maximum (Maximum_in), _local_Blend (Blend_in)
  {
  }
  void Execute (int index)
  {
    // full array typedefs
    //typedef int (&refArray)[SIZE][3][8];
    typedef float (&fullArray1)[SIZE][16];
    typedef float (&fullArray2)[SIZE][16];
    typedef float (&fullArray3)[SIZE][16];
    typedef float (&fullArray4)[SIZE][16];
    typedef float (&fullArray5)[SIZE][16];
    typedef float (&fullArray6)[SIZE][16];
    typedef float (&fullArray7)[SIZE][16];
    typedef float (&fullArray8)[SIZE][16];
    typedef float (&fullArray9)[SIZE][16];
    typedef float (&fullArray10)[SIZE][16];
    typedef float (&fullArray11)[SIZE][16];
    typedef float (&fullArray12)[SIZE][16];
    typedef float (&fullArray13)[SIZE][16];
    typedef float (&fullArray14)[SIZE][16];
    typedef float (&fullArray15)[SIZE][16];
    typedef float (&fullArray16)[SIZE][16];
    typedef float (&fullArray17)[SIZE][16];
    typedef float (&fullArray18)[SIZE][16];
    typedef float (&fullArray19)[SIZE][16];
    typedef float (&fullArray20)[SIZE][16];
    typedef float (&fullArray21)[SIZE][16];
    typedef float (&fullArray22)[SIZE][16];
    typedef float (&fullArray23)[SIZE][16];
    typedef float (&fullArray24)[SIZE][16];
    typedef float (&fullArray25)[SIZE][16];
    typedef float (&fullArray26)[SIZE][16];
    typedef float (&fullArray27)[SIZE][16];
    typedef float (&fullArray28)[SIZE][16];

    // full array extractions
    //refArray _rA = reinterpret_cast < refArray >(*_A);
    fullArray1 _rA = reinterpret_cast < fullArray1 > (*_local_A);
    fullArray2 _rB = reinterpret_cast < fullArray2 > (*_local_B);
    fullArray3 _rC = reinterpret_cast < fullArray3 > (*_local_C);
    fullArray4 _rD = reinterpret_cast < fullArray4 > (*_local_D);
    fullArray5 _rAdd = reinterpret_cast < fullArray5 > (*_local_Add);
    fullArray6 _rMultiply = reinterpret_cast < fullArray6 > (*_local_Multiply);
    fullArray7 _rSubtract = reinterpret_cast < fullArray7 > (*_local_Subtract);
    fullArray8 _rDivide = reinterpret_cast < fullArray8 > (*_local_Divide);
    fullArray9 _rLessThan = reinterpret_cast < fullArray9 > (*_local_LessThan);
    fullArray10 _rGreaterThan =
      reinterpret_cast < fullArray10 > (*_local_GreaterThan);
    fullArray11 _rLessEquals =
      reinterpret_cast < fullArray11 > (*_local_LessEquals);
    fullArray12 _rGreaterEquals =
      reinterpret_cast < fullArray12 > (*_local_GreaterEquals);
    fullArray13 _rEquals = reinterpret_cast < fullArray13 > (*_local_Equals);
    fullArray14 _rAnd = reinterpret_cast < fullArray14 > (*_local_And);
    fullArray15 _rOr = reinterpret_cast < fullArray15 > (*_local_Or);
    fullArray16 _rXor = reinterpret_cast < fullArray16 > (*_local_Xor);
    fullArray17 _rAndNot = reinterpret_cast < fullArray17 > (*_local_AndNot);
    fullArray18 _rNot = reinterpret_cast < fullArray18 > (*_local_Not);
    fullArray19 _rSqrt = reinterpret_cast < fullArray19 > (*_local_Sqrt);
    fullArray20 _rRSqrt = reinterpret_cast < fullArray20 > (*_local_RSqrt);
    fullArray21 _rLog = reinterpret_cast < fullArray21 > (*_local_Log);
    fullArray22 _rExp = reinterpret_cast < fullArray22 > (*_local_Exp);
    fullArray23 _rInverse = reinterpret_cast < fullArray23 > (*_local_Inverse);
    fullArray24 _rAbsoluteValue =
      reinterpret_cast < fullArray24 > (*_local_AbsoluteValue);
    fullArray25 _rSign = reinterpret_cast < fullArray25 > (*_local_Sign);
    fullArray26 _rMinimum = reinterpret_cast < fullArray26 > (*_local_Minimum);
    fullArray27 _rMaximum = reinterpret_cast < fullArray27 > (*_local_Maximum);
    fullArray28 _rBlend = reinterpret_cast < fullArray28 > (*_local_Blend);

    const int ChunkSize = 4;
    // chunk typedef
    //typedef int (&refArray1)[3][8];
    typedef float (&refArray1)[16];
    typedef float (&refArray2)[16];
    typedef float (&refArray3)[16];
    typedef float (&refArray4)[16];
    typedef float (&refArray5)[16];
    typedef float (&refArray6)[16];
    typedef float (&refArray7)[16];
    typedef float (&refArray8)[16];
    typedef float (&refArray9)[16];
    typedef float (&refArray10)[16];
    typedef float (&refArray11)[16];
    typedef float (&refArray12)[16];
    typedef float (&refArray13)[16];
    typedef float (&refArray14)[16];
    typedef float (&refArray15)[16];
    typedef float (&refArray16)[16];
    typedef float (&refArray17)[16];
    typedef float (&refArray18)[16];
    typedef float (&refArray19)[16];
    typedef float (&refArray20)[16];
    typedef float (&refArray21)[16];
    typedef float (&refArray22)[16];
    typedef float (&refArray23)[16];
    typedef float (&refArray24)[16];
    typedef float (&refArray25)[16];
    typedef float (&refArray26)[16];
    typedef float (&refArray27)[16];
    typedef float (&refArray28)[16];

    for (int chunk_offset = 0; chunk_offset < 16; chunk_offset += ChunkSize)
      {
        //refArray1 Ak = reinterpret_cast < refArray1 > (_rA[entry][0][chunk_offset]);
        refArray1 Ak =
          reinterpret_cast < refArray1 > (_rA[index][chunk_offset]);
        refArray2 Bk =
          reinterpret_cast < refArray2 > (_rB[index][chunk_offset]);
        refArray3 Ck =
          reinterpret_cast < refArray3 > (_rC[index][chunk_offset]);
        refArray4 Dk =
          reinterpret_cast < refArray4 > (_rD[index][chunk_offset]);
        refArray5 Addk =
          reinterpret_cast < refArray5 > (_rAdd[index][chunk_offset]);
        refArray6 Multiplyk =
          reinterpret_cast < refArray6 > (_rMultiply[index][chunk_offset]);
        refArray7 Subtractk =
          reinterpret_cast < refArray7 > (_rSubtract[index][chunk_offset]);
        refArray8 Dividek =
          reinterpret_cast < refArray8 > (_rDivide[index][chunk_offset]);
        refArray9 LessThank =
          reinterpret_cast < refArray9 > (_rLessThan[index][chunk_offset]);
        refArray10 GreaterThank =
          reinterpret_cast < refArray10 > (_rGreaterThan[index][chunk_offset]);
        refArray11 LessEqualsk =
          reinterpret_cast < refArray11 > (_rLessEquals[index][chunk_offset]);
        refArray12 GreaterEqualsk =
          reinterpret_cast < refArray12 >
          (_rGreaterEquals[index][chunk_offset]);
        refArray13 Equalsk =
          reinterpret_cast < refArray13 > (_rEquals[index][chunk_offset]);
        refArray14 Andk =
          reinterpret_cast < refArray14 > (_rAnd[index][chunk_offset]);
        refArray15 Ork =
          reinterpret_cast < refArray15 > (_rOr[index][chunk_offset]);
        refArray16 Xork =
          reinterpret_cast < refArray16 > (_rXor[index][chunk_offset]);
        refArray17 AndNotk =
          reinterpret_cast < refArray17 > (_rAndNot[index][chunk_offset]);
        refArray18 Notk =
          reinterpret_cast < refArray18 > (_rNot[index][chunk_offset]);
        refArray19 Sqrtk =
          reinterpret_cast < refArray19 > (_rSqrt[index][chunk_offset]);
        refArray20 RSqrtk =
          reinterpret_cast < refArray20 > (_rRSqrt[index][chunk_offset]);
        refArray21 Logk =
          reinterpret_cast < refArray21 > (_rLog[index][chunk_offset]);
        refArray22 Expk =
          reinterpret_cast < refArray22 > (_rExp[index][chunk_offset]);
        refArray23 Inversek =
          reinterpret_cast < refArray23 > (_rInverse[index][chunk_offset]);
        refArray24 AbsoluteValuek =
          reinterpret_cast < refArray24 >
          (_rAbsoluteValue[index][chunk_offset]);
        refArray25 Signk =
          reinterpret_cast < refArray25 > (_rSign[index][chunk_offset]);
        refArray26 Minimumk =
          reinterpret_cast < refArray26 > (_rMinimum[index][chunk_offset]);
        refArray27 Maximumk =
          reinterpret_cast < refArray27 > (_rMaximum[index][chunk_offset]);
        refArray28 Blendk =
          reinterpret_cast < refArray28 > (_rBlend[index][chunk_offset]);

        NumberTest < float32x4_t, float[16], int[16] > (Ak, Bk, Ck, Dk, Addk,
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
};

#endif

#ifdef ENABLE_MIC_INSTRUCTION_SET

template < int SIZE > class NumberTest_MIC_None
{
private:
  // Generate Variables Here
  float *_local_A;
  float *_local_B;
  float *_local_C;
  float *_local_D;
  float *_local_Add;
  float *_local_Multiply;
  float *_local_Subtract;
  float *_local_Divide;
  float *_local_LessThan;
  float *_local_GreaterThan;
  float *_local_LessEquals;
  float *_local_GreaterEquals;
  float *_local_Equals;
  float *_local_And;
  float *_local_Or;
  float *_local_Xor;
  float *_local_AndNot;
  float *_local_Not;
  float *_local_Sqrt;
  float *_local_RSqrt;
  float *_local_Log;
  float *_local_Exp;
  float *_local_Inverse;
  float *_local_AbsoluteValue;
  float *_local_Sign;
  float *_local_Minimum;
  float *_local_Maximum;
  float *_local_Blend;

public:
    explicit NumberTest_MIC_None (float *A_in, float *B_in, float *C_in,
                                  float *D_in, float *Add_in,
                                  float *Multiply_in, float *Subtract_in,
                                  float *Divide_in, float *LessThan_in,
                                  float *GreaterThan_in, float *LessEquals_in,
                                  float *GreaterEquals_in, float *Equals_in,
                                  float *And_in, float *Or_in, float *Xor_in,
                                  float *AndNot_in, float *Not_in,
                                  float *Sqrt_in, float *RSqrt_in,
                                  float *Log_in, float *Exp_in,
                                  float *Inverse_in, float *AbsoluteValue_in,
                                  float *Sign_in, float *Minimum_in,
                                  float *Maximum_in,
                                  float *Blend_in):_local_A (A_in),
    _local_B (B_in), _local_C (C_in), _local_D (D_in), _local_Add (Add_in),
    _local_Multiply (Multiply_in), _local_Subtract (Subtract_in),
    _local_Divide (Divide_in), _local_LessThan (LessThan_in),
    _local_GreaterThan (GreaterThan_in), _local_LessEquals (LessEquals_in),
    _local_GreaterEquals (GreaterEquals_in), _local_Equals (Equals_in),
    _local_And (And_in), _local_Or (Or_in), _local_Xor (Xor_in),
    _local_AndNot (AndNot_in), _local_Not (Not_in), _local_Sqrt (Sqrt_in),
    _local_RSqrt (RSqrt_in), _local_Log (Log_in), _local_Exp (Exp_in),
    _local_Inverse (Inverse_in), _local_AbsoluteValue (AbsoluteValue_in),
    _local_Sign (Sign_in), _local_Minimum (Minimum_in),
    _local_Maximum (Maximum_in), _local_Blend (Blend_in)
  {
  }
  void Execute (int index)
  {
    // full array typedefs
    //typedef int (&refArray)[SIZE][3][8];
    typedef float (&fullArray1)[SIZE][16];
    typedef float (&fullArray2)[SIZE][16];
    typedef float (&fullArray3)[SIZE][16];
    typedef float (&fullArray4)[SIZE][16];
    typedef float (&fullArray5)[SIZE][16];
    typedef float (&fullArray6)[SIZE][16];
    typedef float (&fullArray7)[SIZE][16];
    typedef float (&fullArray8)[SIZE][16];
    typedef float (&fullArray9)[SIZE][16];
    typedef float (&fullArray10)[SIZE][16];
    typedef float (&fullArray11)[SIZE][16];
    typedef float (&fullArray12)[SIZE][16];
    typedef float (&fullArray13)[SIZE][16];
    typedef float (&fullArray14)[SIZE][16];
    typedef float (&fullArray15)[SIZE][16];
    typedef float (&fullArray16)[SIZE][16];
    typedef float (&fullArray17)[SIZE][16];
    typedef float (&fullArray18)[SIZE][16];
    typedef float (&fullArray19)[SIZE][16];
    typedef float (&fullArray20)[SIZE][16];
    typedef float (&fullArray21)[SIZE][16];
    typedef float (&fullArray22)[SIZE][16];
    typedef float (&fullArray23)[SIZE][16];
    typedef float (&fullArray24)[SIZE][16];
    typedef float (&fullArray25)[SIZE][16];
    typedef float (&fullArray26)[SIZE][16];
    typedef float (&fullArray27)[SIZE][16];
    typedef float (&fullArray28)[SIZE][16];

    // full array extractions
    //refArray _rA = reinterpret_cast < refArray >(*_A);
    fullArray1 _rA = reinterpret_cast < fullArray1 > (*_local_A);
    fullArray2 _rB = reinterpret_cast < fullArray2 > (*_local_B);
    fullArray3 _rC = reinterpret_cast < fullArray3 > (*_local_C);
    fullArray4 _rD = reinterpret_cast < fullArray4 > (*_local_D);
    fullArray5 _rAdd = reinterpret_cast < fullArray5 > (*_local_Add);
    fullArray6 _rMultiply = reinterpret_cast < fullArray6 > (*_local_Multiply);
    fullArray7 _rSubtract = reinterpret_cast < fullArray7 > (*_local_Subtract);
    fullArray8 _rDivide = reinterpret_cast < fullArray8 > (*_local_Divide);
    fullArray9 _rLessThan = reinterpret_cast < fullArray9 > (*_local_LessThan);
    fullArray10 _rGreaterThan =
      reinterpret_cast < fullArray10 > (*_local_GreaterThan);
    fullArray11 _rLessEquals =
      reinterpret_cast < fullArray11 > (*_local_LessEquals);
    fullArray12 _rGreaterEquals =
      reinterpret_cast < fullArray12 > (*_local_GreaterEquals);
    fullArray13 _rEquals = reinterpret_cast < fullArray13 > (*_local_Equals);
    fullArray14 _rAnd = reinterpret_cast < fullArray14 > (*_local_And);
    fullArray15 _rOr = reinterpret_cast < fullArray15 > (*_local_Or);
    fullArray16 _rXor = reinterpret_cast < fullArray16 > (*_local_Xor);
    fullArray17 _rAndNot = reinterpret_cast < fullArray17 > (*_local_AndNot);
    fullArray18 _rNot = reinterpret_cast < fullArray18 > (*_local_Not);
    fullArray19 _rSqrt = reinterpret_cast < fullArray19 > (*_local_Sqrt);
    fullArray20 _rRSqrt = reinterpret_cast < fullArray20 > (*_local_RSqrt);
    fullArray21 _rLog = reinterpret_cast < fullArray21 > (*_local_Log);
    fullArray22 _rExp = reinterpret_cast < fullArray22 > (*_local_Exp);
    fullArray23 _rInverse = reinterpret_cast < fullArray23 > (*_local_Inverse);
    fullArray24 _rAbsoluteValue =
      reinterpret_cast < fullArray24 > (*_local_AbsoluteValue);
    fullArray25 _rSign = reinterpret_cast < fullArray25 > (*_local_Sign);
    fullArray26 _rMinimum = reinterpret_cast < fullArray26 > (*_local_Minimum);
    fullArray27 _rMaximum = reinterpret_cast < fullArray27 > (*_local_Maximum);
    fullArray28 _rBlend = reinterpret_cast < fullArray28 > (*_local_Blend);

    const int ChunkSize = 16;
    // chunk typedef
    //typedef int (&refArray1)[3][8];
    typedef float (&refArray1)[16];
    typedef float (&refArray2)[16];
    typedef float (&refArray3)[16];
    typedef float (&refArray4)[16];
    typedef float (&refArray5)[16];
    typedef float (&refArray6)[16];
    typedef float (&refArray7)[16];
    typedef float (&refArray8)[16];
    typedef float (&refArray9)[16];
    typedef float (&refArray10)[16];
    typedef float (&refArray11)[16];
    typedef float (&refArray12)[16];
    typedef float (&refArray13)[16];
    typedef float (&refArray14)[16];
    typedef float (&refArray15)[16];
    typedef float (&refArray16)[16];
    typedef float (&refArray17)[16];
    typedef float (&refArray18)[16];
    typedef float (&refArray19)[16];
    typedef float (&refArray20)[16];
    typedef float (&refArray21)[16];
    typedef float (&refArray22)[16];
    typedef float (&refArray23)[16];
    typedef float (&refArray24)[16];
    typedef float (&refArray25)[16];
    typedef float (&refArray26)[16];
    typedef float (&refArray27)[16];
    typedef float (&refArray28)[16];

    for (int chunk_offset = 0; chunk_offset < 16; chunk_offset += ChunkSize)
      {
        //refArray1 Ak = reinterpret_cast < refArray1 > (_rA[entry][0][chunk_offset]);
        refArray1 Ak =
          reinterpret_cast < refArray1 > (_rA[index][chunk_offset]);
        refArray2 Bk =
          reinterpret_cast < refArray2 > (_rB[index][chunk_offset]);
        refArray3 Ck =
          reinterpret_cast < refArray3 > (_rC[index][chunk_offset]);
        refArray4 Dk =
          reinterpret_cast < refArray4 > (_rD[index][chunk_offset]);
        refArray5 Addk =
          reinterpret_cast < refArray5 > (_rAdd[index][chunk_offset]);
        refArray6 Multiplyk =
          reinterpret_cast < refArray6 > (_rMultiply[index][chunk_offset]);
        refArray7 Subtractk =
          reinterpret_cast < refArray7 > (_rSubtract[index][chunk_offset]);
        refArray8 Dividek =
          reinterpret_cast < refArray8 > (_rDivide[index][chunk_offset]);
        refArray9 LessThank =
          reinterpret_cast < refArray9 > (_rLessThan[index][chunk_offset]);
        refArray10 GreaterThank =
          reinterpret_cast < refArray10 > (_rGreaterThan[index][chunk_offset]);
        refArray11 LessEqualsk =
          reinterpret_cast < refArray11 > (_rLessEquals[index][chunk_offset]);
        refArray12 GreaterEqualsk =
          reinterpret_cast < refArray12 >
          (_rGreaterEquals[index][chunk_offset]);
        refArray13 Equalsk =
          reinterpret_cast < refArray13 > (_rEquals[index][chunk_offset]);
        refArray14 Andk =
          reinterpret_cast < refArray14 > (_rAnd[index][chunk_offset]);
        refArray15 Ork =
          reinterpret_cast < refArray15 > (_rOr[index][chunk_offset]);
        refArray16 Xork =
          reinterpret_cast < refArray16 > (_rXor[index][chunk_offset]);
        refArray17 AndNotk =
          reinterpret_cast < refArray17 > (_rAndNot[index][chunk_offset]);
        refArray18 Notk =
          reinterpret_cast < refArray18 > (_rNot[index][chunk_offset]);
        refArray19 Sqrtk =
          reinterpret_cast < refArray19 > (_rSqrt[index][chunk_offset]);
        refArray20 RSqrtk =
          reinterpret_cast < refArray20 > (_rRSqrt[index][chunk_offset]);
        refArray21 Logk =
          reinterpret_cast < refArray21 > (_rLog[index][chunk_offset]);
        refArray22 Expk =
          reinterpret_cast < refArray22 > (_rExp[index][chunk_offset]);
        refArray23 Inversek =
          reinterpret_cast < refArray23 > (_rInverse[index][chunk_offset]);
        refArray24 AbsoluteValuek =
          reinterpret_cast < refArray24 >
          (_rAbsoluteValue[index][chunk_offset]);
        refArray25 Signk =
          reinterpret_cast < refArray25 > (_rSign[index][chunk_offset]);
        refArray26 Minimumk =
          reinterpret_cast < refArray26 > (_rMinimum[index][chunk_offset]);
        refArray27 Maximumk =
          reinterpret_cast < refArray27 > (_rMaximum[index][chunk_offset]);
        refArray28 Blendk =
          reinterpret_cast < refArray28 > (_rBlend[index][chunk_offset]);

        NumberTest < __m512, float[16], int[16] > (Ak, Bk, Ck, Dk, Addk,
                                                   Multiplyk, Subtractk,
                                                   Dividek, LessThank,
                                                   GreaterThank, LessEqualsk,
                                                   GreaterEqualsk, Equalsk,
                                                   Andk, Ork, Xork, AndNotk,
                                                   Notk, Sqrtk, RSqrtk, Logk,
                                                   Expk, Inversek,
                                                   AbsoluteValuek, Signk,
                                                   Minimumk, Maximumk, Blendk);
      }

  }
};

#endif

int
main (int argc, char *argv[])
{
  typedef float T;

  int seed = 1;
  int threads = 1;
  int threads_max = 1;
  int passes = 1;
  const int data_size = 1000000;
  if (argc >= 2)
    threads = atoi (argv[1]);
  if (argc >= 3)
    threads_max = atoi (argv[2]);
  if (argc >= 4)
    passes = atoi (argv[3]);
  srand (seed);

  pthread_queue = new PTHREAD_QUEUE (threads_max);

  std::
    cout << "Preparing to Run " << data_size << " of all kernels with " <<
    threads << " threads." << std::endl;



  {
    std::cout << "Running Thread Test for NumberTest " << std::endl;

//=======================================================
//
//        DEFINE ALL VARIABLES USED BY KERNEL
//
//=======================================================
    std::cout << "\nAllocating all data: ";
    std::cout.flush ();

    start_timer ();
    typedef T (&A_type)[data_size][16];
    A_type A =
      reinterpret_cast < A_type >
      (*((T *) (_mm_malloc (data_size * 16 * sizeof (T), 64))));
    typedef T (&B_type)[data_size][16];
    B_type B =
      reinterpret_cast < B_type >
      (*((T *) (_mm_malloc (data_size * 16 * sizeof (T), 64))));
    typedef T (&C_type)[data_size][16];
    C_type C =
      reinterpret_cast < C_type >
      (*((T *) (_mm_malloc (data_size * 16 * sizeof (T), 64))));
    typedef T (&D_type)[data_size][16];
    D_type D =
      reinterpret_cast < D_type >
      (*((T *) (_mm_malloc (data_size * 16 * sizeof (T), 64))));
    typedef T (&Add_type)[data_size][16];
    Add_type Add =
      reinterpret_cast < Add_type >
      (*((T *) (_mm_malloc (data_size * 16 * sizeof (T), 64))));
    typedef T (&Multiply_type)[data_size][16];
    Multiply_type Multiply =
      reinterpret_cast < Multiply_type >
      (*((T *) (_mm_malloc (data_size * 16 * sizeof (T), 64))));
    typedef T (&Subtract_type)[data_size][16];
    Subtract_type Subtract =
      reinterpret_cast < Subtract_type >
      (*((T *) (_mm_malloc (data_size * 16 * sizeof (T), 64))));
    typedef T (&Divide_type)[data_size][16];
    Divide_type Divide =
      reinterpret_cast < Divide_type >
      (*((T *) (_mm_malloc (data_size * 16 * sizeof (T), 64))));
    typedef T (&LessThan_type)[data_size][16];
    LessThan_type LessThan =
      reinterpret_cast < LessThan_type >
      (*((T *) (_mm_malloc (data_size * 16 * sizeof (T), 64))));
    typedef T (&GreaterThan_type)[data_size][16];
    GreaterThan_type GreaterThan =
      reinterpret_cast < GreaterThan_type >
      (*((T *) (_mm_malloc (data_size * 16 * sizeof (T), 64))));
    typedef T (&LessEquals_type)[data_size][16];
    LessEquals_type LessEquals =
      reinterpret_cast < LessEquals_type >
      (*((T *) (_mm_malloc (data_size * 16 * sizeof (T), 64))));
    typedef T (&GreaterEquals_type)[data_size][16];
    GreaterEquals_type GreaterEquals =
      reinterpret_cast < GreaterEquals_type >
      (*((T *) (_mm_malloc (data_size * 16 * sizeof (T), 64))));
    typedef T (&Equals_type)[data_size][16];
    Equals_type Equals =
      reinterpret_cast < Equals_type >
      (*((T *) (_mm_malloc (data_size * 16 * sizeof (T), 64))));
    typedef T (&And_type)[data_size][16];
    And_type And =
      reinterpret_cast < And_type >
      (*((T *) (_mm_malloc (data_size * 16 * sizeof (T), 64))));
    typedef T (&Or_type)[data_size][16];
    Or_type Or =
      reinterpret_cast < Or_type >
      (*((T *) (_mm_malloc (data_size * 16 * sizeof (T), 64))));
    typedef T (&Xor_type)[data_size][16];
    Xor_type Xor =
      reinterpret_cast < Xor_type >
      (*((T *) (_mm_malloc (data_size * 16 * sizeof (T), 64))));
    typedef T (&AndNot_type)[data_size][16];
    AndNot_type AndNot =
      reinterpret_cast < AndNot_type >
      (*((T *) (_mm_malloc (data_size * 16 * sizeof (T), 64))));
    typedef T (&Not_type)[data_size][16];
    Not_type Not =
      reinterpret_cast < Not_type >
      (*((T *) (_mm_malloc (data_size * 16 * sizeof (T), 64))));
    typedef T (&Sqrt_type)[data_size][16];
    Sqrt_type Sqrt =
      reinterpret_cast < Sqrt_type >
      (*((T *) (_mm_malloc (data_size * 16 * sizeof (T), 64))));
    typedef T (&RSqrt_type)[data_size][16];
    RSqrt_type RSqrt =
      reinterpret_cast < RSqrt_type >
      (*((T *) (_mm_malloc (data_size * 16 * sizeof (T), 64))));
    typedef T (&Log_type)[data_size][16];
    Log_type Log =
      reinterpret_cast < Log_type >
      (*((T *) (_mm_malloc (data_size * 16 * sizeof (T), 64))));
    typedef T (&Exp_type)[data_size][16];
    Exp_type Exp =
      reinterpret_cast < Exp_type >
      (*((T *) (_mm_malloc (data_size * 16 * sizeof (T), 64))));
    typedef T (&Inverse_type)[data_size][16];
    Inverse_type Inverse =
      reinterpret_cast < Inverse_type >
      (*((T *) (_mm_malloc (data_size * 16 * sizeof (T), 64))));
    typedef T (&AbsoluteValue_type)[data_size][16];
    AbsoluteValue_type AbsoluteValue =
      reinterpret_cast < AbsoluteValue_type >
      (*((T *) (_mm_malloc (data_size * 16 * sizeof (T), 64))));
    typedef T (&Sign_type)[data_size][16];
    Sign_type Sign =
      reinterpret_cast < Sign_type >
      (*((T *) (_mm_malloc (data_size * 16 * sizeof (T), 64))));
    typedef T (&Minimum_type)[data_size][16];
    Minimum_type Minimum =
      reinterpret_cast < Minimum_type >
      (*((T *) (_mm_malloc (data_size * 16 * sizeof (T), 64))));
    typedef T (&Maximum_type)[data_size][16];
    Maximum_type Maximum =
      reinterpret_cast < Maximum_type >
      (*((T *) (_mm_malloc (data_size * 16 * sizeof (T), 64))));
    typedef T (&Blend_type)[data_size][16];
    Blend_type Blend =
      reinterpret_cast < Blend_type >
      (*((T *) (_mm_malloc (data_size * 16 * sizeof (T), 64))));


    for (int __a = 0; __a < data_size; __a++)
      for (int __b = 0; __b < 16; __b++)
        {
          A[__a][__b] = Get_Random < float >();
    } for (int __a = 0; __a < data_size; __a++)
      for (int __b = 0; __b < 16; __b++)
        {
          B[__a][__b] = Get_Random < float >();
    } for (int __a = 0; __a < data_size; __a++)
      for (int __b = 0; __b < 16; __b++)
        {
          C[__a][__b] = Get_Random < float >();
    } for (int __a = 0; __a < data_size; __a++)
      for (int __b = 0; __b < 16; __b++)
        {
          D[__a][__b] = Get_Random < float >();
    } for (int __a = 0; __a < data_size; __a++)
      for (int __b = 0; __b < 16; __b++)
        {
          Add[__a][__b] = Get_Random < float >();
    } for (int __a = 0; __a < data_size; __a++)
      for (int __b = 0; __b < 16; __b++)
        {
          Multiply[__a][__b] = Get_Random < float >();
    } for (int __a = 0; __a < data_size; __a++)
      for (int __b = 0; __b < 16; __b++)
        {
          Subtract[__a][__b] = Get_Random < float >();
    } for (int __a = 0; __a < data_size; __a++)
      for (int __b = 0; __b < 16; __b++)
        {
          Divide[__a][__b] = Get_Random < float >();
    } for (int __a = 0; __a < data_size; __a++)
      for (int __b = 0; __b < 16; __b++)
        {
          LessThan[__a][__b] = Get_Random < float >();
    } for (int __a = 0; __a < data_size; __a++)
      for (int __b = 0; __b < 16; __b++)
        {
          GreaterThan[__a][__b] = Get_Random < float >();
    } for (int __a = 0; __a < data_size; __a++)
      for (int __b = 0; __b < 16; __b++)
        {
          LessEquals[__a][__b] = Get_Random < float >();
    } for (int __a = 0; __a < data_size; __a++)
      for (int __b = 0; __b < 16; __b++)
        {
          GreaterEquals[__a][__b] = Get_Random < float >();
    } for (int __a = 0; __a < data_size; __a++)
      for (int __b = 0; __b < 16; __b++)
        {
          Equals[__a][__b] = Get_Random < float >();
    } for (int __a = 0; __a < data_size; __a++)
      for (int __b = 0; __b < 16; __b++)
        {
          And[__a][__b] = Get_Random < float >();
    } for (int __a = 0; __a < data_size; __a++)
      for (int __b = 0; __b < 16; __b++)
        {
          Or[__a][__b] = Get_Random < float >();
    } for (int __a = 0; __a < data_size; __a++)
      for (int __b = 0; __b < 16; __b++)
        {
          Xor[__a][__b] = Get_Random < float >();
    } for (int __a = 0; __a < data_size; __a++)
      for (int __b = 0; __b < 16; __b++)
        {
          AndNot[__a][__b] = Get_Random < float >();
    } for (int __a = 0; __a < data_size; __a++)
      for (int __b = 0; __b < 16; __b++)
        {
          Not[__a][__b] = Get_Random < float >();
    } for (int __a = 0; __a < data_size; __a++)
      for (int __b = 0; __b < 16; __b++)
        {
          Sqrt[__a][__b] = Get_Random < float >();
    } for (int __a = 0; __a < data_size; __a++)
      for (int __b = 0; __b < 16; __b++)
        {
          RSqrt[__a][__b] = Get_Random < float >();
    } for (int __a = 0; __a < data_size; __a++)
      for (int __b = 0; __b < 16; __b++)
        {
          Log[__a][__b] = Get_Random < float >();
    } for (int __a = 0; __a < data_size; __a++)
      for (int __b = 0; __b < 16; __b++)
        {
          Exp[__a][__b] = Get_Random < float >();
    } for (int __a = 0; __a < data_size; __a++)
      for (int __b = 0; __b < 16; __b++)
        {
          Inverse[__a][__b] = Get_Random < float >();
    } for (int __a = 0; __a < data_size; __a++)
      for (int __b = 0; __b < 16; __b++)
        {
          AbsoluteValue[__a][__b] = Get_Random < float >();
    } for (int __a = 0; __a < data_size; __a++)
      for (int __b = 0; __b < 16; __b++)
        {
          Sign[__a][__b] = Get_Random < float >();
    } for (int __a = 0; __a < data_size; __a++)
      for (int __b = 0; __b < 16; __b++)
        {
          Minimum[__a][__b] = Get_Random < float >();
    } for (int __a = 0; __a < data_size; __a++)
      for (int __b = 0; __b < 16; __b++)
        {
          Maximum[__a][__b] = Get_Random < float >();
    } for (int __a = 0; __a < data_size; __a++)
      for (int __b = 0; __b < 16; __b++)
        {
          Blend[__a][__b] = Get_Random < float >();
        }
    stop_timer ();

    std::cout << get_time () << "s\n\n" << std::endl;


//=======================================================
//
//             COMPUTE SCALAR RESULTS
//
//=======================================================

    {
      std::cout << "	Running " << data_size << " of SCALAR :  " << std::endl;


      NumberTest_SCALAR_None < data_size > op ((float *) &A, (float *) &B,
                                               (float *) &C, (float *) &D,
                                               (float *) &Add,
                                               (float *) &Multiply,
                                               (float *) &Subtract,
                                               (float *) &Divide,
                                               (float *) &LessThan,
                                               (float *) &GreaterThan,
                                               (float *) &LessEquals,
                                               (float *) &GreaterEquals,
                                               (float *) &Equals,
                                               (float *) &And, (float *) &Or,
                                               (float *) &Xor,
                                               (float *) &AndNot,
                                               (float *) &Not, (float *) &Sqrt,
                                               (float *) &RSqrt, (float *) &Log,
                                               (float *) &Exp,
                                               (float *) &Inverse,
                                               (float *) &AbsoluteValue,
                                               (float *) &Sign,
                                               (float *) &Minimum,
                                               (float *) &Maximum,
                                               (float *) &Blend);

      for (int t = threads; t <= threads_max;
           t += std::max < int >(((threads_max - threads) / 30), 1))
        {
          std::cout << "Running Test with " << t << " threads." << std::endl;
          MT_Streaming_Kernels::Kernel_Serial_Base_Helper <
            NumberTest_SCALAR_None < data_size > >helper (op, data_size, t);

          double min_time = 10000000;
          double max_time = -1;
          double avg_time = 0;

          for (int i = 0; i < passes; i++)
            {
              start_timer ();
              helper.Run_Parallel ();
              stop_timer ();
              std::cout << get_time () << "s" << std::endl;
              min_time = std::min < double >(min_time, get_time ());
              max_time = std::max < double >(max_time, get_time ());
              avg_time += get_time ();
            }
          avg_time = avg_time / passes;
          std::cout << "Min pass time: " << min_time << std::endl;
          std::cout << "Max pass time: " << max_time << std::endl;
          std::cout << "Avg pass time: " << avg_time << std::endl;
        }




    }

//=======================================================
//
//             COMPUTE SSE RESULTS
//
//=======================================================
#ifdef ENABLE_SSE_INSTRUCTION_SET
    {
      std::cout << "	Running " << data_size << " of SSE :  " << std::endl;


      NumberTest_SSE_None < data_size > op ((float *) &A, (float *) &B,
                                            (float *) &C, (float *) &D,
                                            (float *) &Add, (float *) &Multiply,
                                            (float *) &Subtract,
                                            (float *) &Divide,
                                            (float *) &LessThan,
                                            (float *) &GreaterThan,
                                            (float *) &LessEquals,
                                            (float *) &GreaterEquals,
                                            (float *) &Equals, (float *) &And,
                                            (float *) &Or, (float *) &Xor,
                                            (float *) &AndNot, (float *) &Not,
                                            (float *) &Sqrt, (float *) &RSqrt,
                                            (float *) &Log, (float *) &Exp,
                                            (float *) &Inverse,
                                            (float *) &AbsoluteValue,
                                            (float *) &Sign, (float *) &Minimum,
                                            (float *) &Maximum,
                                            (float *) &Blend);

      for (int t = threads; t <= threads_max;
           t += std::max < int >(((threads_max - threads) / 30), 1))
        {
          std::cout << "Running Test with " << t << " threads." << std::endl;
          MT_Streaming_Kernels::Kernel_Serial_Base_Helper <
            NumberTest_SSE_None < data_size > >helper (op, data_size, t);

          double min_time = 10000000;
          double max_time = -1;
          double avg_time = 0;

          for (int i = 0; i < passes; i++)
            {
              start_timer ();
              helper.Run_Parallel ();
              stop_timer ();
              std::cout << get_time () << "s" << std::endl;
              min_time = std::min < double >(min_time, get_time ());
              max_time = std::max < double >(max_time, get_time ());
              avg_time += get_time ();
            }
          avg_time = avg_time / passes;
          std::cout << "Min pass time: " << min_time << std::endl;
          std::cout << "Max pass time: " << max_time << std::endl;
          std::cout << "Avg pass time: " << avg_time << std::endl;
        }




    }
#endif

//=======================================================
//
//             COMPUTE AVX RESULTS
//
//=======================================================
#ifdef ENABLE_AVX_INSTRUCTION_SET
    {
      std::cout << "	Running " << data_size << " of AVX :  " << std::endl;


      NumberTest_AVX_None < data_size > op ((float *) &A, (float *) &B,
                                            (float *) &C, (float *) &D,
                                            (float *) &Add, (float *) &Multiply,
                                            (float *) &Subtract,
                                            (float *) &Divide,
                                            (float *) &LessThan,
                                            (float *) &GreaterThan,
                                            (float *) &LessEquals,
                                            (float *) &GreaterEquals,
                                            (float *) &Equals, (float *) &And,
                                            (float *) &Or, (float *) &Xor,
                                            (float *) &AndNot, (float *) &Not,
                                            (float *) &Sqrt, (float *) &RSqrt,
                                            (float *) &Log, (float *) &Exp,
                                            (float *) &Inverse,
                                            (float *) &AbsoluteValue,
                                            (float *) &Sign, (float *) &Minimum,
                                            (float *) &Maximum,
                                            (float *) &Blend);

      for (int t = threads; t <= threads_max;
           t += std::max < int >(((threads_max - threads) / 30), 1))
        {
          std::cout << "Running Test with " << t << " threads." << std::endl;
          MT_Streaming_Kernels::Kernel_Serial_Base_Helper <
            NumberTest_AVX_None < data_size > >helper (op, data_size, t);

          double min_time = 10000000;
          double max_time = -1;
          double avg_time = 0;

          for (int i = 0; i < passes; i++)
            {
              start_timer ();
              helper.Run_Parallel ();
              stop_timer ();
              std::cout << get_time () << "s" << std::endl;
              min_time = std::min < double >(min_time, get_time ());
              max_time = std::max < double >(max_time, get_time ());
              avg_time += get_time ();
            }
          avg_time = avg_time / passes;
          std::cout << "Min pass time: " << min_time << std::endl;
          std::cout << "Max pass time: " << max_time << std::endl;
          std::cout << "Avg pass time: " << avg_time << std::endl;
        }




    }
#endif

//=======================================================
//
//             COMPUTE NEON RESULTS
//
//=======================================================
#ifdef ENABLE_NEON_INSTRUCTION_SET
    {
      std::cout << "	Running " << data_size << " of NEON :  " << std::endl;


      NumberTest_NEON_None < data_size > op ((float *) &A, (float *) &B,
                                             (float *) &C, (float *) &D,
                                             (float *) &Add,
                                             (float *) &Multiply,
                                             (float *) &Subtract,
                                             (float *) &Divide,
                                             (float *) &LessThan,
                                             (float *) &GreaterThan,
                                             (float *) &LessEquals,
                                             (float *) &GreaterEquals,
                                             (float *) &Equals, (float *) &And,
                                             (float *) &Or, (float *) &Xor,
                                             (float *) &AndNot, (float *) &Not,
                                             (float *) &Sqrt, (float *) &RSqrt,
                                             (float *) &Log, (float *) &Exp,
                                             (float *) &Inverse,
                                             (float *) &AbsoluteValue,
                                             (float *) &Sign,
                                             (float *) &Minimum,
                                             (float *) &Maximum,
                                             (float *) &Blend);

      for (int t = threads; t <= threads_max;
           t += std::max < int >(((threads_max - threads) / 30), 1))
        {
          std::cout << "Running Test with " << t << " threads." << std::endl;
          MT_Streaming_Kernels::Kernel_Serial_Base_Helper <
            NumberTest_NEON_None < data_size > >helper (op, data_size, t);

          double min_time = 10000000;
          double max_time = -1;
          double avg_time = 0;

          for (int i = 0; i < passes; i++)
            {
              start_timer ();
              helper.Run_Parallel ();
              stop_timer ();
              std::cout << get_time () << "s" << std::endl;
              min_time = std::min < double >(min_time, get_time ());
              max_time = std::max < double >(max_time, get_time ());
              avg_time += get_time ();
            }
          avg_time = avg_time / passes;
          std::cout << "Min pass time: " << min_time << std::endl;
          std::cout << "Max pass time: " << max_time << std::endl;
          std::cout << "Avg pass time: " << avg_time << std::endl;
        }




    }
#endif

//=======================================================
//
//             COMPUTE MIC RESULTS
//
//=======================================================
#ifdef ENABLE_MIC_INSTRUCTION_SET
    {
      std::cout << "	Running " << data_size << " of MIC :  " << std::endl;


      NumberTest_MIC_None < data_size > op ((float *) &A, (float *) &B,
                                            (float *) &C, (float *) &D,
                                            (float *) &Add, (float *) &Multiply,
                                            (float *) &Subtract,
                                            (float *) &Divide,
                                            (float *) &LessThan,
                                            (float *) &GreaterThan,
                                            (float *) &LessEquals,
                                            (float *) &GreaterEquals,
                                            (float *) &Equals, (float *) &And,
                                            (float *) &Or, (float *) &Xor,
                                            (float *) &AndNot, (float *) &Not,
                                            (float *) &Sqrt, (float *) &RSqrt,
                                            (float *) &Log, (float *) &Exp,
                                            (float *) &Inverse,
                                            (float *) &AbsoluteValue,
                                            (float *) &Sign, (float *) &Minimum,
                                            (float *) &Maximum,
                                            (float *) &Blend);

      for (int t = threads; t <= threads_max;
           t += std::max < int >(((threads_max - threads) / 30), 1))
        {
          std::cout << "Running Test with " << t << " threads." << std::endl;
          MT_Streaming_Kernels::Kernel_Serial_Base_Helper <
            NumberTest_MIC_None < data_size > >helper (op, data_size, t);

          double min_time = 10000000;
          double max_time = -1;
          double avg_time = 0;

          for (int i = 0; i < passes; i++)
            {
              start_timer ();
              helper.Run_Parallel ();
              stop_timer ();
              std::cout << get_time () << "s" << std::endl;
              min_time = std::min < double >(min_time, get_time ());
              max_time = std::max < double >(max_time, get_time ());
              avg_time += get_time ();
            }
          avg_time = avg_time / passes;
          std::cout << "Min pass time: " << min_time << std::endl;
          std::cout << "Max pass time: " << max_time << std::endl;
          std::cout << "Avg pass time: " << avg_time << std::endl;
        }




    }
#endif

//=======================================================
//
//        FREE MEMORY USED BY ALL VARIABLES
//
//=======================================================
    std::cout << "\nFreeing all data: " << std::endl;
    std::cout.flush ();

    _mm_free (reinterpret_cast < void *>(A));
    _mm_free (reinterpret_cast < void *>(B));
    _mm_free (reinterpret_cast < void *>(C));
    _mm_free (reinterpret_cast < void *>(D));
    _mm_free (reinterpret_cast < void *>(Add));
    _mm_free (reinterpret_cast < void *>(Multiply));
    _mm_free (reinterpret_cast < void *>(Subtract));
    _mm_free (reinterpret_cast < void *>(Divide));
    _mm_free (reinterpret_cast < void *>(LessThan));
    _mm_free (reinterpret_cast < void *>(GreaterThan));
    _mm_free (reinterpret_cast < void *>(LessEquals));
    _mm_free (reinterpret_cast < void *>(GreaterEquals));
    _mm_free (reinterpret_cast < void *>(Equals));
    _mm_free (reinterpret_cast < void *>(And));
    _mm_free (reinterpret_cast < void *>(Or));
    _mm_free (reinterpret_cast < void *>(Xor));
    _mm_free (reinterpret_cast < void *>(AndNot));
    _mm_free (reinterpret_cast < void *>(Not));
    _mm_free (reinterpret_cast < void *>(Sqrt));
    _mm_free (reinterpret_cast < void *>(RSqrt));
    _mm_free (reinterpret_cast < void *>(Log));
    _mm_free (reinterpret_cast < void *>(Exp));
    _mm_free (reinterpret_cast < void *>(Inverse));
    _mm_free (reinterpret_cast < void *>(AbsoluteValue));
    _mm_free (reinterpret_cast < void *>(Sign));
    _mm_free (reinterpret_cast < void *>(Minimum));
    _mm_free (reinterpret_cast < void *>(Maximum));
    _mm_free (reinterpret_cast < void *>(Blend));


  }


  return 0;

}

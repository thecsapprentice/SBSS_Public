
#include <cstdlib>
#include <iostream>
#include "KernelCommon.h"

struct NEOHOOKEAN_TAG;
struct COROTATED_TAG;

#include "NumberTest.h"
#include "NumberTest_Reference.h"

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
    std::cout << "Running SIMD Test for NumberTest " << std::endl;


//=======================================================
//
//        DEFINE ALL VARIABLES USED BY KERNEL
//
//=======================================================

    T               A[16] __attribute__ ((aligned (64)));
    T               B[16] __attribute__ ((aligned (64)));
    T               C[16] __attribute__ ((aligned (64)));
    T               D[16] __attribute__ ((aligned (64)));
    T               Add[16] __attribute__ ((aligned (64)));
    T               Add_reference[16] __attribute__ ((aligned (64)));
    T               Add_original[16] __attribute__ ((aligned (64)));
    T               Multiply[16] __attribute__ ((aligned (64)));
    T               Multiply_reference[16] __attribute__ ((aligned (64)));
    T               Multiply_original[16] __attribute__ ((aligned (64)));
    T               Subtract[16] __attribute__ ((aligned (64)));
    T               Subtract_reference[16] __attribute__ ((aligned (64)));
    T               Subtract_original[16] __attribute__ ((aligned (64)));
    T               Divide[16] __attribute__ ((aligned (64)));
    T               Divide_reference[16] __attribute__ ((aligned (64)));
    T               Divide_original[16] __attribute__ ((aligned (64)));
    T               LessThan[16] __attribute__ ((aligned (64)));
    T               LessThan_reference[16] __attribute__ ((aligned (64)));
    T               LessThan_original[16] __attribute__ ((aligned (64)));
    T               GreaterThan[16] __attribute__ ((aligned (64)));
    T               GreaterThan_reference[16] __attribute__ ((aligned (64)));
    T               GreaterThan_original[16] __attribute__ ((aligned (64)));
    T               LessEquals[16] __attribute__ ((aligned (64)));
    T               LessEquals_reference[16] __attribute__ ((aligned (64)));
    T               LessEquals_original[16] __attribute__ ((aligned (64)));
    T               GreaterEquals[16] __attribute__ ((aligned (64)));
    T               GreaterEquals_reference[16] __attribute__ ((aligned (64)));
    T               GreaterEquals_original[16] __attribute__ ((aligned (64)));
    T               Equals[16] __attribute__ ((aligned (64)));
    T               Equals_reference[16] __attribute__ ((aligned (64)));
    T               Equals_original[16] __attribute__ ((aligned (64)));
    T               And[16] __attribute__ ((aligned (64)));
    T               And_reference[16] __attribute__ ((aligned (64)));
    T               And_original[16] __attribute__ ((aligned (64)));
    T               Or[16] __attribute__ ((aligned (64)));
    T               Or_reference[16] __attribute__ ((aligned (64)));
    T               Or_original[16] __attribute__ ((aligned (64)));
    T               Xor[16] __attribute__ ((aligned (64)));
    T               Xor_reference[16] __attribute__ ((aligned (64)));
    T               Xor_original[16] __attribute__ ((aligned (64)));
    T               AndNot[16] __attribute__ ((aligned (64)));
    T               AndNot_reference[16] __attribute__ ((aligned (64)));
    T               AndNot_original[16] __attribute__ ((aligned (64)));
    T               Not[16] __attribute__ ((aligned (64)));
    T               Not_reference[16] __attribute__ ((aligned (64)));
    T               Not_original[16] __attribute__ ((aligned (64)));
    T               Sqrt[16] __attribute__ ((aligned (64)));
    T               Sqrt_reference[16] __attribute__ ((aligned (64)));
    T               Sqrt_original[16] __attribute__ ((aligned (64)));
    T               RSqrt[16] __attribute__ ((aligned (64)));
    T               RSqrt_reference[16] __attribute__ ((aligned (64)));
    T               RSqrt_original[16] __attribute__ ((aligned (64)));
    T               Log[16] __attribute__ ((aligned (64)));
    T               Log_reference[16] __attribute__ ((aligned (64)));
    T               Log_original[16] __attribute__ ((aligned (64)));
    T               Exp[16] __attribute__ ((aligned (64)));
    T               Exp_reference[16] __attribute__ ((aligned (64)));
    T               Exp_original[16] __attribute__ ((aligned (64)));
    T               Inverse[16] __attribute__ ((aligned (64)));
    T               Inverse_reference[16] __attribute__ ((aligned (64)));
    T               Inverse_original[16] __attribute__ ((aligned (64)));
    T               AbsoluteValue[16] __attribute__ ((aligned (64)));
    T               AbsoluteValue_reference[16] __attribute__ ((aligned (64)));
    T               AbsoluteValue_original[16] __attribute__ ((aligned (64)));
    T               Sign[16] __attribute__ ((aligned (64)));
    T               Sign_reference[16] __attribute__ ((aligned (64)));
    T               Sign_original[16] __attribute__ ((aligned (64)));
    T               Minimum[16] __attribute__ ((aligned (64)));
    T               Minimum_reference[16] __attribute__ ((aligned (64)));
    T               Minimum_original[16] __attribute__ ((aligned (64)));
    T               Maximum[16] __attribute__ ((aligned (64)));
    T               Maximum_reference[16] __attribute__ ((aligned (64)));
    T               Maximum_original[16] __attribute__ ((aligned (64)));
    T               Blend[16] __attribute__ ((aligned (64)));
    T               Blend_reference[16] __attribute__ ((aligned (64)));
    T               Blend_original[16] __attribute__ ((aligned (64)));


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
//             COMPUTE REFERENCE RESULTS
//
//=======================================================

    T __mA __attribute__ ((aligned (4)));
    T __mB __attribute__ ((aligned (4)));
    T __mC __attribute__ ((aligned (4)));
    T __mD __attribute__ ((aligned (4)));
    T __mAdd __attribute__ ((aligned (4)));
    T __mAdd_reference __attribute__ ((aligned (4)));
    T __mAdd_original __attribute__ ((aligned (4)));
    T __mMultiply __attribute__ ((aligned (4)));
    T __mMultiply_reference __attribute__ ((aligned (4)));
    T __mMultiply_original __attribute__ ((aligned (4)));
    T __mSubtract __attribute__ ((aligned (4)));
    T __mSubtract_reference __attribute__ ((aligned (4)));
    T __mSubtract_original __attribute__ ((aligned (4)));
    T __mDivide __attribute__ ((aligned (4)));
    T __mDivide_reference __attribute__ ((aligned (4)));
    T __mDivide_original __attribute__ ((aligned (4)));
    T __mLessThan __attribute__ ((aligned (4)));
    T __mLessThan_reference __attribute__ ((aligned (4)));
    T __mLessThan_original __attribute__ ((aligned (4)));
    T __mGreaterThan __attribute__ ((aligned (4)));
    T __mGreaterThan_reference __attribute__ ((aligned (4)));
    T __mGreaterThan_original __attribute__ ((aligned (4)));
    T __mLessEquals __attribute__ ((aligned (4)));
    T __mLessEquals_reference __attribute__ ((aligned (4)));
    T __mLessEquals_original __attribute__ ((aligned (4)));
    T __mGreaterEquals __attribute__ ((aligned (4)));
    T __mGreaterEquals_reference __attribute__ ((aligned (4)));
    T __mGreaterEquals_original __attribute__ ((aligned (4)));
    T __mEquals __attribute__ ((aligned (4)));
    T __mEquals_reference __attribute__ ((aligned (4)));
    T __mEquals_original __attribute__ ((aligned (4)));
    T __mAnd __attribute__ ((aligned (4)));
    T __mAnd_reference __attribute__ ((aligned (4)));
    T __mAnd_original __attribute__ ((aligned (4)));
    T __mOr __attribute__ ((aligned (4)));
    T __mOr_reference __attribute__ ((aligned (4)));
    T __mOr_original __attribute__ ((aligned (4)));
    T __mXor __attribute__ ((aligned (4)));
    T __mXor_reference __attribute__ ((aligned (4)));
    T __mXor_original __attribute__ ((aligned (4)));
    T __mAndNot __attribute__ ((aligned (4)));
    T __mAndNot_reference __attribute__ ((aligned (4)));
    T __mAndNot_original __attribute__ ((aligned (4)));
    T __mNot __attribute__ ((aligned (4)));
    T __mNot_reference __attribute__ ((aligned (4)));
    T __mNot_original __attribute__ ((aligned (4)));
    T __mSqrt __attribute__ ((aligned (4)));
    T __mSqrt_reference __attribute__ ((aligned (4)));
    T __mSqrt_original __attribute__ ((aligned (4)));
    T __mRSqrt __attribute__ ((aligned (4)));
    T __mRSqrt_reference __attribute__ ((aligned (4)));
    T __mRSqrt_original __attribute__ ((aligned (4)));
    T __mLog __attribute__ ((aligned (4)));
    T __mLog_reference __attribute__ ((aligned (4)));
    T __mLog_original __attribute__ ((aligned (4)));
    T __mExp __attribute__ ((aligned (4)));
    T __mExp_reference __attribute__ ((aligned (4)));
    T __mExp_original __attribute__ ((aligned (4)));
    T __mInverse __attribute__ ((aligned (4)));
    T __mInverse_reference __attribute__ ((aligned (4)));
    T __mInverse_original __attribute__ ((aligned (4)));
    T __mAbsoluteValue __attribute__ ((aligned (4)));
    T __mAbsoluteValue_reference __attribute__ ((aligned (4)));
    T __mAbsoluteValue_original __attribute__ ((aligned (4)));
    T __mSign __attribute__ ((aligned (4)));
    T __mSign_reference __attribute__ ((aligned (4)));
    T __mSign_original __attribute__ ((aligned (4)));
    T __mMinimum __attribute__ ((aligned (4)));
    T __mMinimum_reference __attribute__ ((aligned (4)));
    T __mMinimum_original __attribute__ ((aligned (4)));
    T __mMaximum __attribute__ ((aligned (4)));
    T __mMaximum_reference __attribute__ ((aligned (4)));
    T __mMaximum_original __attribute__ ((aligned (4)));
    T __mBlend __attribute__ ((aligned (4)));
    T __mBlend_reference __attribute__ ((aligned (4)));
    T __mBlend_original __attribute__ ((aligned (4)));
    for (int k = 0; k < 16; k++)
    {
      __mA = A[k];
      __mB = B[k];
      __mC = C[k];
      __mD = D[k];
      __mAdd_reference = Add_reference[k];

      __mMultiply_reference = Multiply_reference[k];

      __mSubtract_reference = Subtract_reference[k];

      __mDivide_reference = Divide_reference[k];

      __mLessThan_reference = LessThan_reference[k];

      __mGreaterThan_reference = GreaterThan_reference[k];

      __mLessEquals_reference = LessEquals_reference[k];

      __mGreaterEquals_reference = GreaterEquals_reference[k];

      __mEquals_reference = Equals_reference[k];

      __mAnd_reference = And_reference[k];

      __mOr_reference = Or_reference[k];

      __mXor_reference = Xor_reference[k];

      __mAndNot_reference = AndNot_reference[k];

      __mNot_reference = Not_reference[k];

      __mSqrt_reference = Sqrt_reference[k];

      __mRSqrt_reference = RSqrt_reference[k];

      __mLog_reference = Log_reference[k];

      __mExp_reference = Exp_reference[k];

      __mInverse_reference = Inverse_reference[k];

      __mAbsoluteValue_reference = AbsoluteValue_reference[k];

      __mSign_reference = Sign_reference[k];

      __mMinimum_reference = Minimum_reference[k];

      __mMaximum_reference = Maximum_reference[k];

      __mBlend_reference = Blend_reference[k];
      NumberTest < float, float, int >(__mA, __mB, __mC, __mD, __mAdd_reference,
                                       __mMultiply_reference,
                                       __mSubtract_reference,
                                       __mDivide_reference,
                                       __mLessThan_reference,
                                       __mGreaterThan_reference,
                                       __mLessEquals_reference,
                                       __mGreaterEquals_reference,
                                       __mEquals_reference, __mAnd_reference,
                                       __mOr_reference, __mXor_reference,
                                       __mAndNot_reference, __mNot_reference,
                                       __mSqrt_reference, __mRSqrt_reference,
                                       __mLog_reference, __mExp_reference,
                                       __mInverse_reference,
                                       __mAbsoluteValue_reference,
                                       __mSign_reference, __mMinimum_reference,
                                       __mMaximum_reference,
                                       __mBlend_reference);
      Add_reference[k] = __mAdd_reference;
      Multiply_reference[k] = __mMultiply_reference;
      Subtract_reference[k] = __mSubtract_reference;
      Divide_reference[k] = __mDivide_reference;
      LessThan_reference[k] = __mLessThan_reference;
      GreaterThan_reference[k] = __mGreaterThan_reference;
      LessEquals_reference[k] = __mLessEquals_reference;
      GreaterEquals_reference[k] = __mGreaterEquals_reference;
      Equals_reference[k] = __mEquals_reference;
      And_reference[k] = __mAnd_reference;
      Or_reference[k] = __mOr_reference;
      Xor_reference[k] = __mXor_reference;
      AndNot_reference[k] = __mAndNot_reference;
      Not_reference[k] = __mNot_reference;
      Sqrt_reference[k] = __mSqrt_reference;
      RSqrt_reference[k] = __mRSqrt_reference;
      Log_reference[k] = __mLog_reference;
      Exp_reference[k] = __mExp_reference;
      Inverse_reference[k] = __mInverse_reference;
      AbsoluteValue_reference[k] = __mAbsoluteValue_reference;
      Sign_reference[k] = __mSign_reference;
      Minimum_reference[k] = __mMinimum_reference;
      Maximum_reference[k] = __mMaximum_reference;
      Blend_reference[k] = __mBlend_reference;
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
      typedef         T (&refArray6)[16];
      typedef         T (&refArray7)[16];
      typedef         T (&refArray8)[16];
      typedef         T (&refArray9)[16];
      typedef         T (&refArray10)[16];
      typedef         T (&refArray11)[16];
      typedef         T (&refArray12)[16];
      typedef         T (&refArray13)[16];
      typedef         T (&refArray14)[16];
      typedef         T (&refArray15)[16];
      typedef         T (&refArray16)[16];
      typedef         T (&refArray17)[16];
      typedef         T (&refArray18)[16];
      typedef         T (&refArray19)[16];
      typedef         T (&refArray20)[16];
      typedef         T (&refArray21)[16];
      typedef         T (&refArray22)[16];
      typedef         T (&refArray23)[16];
      typedef         T (&refArray24)[16];
      typedef         T (&refArray25)[16];
      typedef         T (&refArray26)[16];
      typedef         T (&refArray27)[16];
      typedef         T (&refArray28)[16];
      for (int __a = 0; __a < 16; __a++)
        Add[__a] = Add_original[__a];
      for (int __a = 0; __a < 16; __a++)
        Multiply[__a] = Multiply_original[__a];
      for (int __a = 0; __a < 16; __a++)
        Subtract[__a] = Subtract_original[__a];
      for (int __a = 0; __a < 16; __a++)
        Divide[__a] = Divide_original[__a];
      for (int __a = 0; __a < 16; __a++)
        LessThan[__a] = LessThan_original[__a];
      for (int __a = 0; __a < 16; __a++)
        GreaterThan[__a] = GreaterThan_original[__a];
      for (int __a = 0; __a < 16; __a++)
        LessEquals[__a] = LessEquals_original[__a];
      for (int __a = 0; __a < 16; __a++)
        GreaterEquals[__a] = GreaterEquals_original[__a];
      for (int __a = 0; __a < 16; __a++)
        Equals[__a] = Equals_original[__a];
      for (int __a = 0; __a < 16; __a++)
        And[__a] = And_original[__a];
      for (int __a = 0; __a < 16; __a++)
        Or[__a] = Or_original[__a];
      for (int __a = 0; __a < 16; __a++)
        Xor[__a] = Xor_original[__a];
      for (int __a = 0; __a < 16; __a++)
        AndNot[__a] = AndNot_original[__a];
      for (int __a = 0; __a < 16; __a++)
        Not[__a] = Not_original[__a];
      for (int __a = 0; __a < 16; __a++)
        Sqrt[__a] = Sqrt_original[__a];
      for (int __a = 0; __a < 16; __a++)
        RSqrt[__a] = RSqrt_original[__a];
      for (int __a = 0; __a < 16; __a++)
        Log[__a] = Log_original[__a];
      for (int __a = 0; __a < 16; __a++)
        Exp[__a] = Exp_original[__a];
      for (int __a = 0; __a < 16; __a++)
        Inverse[__a] = Inverse_original[__a];
      for (int __a = 0; __a < 16; __a++)
        AbsoluteValue[__a] = AbsoluteValue_original[__a];
      for (int __a = 0; __a < 16; __a++)
        Sign[__a] = Sign_original[__a];
      for (int __a = 0; __a < 16; __a++)
        Minimum[__a] = Minimum_original[__a];
      for (int __a = 0; __a < 16; __a++)
        Maximum[__a] = Maximum_original[__a];
      for (int __a = 0; __a < 16; __a++)
        Blend[__a] = Blend_original[__a];
      for (int i = 0; i < 16; i += 1)
      {
        refArray1       Ak = reinterpret_cast < refArray1 > (A[i]);
        refArray2       Bk = reinterpret_cast < refArray2 > (B[i]);
        refArray3       Ck = reinterpret_cast < refArray3 > (C[i]);
        refArray4       Dk = reinterpret_cast < refArray4 > (D[i]);
        refArray5       Addk = reinterpret_cast < refArray5 > (Add[i]);
        refArray6       Multiplyk =
          reinterpret_cast < refArray6 > (Multiply[i]);
        refArray7       Subtractk =
          reinterpret_cast < refArray7 > (Subtract[i]);
        refArray8       Dividek = reinterpret_cast < refArray8 > (Divide[i]);
        refArray9       LessThank =
          reinterpret_cast < refArray9 > (LessThan[i]);
        refArray10      GreaterThank =
          reinterpret_cast < refArray10 > (GreaterThan[i]);
        refArray11      LessEqualsk =
          reinterpret_cast < refArray11 > (LessEquals[i]);
        refArray12      GreaterEqualsk =
          reinterpret_cast < refArray12 > (GreaterEquals[i]);
        refArray13      Equalsk = reinterpret_cast < refArray13 > (Equals[i]);
        refArray14      Andk = reinterpret_cast < refArray14 > (And[i]);
        refArray15      Ork = reinterpret_cast < refArray15 > (Or[i]);
        refArray16      Xork = reinterpret_cast < refArray16 > (Xor[i]);
        refArray17      AndNotk = reinterpret_cast < refArray17 > (AndNot[i]);
        refArray18      Notk = reinterpret_cast < refArray18 > (Not[i]);
        refArray19      Sqrtk = reinterpret_cast < refArray19 > (Sqrt[i]);
        refArray20      RSqrtk = reinterpret_cast < refArray20 > (RSqrt[i]);
        refArray21      Logk = reinterpret_cast < refArray21 > (Log[i]);
        refArray22      Expk = reinterpret_cast < refArray22 > (Exp[i]);
        refArray23      Inversek = reinterpret_cast < refArray23 > (Inverse[i]);
        refArray24      AbsoluteValuek =
          reinterpret_cast < refArray24 > (AbsoluteValue[i]);
        refArray25      Signk = reinterpret_cast < refArray25 > (Sign[i]);
        refArray26      Minimumk = reinterpret_cast < refArray26 > (Minimum[i]);
        refArray27      Maximumk = reinterpret_cast < refArray27 > (Maximum[i]);
        refArray28      Blendk = reinterpret_cast < refArray28 > (Blend[i]);
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
      for (int __a = 0; __a < 16; __a++)
        if (std::abs ((Add[__a] - Add_reference[__a]) / (Add_reference[__a])) >
            1)
        {
          std::cerr << "Mismatch detected in SCALAR implementation" << std::
            endl;
          std::cerr << "Variable Add:" << std::endl;
          std::cerr << "seed=" << seed << ", __a=" << __a << std::endl;
          std::cerr << "Add SCALAR=  " << Add[__a] << std::endl;
          std::cerr << "Add Reference=  " << Add_reference[__a] << std::endl;
          std::cerr << "Add Rel Difference=  " << std::
            abs ((Add[__a] -
                  Add_reference[__a]) / (Add_reference[__a])) << std::endl;
          std::cerr << "Add Abs Difference=  " << std::abs (Add[__a] -
                                                            Add_reference[__a])
            << std::endl;
          return 1;
        }
      for (int __a = 0; __a < 16; __a++)
        if (std::
            abs ((Multiply[__a] -
                  Multiply_reference[__a]) / (Multiply_reference[__a])) > 1)
        {
          std::cerr << "Mismatch detected in SCALAR implementation" << std::
            endl;
          std::cerr << "Variable Multiply:" << std::endl;
          std::cerr << "seed=" << seed << ", __a=" << __a << std::endl;
          std::cerr << "Multiply SCALAR=  " << Multiply[__a] << std::endl;
          std::
            cerr << "Multiply Reference=  " << Multiply_reference[__a] << std::
            endl;
          std::cerr << "Multiply Rel Difference=  " << std::
            abs ((Multiply[__a] -
                  Multiply_reference[__a]) /
                 (Multiply_reference[__a])) << std::endl;
          std::cerr << "Multiply Abs Difference=  " << std::abs (Multiply[__a] -
                                                                 Multiply_reference
                                                                 [__a]) << std::
            endl;
          return 1;
        }
      for (int __a = 0; __a < 16; __a++)
        if (std::
            abs ((Subtract[__a] -
                  Subtract_reference[__a]) / (Subtract_reference[__a])) > 1)
        {
          std::cerr << "Mismatch detected in SCALAR implementation" << std::
            endl;
          std::cerr << "Variable Subtract:" << std::endl;
          std::cerr << "seed=" << seed << ", __a=" << __a << std::endl;
          std::cerr << "Subtract SCALAR=  " << Subtract[__a] << std::endl;
          std::
            cerr << "Subtract Reference=  " << Subtract_reference[__a] << std::
            endl;
          std::cerr << "Subtract Rel Difference=  " << std::
            abs ((Subtract[__a] -
                  Subtract_reference[__a]) /
                 (Subtract_reference[__a])) << std::endl;
          std::cerr << "Subtract Abs Difference=  " << std::abs (Subtract[__a] -
                                                                 Subtract_reference
                                                                 [__a]) << std::
            endl;
          return 1;
        }
      for (int __a = 0; __a < 16; __a++)
        if (std::
            abs ((Divide[__a] -
                  Divide_reference[__a]) / (Divide_reference[__a])) > 1)
        {
          std::cerr << "Mismatch detected in SCALAR implementation" << std::
            endl;
          std::cerr << "Variable Divide:" << std::endl;
          std::cerr << "seed=" << seed << ", __a=" << __a << std::endl;
          std::cerr << "Divide SCALAR=  " << Divide[__a] << std::endl;
          std::cerr << "Divide Reference=  " << Divide_reference[__a] << std::
            endl;
          std::cerr << "Divide Rel Difference=  " << std::
            abs ((Divide[__a] -
                  Divide_reference[__a]) /
                 (Divide_reference[__a])) << std::endl;
          std::cerr << "Divide Abs Difference=  " << std::abs (Divide[__a] -
                                                               Divide_reference
                                                               [__a]) << std::
            endl;
          return 1;
        }
      for (int __a = 0; __a < 16; __a++)
        if (std::
            abs ((LessThan[__a] -
                  LessThan_reference[__a]) / (LessThan_reference[__a])) > 1)
        {
          std::cerr << "Mismatch detected in SCALAR implementation" << std::
            endl;
          std::cerr << "Variable LessThan:" << std::endl;
          std::cerr << "seed=" << seed << ", __a=" << __a << std::endl;
          std::cerr << "LessThan SCALAR=  " << LessThan[__a] << std::endl;
          std::
            cerr << "LessThan Reference=  " << LessThan_reference[__a] << std::
            endl;
          std::cerr << "LessThan Rel Difference=  " << std::
            abs ((LessThan[__a] -
                  LessThan_reference[__a]) /
                 (LessThan_reference[__a])) << std::endl;
          std::cerr << "LessThan Abs Difference=  " << std::abs (LessThan[__a] -
                                                                 LessThan_reference
                                                                 [__a]) << std::
            endl;
          return 1;
        }
      for (int __a = 0; __a < 16; __a++)
        if (std::
            abs ((GreaterThan[__a] -
                  GreaterThan_reference[__a]) / (GreaterThan_reference[__a])) >
            1)
        {
          std::cerr << "Mismatch detected in SCALAR implementation" << std::
            endl;
          std::cerr << "Variable GreaterThan:" << std::endl;
          std::cerr << "seed=" << seed << ", __a=" << __a << std::endl;
          std::cerr << "GreaterThan SCALAR=  " << GreaterThan[__a] << std::endl;
          std::
            cerr << "GreaterThan Reference=  " << GreaterThan_reference[__a] <<
            std::endl;
          std::cerr << "GreaterThan Rel Difference=  " << std::
            abs ((GreaterThan[__a] -
                  GreaterThan_reference[__a]) /
                 (GreaterThan_reference[__a])) << std::endl;
          std::cerr << "GreaterThan Abs Difference=  " << std::
            abs (GreaterThan[__a] - GreaterThan_reference[__a]) << std::endl;
          return 1;
        }
      for (int __a = 0; __a < 16; __a++)
        if (std::
            abs ((LessEquals[__a] -
                  LessEquals_reference[__a]) / (LessEquals_reference[__a])) > 1)
        {
          std::cerr << "Mismatch detected in SCALAR implementation" << std::
            endl;
          std::cerr << "Variable LessEquals:" << std::endl;
          std::cerr << "seed=" << seed << ", __a=" << __a << std::endl;
          std::cerr << "LessEquals SCALAR=  " << LessEquals[__a] << std::endl;
          std::
            cerr << "LessEquals Reference=  " << LessEquals_reference[__a] <<
            std::endl;
          std::cerr << "LessEquals Rel Difference=  " << std::
            abs ((LessEquals[__a] -
                  LessEquals_reference[__a]) /
                 (LessEquals_reference[__a])) << std::endl;
          std::cerr << "LessEquals Abs Difference=  " << std::
            abs (LessEquals[__a] - LessEquals_reference[__a]) << std::endl;
          return 1;
        }
      for (int __a = 0; __a < 16; __a++)
        if (std::
            abs ((GreaterEquals[__a] -
                  GreaterEquals_reference[__a]) /
                 (GreaterEquals_reference[__a])) > 1)
        {
          std::cerr << "Mismatch detected in SCALAR implementation" << std::
            endl;
          std::cerr << "Variable GreaterEquals:" << std::endl;
          std::cerr << "seed=" << seed << ", __a=" << __a << std::endl;
          std::cerr << "GreaterEquals SCALAR=  " << GreaterEquals[__a] << std::
            endl;
          std::
            cerr << "GreaterEquals Reference=  " << GreaterEquals_reference[__a]
            << std::endl;
          std::cerr << "GreaterEquals Rel Difference=  " << std::
            abs ((GreaterEquals[__a] -
                  GreaterEquals_reference[__a]) /
                 (GreaterEquals_reference[__a])) << std::endl;
          std::cerr << "GreaterEquals Abs Difference=  " << std::
            abs (GreaterEquals[__a] -
                 GreaterEquals_reference[__a]) << std::endl;
          return 1;
        }
      for (int __a = 0; __a < 16; __a++)
        if (std::
            abs ((Equals[__a] -
                  Equals_reference[__a]) / (Equals_reference[__a])) > 1)
        {
          std::cerr << "Mismatch detected in SCALAR implementation" << std::
            endl;
          std::cerr << "Variable Equals:" << std::endl;
          std::cerr << "seed=" << seed << ", __a=" << __a << std::endl;
          std::cerr << "Equals SCALAR=  " << Equals[__a] << std::endl;
          std::cerr << "Equals Reference=  " << Equals_reference[__a] << std::
            endl;
          std::cerr << "Equals Rel Difference=  " << std::
            abs ((Equals[__a] -
                  Equals_reference[__a]) /
                 (Equals_reference[__a])) << std::endl;
          std::cerr << "Equals Abs Difference=  " << std::abs (Equals[__a] -
                                                               Equals_reference
                                                               [__a]) << std::
            endl;
          return 1;
        }
      for (int __a = 0; __a < 16; __a++)
        if (std::abs ((And[__a] - And_reference[__a]) / (And_reference[__a])) >
            1)
        {
          std::cerr << "Mismatch detected in SCALAR implementation" << std::
            endl;
          std::cerr << "Variable And:" << std::endl;
          std::cerr << "seed=" << seed << ", __a=" << __a << std::endl;
          std::cerr << "And SCALAR=  " << And[__a] << std::endl;
          std::cerr << "And Reference=  " << And_reference[__a] << std::endl;
          std::cerr << "And Rel Difference=  " << std::
            abs ((And[__a] -
                  And_reference[__a]) / (And_reference[__a])) << std::endl;
          std::cerr << "And Abs Difference=  " << std::abs (And[__a] -
                                                            And_reference[__a])
            << std::endl;
          return 1;
        }
      for (int __a = 0; __a < 16; __a++)
        if (std::abs ((Or[__a] - Or_reference[__a]) / (Or_reference[__a])) > 1)
        {
          std::cerr << "Mismatch detected in SCALAR implementation" << std::
            endl;
          std::cerr << "Variable Or:" << std::endl;
          std::cerr << "seed=" << seed << ", __a=" << __a << std::endl;
          std::cerr << "Or SCALAR=  " << Or[__a] << std::endl;
          std::cerr << "Or Reference=  " << Or_reference[__a] << std::endl;
          std::cerr << "Or Rel Difference=  " << std::
            abs ((Or[__a] -
                  Or_reference[__a]) / (Or_reference[__a])) << std::endl;
          std::cerr << "Or Abs Difference=  " << std::abs (Or[__a] -
                                                           Or_reference[__a]) <<
            std::endl;
          return 1;
        }
      for (int __a = 0; __a < 16; __a++)
        if (std::abs ((Xor[__a] - Xor_reference[__a]) / (Xor_reference[__a])) >
            1)
        {
          std::cerr << "Mismatch detected in SCALAR implementation" << std::
            endl;
          std::cerr << "Variable Xor:" << std::endl;
          std::cerr << "seed=" << seed << ", __a=" << __a << std::endl;
          std::cerr << "Xor SCALAR=  " << Xor[__a] << std::endl;
          std::cerr << "Xor Reference=  " << Xor_reference[__a] << std::endl;
          std::cerr << "Xor Rel Difference=  " << std::
            abs ((Xor[__a] -
                  Xor_reference[__a]) / (Xor_reference[__a])) << std::endl;
          std::cerr << "Xor Abs Difference=  " << std::abs (Xor[__a] -
                                                            Xor_reference[__a])
            << std::endl;
          return 1;
        }
      for (int __a = 0; __a < 16; __a++)
        if (std::
            abs ((AndNot[__a] -
                  AndNot_reference[__a]) / (AndNot_reference[__a])) > 1)
        {
          std::cerr << "Mismatch detected in SCALAR implementation" << std::
            endl;
          std::cerr << "Variable AndNot:" << std::endl;
          std::cerr << "seed=" << seed << ", __a=" << __a << std::endl;
          std::cerr << "AndNot SCALAR=  " << AndNot[__a] << std::endl;
          std::cerr << "AndNot Reference=  " << AndNot_reference[__a] << std::
            endl;
          std::cerr << "AndNot Rel Difference=  " << std::
            abs ((AndNot[__a] -
                  AndNot_reference[__a]) /
                 (AndNot_reference[__a])) << std::endl;
          std::cerr << "AndNot Abs Difference=  " << std::abs (AndNot[__a] -
                                                               AndNot_reference
                                                               [__a]) << std::
            endl;
          return 1;
        }
      for (int __a = 0; __a < 16; __a++)
        if (std::abs ((Not[__a] - Not_reference[__a]) / (Not_reference[__a])) >
            1)
        {
          std::cerr << "Mismatch detected in SCALAR implementation" << std::
            endl;
          std::cerr << "Variable Not:" << std::endl;
          std::cerr << "seed=" << seed << ", __a=" << __a << std::endl;
          std::cerr << "Not SCALAR=  " << Not[__a] << std::endl;
          std::cerr << "Not Reference=  " << Not_reference[__a] << std::endl;
          std::cerr << "Not Rel Difference=  " << std::
            abs ((Not[__a] -
                  Not_reference[__a]) / (Not_reference[__a])) << std::endl;
          std::cerr << "Not Abs Difference=  " << std::abs (Not[__a] -
                                                            Not_reference[__a])
            << std::endl;
          return 1;
        }
      for (int __a = 0; __a < 16; __a++)
        if (std::
            abs ((Sqrt[__a] - Sqrt_reference[__a]) / (Sqrt_reference[__a])) > 1)
        {
          std::cerr << "Mismatch detected in SCALAR implementation" << std::
            endl;
          std::cerr << "Variable Sqrt:" << std::endl;
          std::cerr << "seed=" << seed << ", __a=" << __a << std::endl;
          std::cerr << "Sqrt SCALAR=  " << Sqrt[__a] << std::endl;
          std::cerr << "Sqrt Reference=  " << Sqrt_reference[__a] << std::endl;
          std::cerr << "Sqrt Rel Difference=  " << std::
            abs ((Sqrt[__a] -
                  Sqrt_reference[__a]) / (Sqrt_reference[__a])) << std::endl;
          std::cerr << "Sqrt Abs Difference=  " << std::abs (Sqrt[__a] -
                                                             Sqrt_reference
                                                             [__a]) << std::
            endl;
          return 1;
        }
      for (int __a = 0; __a < 16; __a++)
        if (std::
            abs ((RSqrt[__a] - RSqrt_reference[__a]) / (RSqrt_reference[__a])) >
            1)
        {
          std::cerr << "Mismatch detected in SCALAR implementation" << std::
            endl;
          std::cerr << "Variable RSqrt:" << std::endl;
          std::cerr << "seed=" << seed << ", __a=" << __a << std::endl;
          std::cerr << "RSqrt SCALAR=  " << RSqrt[__a] << std::endl;
          std::cerr << "RSqrt Reference=  " << RSqrt_reference[__a] << std::
            endl;
          std::cerr << "RSqrt Rel Difference=  " << std::
            abs ((RSqrt[__a] -
                  RSqrt_reference[__a]) / (RSqrt_reference[__a])) << std::endl;
          std::cerr << "RSqrt Abs Difference=  " << std::abs (RSqrt[__a] -
                                                              RSqrt_reference
                                                              [__a]) << std::
            endl;
          return 1;
        }
      for (int __a = 0; __a < 16; __a++)
        if (std::abs ((Log[__a] - Log_reference[__a]) / (Log_reference[__a])) >
            1)
        {
          std::cerr << "Mismatch detected in SCALAR implementation" << std::
            endl;
          std::cerr << "Variable Log:" << std::endl;
          std::cerr << "seed=" << seed << ", __a=" << __a << std::endl;
          std::cerr << "Log SCALAR=  " << Log[__a] << std::endl;
          std::cerr << "Log Reference=  " << Log_reference[__a] << std::endl;
          std::cerr << "Log Rel Difference=  " << std::
            abs ((Log[__a] -
                  Log_reference[__a]) / (Log_reference[__a])) << std::endl;
          std::cerr << "Log Abs Difference=  " << std::abs (Log[__a] -
                                                            Log_reference[__a])
            << std::endl;
          return 1;
        }
      for (int __a = 0; __a < 16; __a++)
        if (std::abs ((Exp[__a] - Exp_reference[__a]) / (Exp_reference[__a])) >
            1)
        {
          std::cerr << "Mismatch detected in SCALAR implementation" << std::
            endl;
          std::cerr << "Variable Exp:" << std::endl;
          std::cerr << "seed=" << seed << ", __a=" << __a << std::endl;
          std::cerr << "Exp SCALAR=  " << Exp[__a] << std::endl;
          std::cerr << "Exp Reference=  " << Exp_reference[__a] << std::endl;
          std::cerr << "Exp Rel Difference=  " << std::
            abs ((Exp[__a] -
                  Exp_reference[__a]) / (Exp_reference[__a])) << std::endl;
          std::cerr << "Exp Abs Difference=  " << std::abs (Exp[__a] -
                                                            Exp_reference[__a])
            << std::endl;
          return 1;
        }
      for (int __a = 0; __a < 16; __a++)
        if (std::
            abs ((Inverse[__a] -
                  Inverse_reference[__a]) / (Inverse_reference[__a])) > 1)
        {
          std::cerr << "Mismatch detected in SCALAR implementation" << std::
            endl;
          std::cerr << "Variable Inverse:" << std::endl;
          std::cerr << "seed=" << seed << ", __a=" << __a << std::endl;
          std::cerr << "Inverse SCALAR=  " << Inverse[__a] << std::endl;
          std::cerr << "Inverse Reference=  " << Inverse_reference[__a] << std::
            endl;
          std::cerr << "Inverse Rel Difference=  " << std::
            abs ((Inverse[__a] -
                  Inverse_reference[__a]) /
                 (Inverse_reference[__a])) << std::endl;
          std::cerr << "Inverse Abs Difference=  " << std::abs (Inverse[__a] -
                                                                Inverse_reference
                                                                [__a]) << std::
            endl;
          return 1;
        }
      for (int __a = 0; __a < 16; __a++)
        if (std::
            abs ((AbsoluteValue[__a] -
                  AbsoluteValue_reference[__a]) /
                 (AbsoluteValue_reference[__a])) > 1)
        {
          std::cerr << "Mismatch detected in SCALAR implementation" << std::
            endl;
          std::cerr << "Variable AbsoluteValue:" << std::endl;
          std::cerr << "seed=" << seed << ", __a=" << __a << std::endl;
          std::cerr << "AbsoluteValue SCALAR=  " << AbsoluteValue[__a] << std::
            endl;
          std::
            cerr << "AbsoluteValue Reference=  " << AbsoluteValue_reference[__a]
            << std::endl;
          std::cerr << "AbsoluteValue Rel Difference=  " << std::
            abs ((AbsoluteValue[__a] -
                  AbsoluteValue_reference[__a]) /
                 (AbsoluteValue_reference[__a])) << std::endl;
          std::cerr << "AbsoluteValue Abs Difference=  " << std::
            abs (AbsoluteValue[__a] -
                 AbsoluteValue_reference[__a]) << std::endl;
          return 1;
        }
      for (int __a = 0; __a < 16; __a++)
        if (std::
            abs ((Sign[__a] - Sign_reference[__a]) / (Sign_reference[__a])) > 1)
        {
          std::cerr << "Mismatch detected in SCALAR implementation" << std::
            endl;
          std::cerr << "Variable Sign:" << std::endl;
          std::cerr << "seed=" << seed << ", __a=" << __a << std::endl;
          std::cerr << "Sign SCALAR=  " << Sign[__a] << std::endl;
          std::cerr << "Sign Reference=  " << Sign_reference[__a] << std::endl;
          std::cerr << "Sign Rel Difference=  " << std::
            abs ((Sign[__a] -
                  Sign_reference[__a]) / (Sign_reference[__a])) << std::endl;
          std::cerr << "Sign Abs Difference=  " << std::abs (Sign[__a] -
                                                             Sign_reference
                                                             [__a]) << std::
            endl;
          return 1;
        }
      for (int __a = 0; __a < 16; __a++)
        if (std::
            abs ((Minimum[__a] -
                  Minimum_reference[__a]) / (Minimum_reference[__a])) > 1)
        {
          std::cerr << "Mismatch detected in SCALAR implementation" << std::
            endl;
          std::cerr << "Variable Minimum:" << std::endl;
          std::cerr << "seed=" << seed << ", __a=" << __a << std::endl;
          std::cerr << "Minimum SCALAR=  " << Minimum[__a] << std::endl;
          std::cerr << "Minimum Reference=  " << Minimum_reference[__a] << std::
            endl;
          std::cerr << "Minimum Rel Difference=  " << std::
            abs ((Minimum[__a] -
                  Minimum_reference[__a]) /
                 (Minimum_reference[__a])) << std::endl;
          std::cerr << "Minimum Abs Difference=  " << std::abs (Minimum[__a] -
                                                                Minimum_reference
                                                                [__a]) << std::
            endl;
          return 1;
        }
      for (int __a = 0; __a < 16; __a++)
        if (std::
            abs ((Maximum[__a] -
                  Maximum_reference[__a]) / (Maximum_reference[__a])) > 1)
        {
          std::cerr << "Mismatch detected in SCALAR implementation" << std::
            endl;
          std::cerr << "Variable Maximum:" << std::endl;
          std::cerr << "seed=" << seed << ", __a=" << __a << std::endl;
          std::cerr << "Maximum SCALAR=  " << Maximum[__a] << std::endl;
          std::cerr << "Maximum Reference=  " << Maximum_reference[__a] << std::
            endl;
          std::cerr << "Maximum Rel Difference=  " << std::
            abs ((Maximum[__a] -
                  Maximum_reference[__a]) /
                 (Maximum_reference[__a])) << std::endl;
          std::cerr << "Maximum Abs Difference=  " << std::abs (Maximum[__a] -
                                                                Maximum_reference
                                                                [__a]) << std::
            endl;
          return 1;
        }
      for (int __a = 0; __a < 16; __a++)
        if (std::
            abs ((Blend[__a] - Blend_reference[__a]) / (Blend_reference[__a])) >
            1)
        {
          std::cerr << "Mismatch detected in SCALAR implementation" << std::
            endl;
          std::cerr << "Variable Blend:" << std::endl;
          std::cerr << "seed=" << seed << ", __a=" << __a << std::endl;
          std::cerr << "Blend SCALAR=  " << Blend[__a] << std::endl;
          std::cerr << "Blend Reference=  " << Blend_reference[__a] << std::
            endl;
          std::cerr << "Blend Rel Difference=  " << std::
            abs ((Blend[__a] -
                  Blend_reference[__a]) / (Blend_reference[__a])) << std::endl;
          std::cerr << "Blend Abs Difference=  " << std::abs (Blend[__a] -
                                                              Blend_reference
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
      typedef         T (&refArray6)[16];
      typedef         T (&refArray7)[16];
      typedef         T (&refArray8)[16];
      typedef         T (&refArray9)[16];
      typedef         T (&refArray10)[16];
      typedef         T (&refArray11)[16];
      typedef         T (&refArray12)[16];
      typedef         T (&refArray13)[16];
      typedef         T (&refArray14)[16];
      typedef         T (&refArray15)[16];
      typedef         T (&refArray16)[16];
      typedef         T (&refArray17)[16];
      typedef         T (&refArray18)[16];
      typedef         T (&refArray19)[16];
      typedef         T (&refArray20)[16];
      typedef         T (&refArray21)[16];
      typedef         T (&refArray22)[16];
      typedef         T (&refArray23)[16];
      typedef         T (&refArray24)[16];
      typedef         T (&refArray25)[16];
      typedef         T (&refArray26)[16];
      typedef         T (&refArray27)[16];
      typedef         T (&refArray28)[16];
      for (int __a = 0; __a < 16; __a++)
        Add[__a] = Add_original[__a];
      for (int __a = 0; __a < 16; __a++)
        Multiply[__a] = Multiply_original[__a];
      for (int __a = 0; __a < 16; __a++)
        Subtract[__a] = Subtract_original[__a];
      for (int __a = 0; __a < 16; __a++)
        Divide[__a] = Divide_original[__a];
      for (int __a = 0; __a < 16; __a++)
        LessThan[__a] = LessThan_original[__a];
      for (int __a = 0; __a < 16; __a++)
        GreaterThan[__a] = GreaterThan_original[__a];
      for (int __a = 0; __a < 16; __a++)
        LessEquals[__a] = LessEquals_original[__a];
      for (int __a = 0; __a < 16; __a++)
        GreaterEquals[__a] = GreaterEquals_original[__a];
      for (int __a = 0; __a < 16; __a++)
        Equals[__a] = Equals_original[__a];
      for (int __a = 0; __a < 16; __a++)
        And[__a] = And_original[__a];
      for (int __a = 0; __a < 16; __a++)
        Or[__a] = Or_original[__a];
      for (int __a = 0; __a < 16; __a++)
        Xor[__a] = Xor_original[__a];
      for (int __a = 0; __a < 16; __a++)
        AndNot[__a] = AndNot_original[__a];
      for (int __a = 0; __a < 16; __a++)
        Not[__a] = Not_original[__a];
      for (int __a = 0; __a < 16; __a++)
        Sqrt[__a] = Sqrt_original[__a];
      for (int __a = 0; __a < 16; __a++)
        RSqrt[__a] = RSqrt_original[__a];
      for (int __a = 0; __a < 16; __a++)
        Log[__a] = Log_original[__a];
      for (int __a = 0; __a < 16; __a++)
        Exp[__a] = Exp_original[__a];
      for (int __a = 0; __a < 16; __a++)
        Inverse[__a] = Inverse_original[__a];
      for (int __a = 0; __a < 16; __a++)
        AbsoluteValue[__a] = AbsoluteValue_original[__a];
      for (int __a = 0; __a < 16; __a++)
        Sign[__a] = Sign_original[__a];
      for (int __a = 0; __a < 16; __a++)
        Minimum[__a] = Minimum_original[__a];
      for (int __a = 0; __a < 16; __a++)
        Maximum[__a] = Maximum_original[__a];
      for (int __a = 0; __a < 16; __a++)
        Blend[__a] = Blend_original[__a];
      for (int i = 0; i < 16; i += 4)
      {
        refArray1       Ak = reinterpret_cast < refArray1 > (A[i]);
        refArray2       Bk = reinterpret_cast < refArray2 > (B[i]);
        refArray3       Ck = reinterpret_cast < refArray3 > (C[i]);
        refArray4       Dk = reinterpret_cast < refArray4 > (D[i]);
        refArray5       Addk = reinterpret_cast < refArray5 > (Add[i]);
        refArray6       Multiplyk =
          reinterpret_cast < refArray6 > (Multiply[i]);
        refArray7       Subtractk =
          reinterpret_cast < refArray7 > (Subtract[i]);
        refArray8       Dividek = reinterpret_cast < refArray8 > (Divide[i]);
        refArray9       LessThank =
          reinterpret_cast < refArray9 > (LessThan[i]);
        refArray10      GreaterThank =
          reinterpret_cast < refArray10 > (GreaterThan[i]);
        refArray11      LessEqualsk =
          reinterpret_cast < refArray11 > (LessEquals[i]);
        refArray12      GreaterEqualsk =
          reinterpret_cast < refArray12 > (GreaterEquals[i]);
        refArray13      Equalsk = reinterpret_cast < refArray13 > (Equals[i]);
        refArray14      Andk = reinterpret_cast < refArray14 > (And[i]);
        refArray15      Ork = reinterpret_cast < refArray15 > (Or[i]);
        refArray16      Xork = reinterpret_cast < refArray16 > (Xor[i]);
        refArray17      AndNotk = reinterpret_cast < refArray17 > (AndNot[i]);
        refArray18      Notk = reinterpret_cast < refArray18 > (Not[i]);
        refArray19      Sqrtk = reinterpret_cast < refArray19 > (Sqrt[i]);
        refArray20      RSqrtk = reinterpret_cast < refArray20 > (RSqrt[i]);
        refArray21      Logk = reinterpret_cast < refArray21 > (Log[i]);
        refArray22      Expk = reinterpret_cast < refArray22 > (Exp[i]);
        refArray23      Inversek = reinterpret_cast < refArray23 > (Inverse[i]);
        refArray24      AbsoluteValuek =
          reinterpret_cast < refArray24 > (AbsoluteValue[i]);
        refArray25      Signk = reinterpret_cast < refArray25 > (Sign[i]);
        refArray26      Minimumk = reinterpret_cast < refArray26 > (Minimum[i]);
        refArray27      Maximumk = reinterpret_cast < refArray27 > (Maximum[i]);
        refArray28      Blendk = reinterpret_cast < refArray28 > (Blend[i]);
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
      for (int __a = 0; __a < 16; __a++)
        if (std::abs ((Add[__a] - Add_reference[__a]) / (Add_reference[__a])) >
            1)
        {
          std::cerr << "Mismatch detected in SSE implementation" << std::endl;
          std::cerr << "Variable Add:" << std::endl;
          std::cerr << "seed=" << seed << ", __a=" << __a << std::endl;
          std::cerr << "Add SSE=  " << Add[__a] << std::endl;
          std::cerr << "Add Reference=  " << Add_reference[__a] << std::endl;
          std::cerr << "Add Rel Difference=  " << std::
            abs ((Add[__a] -
                  Add_reference[__a]) / (Add_reference[__a])) << std::endl;
          std::cerr << "Add Abs Difference=  " << std::abs (Add[__a] -
                                                            Add_reference[__a])
            << std::endl;
          return 1;
        }
      for (int __a = 0; __a < 16; __a++)
        if (std::
            abs ((Multiply[__a] -
                  Multiply_reference[__a]) / (Multiply_reference[__a])) > 1)
        {
          std::cerr << "Mismatch detected in SSE implementation" << std::endl;
          std::cerr << "Variable Multiply:" << std::endl;
          std::cerr << "seed=" << seed << ", __a=" << __a << std::endl;
          std::cerr << "Multiply SSE=  " << Multiply[__a] << std::endl;
          std::
            cerr << "Multiply Reference=  " << Multiply_reference[__a] << std::
            endl;
          std::cerr << "Multiply Rel Difference=  " << std::
            abs ((Multiply[__a] -
                  Multiply_reference[__a]) /
                 (Multiply_reference[__a])) << std::endl;
          std::cerr << "Multiply Abs Difference=  " << std::abs (Multiply[__a] -
                                                                 Multiply_reference
                                                                 [__a]) << std::
            endl;
          return 1;
        }
      for (int __a = 0; __a < 16; __a++)
        if (std::
            abs ((Subtract[__a] -
                  Subtract_reference[__a]) / (Subtract_reference[__a])) > 1)
        {
          std::cerr << "Mismatch detected in SSE implementation" << std::endl;
          std::cerr << "Variable Subtract:" << std::endl;
          std::cerr << "seed=" << seed << ", __a=" << __a << std::endl;
          std::cerr << "Subtract SSE=  " << Subtract[__a] << std::endl;
          std::
            cerr << "Subtract Reference=  " << Subtract_reference[__a] << std::
            endl;
          std::cerr << "Subtract Rel Difference=  " << std::
            abs ((Subtract[__a] -
                  Subtract_reference[__a]) /
                 (Subtract_reference[__a])) << std::endl;
          std::cerr << "Subtract Abs Difference=  " << std::abs (Subtract[__a] -
                                                                 Subtract_reference
                                                                 [__a]) << std::
            endl;
          return 1;
        }
      for (int __a = 0; __a < 16; __a++)
        if (std::
            abs ((Divide[__a] -
                  Divide_reference[__a]) / (Divide_reference[__a])) > 1)
        {
          std::cerr << "Mismatch detected in SSE implementation" << std::endl;
          std::cerr << "Variable Divide:" << std::endl;
          std::cerr << "seed=" << seed << ", __a=" << __a << std::endl;
          std::cerr << "Divide SSE=  " << Divide[__a] << std::endl;
          std::cerr << "Divide Reference=  " << Divide_reference[__a] << std::
            endl;
          std::cerr << "Divide Rel Difference=  " << std::
            abs ((Divide[__a] -
                  Divide_reference[__a]) /
                 (Divide_reference[__a])) << std::endl;
          std::cerr << "Divide Abs Difference=  " << std::abs (Divide[__a] -
                                                               Divide_reference
                                                               [__a]) << std::
            endl;
          return 1;
        }
      for (int __a = 0; __a < 16; __a++)
        if (std::
            abs ((LessThan[__a] -
                  LessThan_reference[__a]) / (LessThan_reference[__a])) > 1)
        {
          std::cerr << "Mismatch detected in SSE implementation" << std::endl;
          std::cerr << "Variable LessThan:" << std::endl;
          std::cerr << "seed=" << seed << ", __a=" << __a << std::endl;
          std::cerr << "LessThan SSE=  " << LessThan[__a] << std::endl;
          std::
            cerr << "LessThan Reference=  " << LessThan_reference[__a] << std::
            endl;
          std::cerr << "LessThan Rel Difference=  " << std::
            abs ((LessThan[__a] -
                  LessThan_reference[__a]) /
                 (LessThan_reference[__a])) << std::endl;
          std::cerr << "LessThan Abs Difference=  " << std::abs (LessThan[__a] -
                                                                 LessThan_reference
                                                                 [__a]) << std::
            endl;
          return 1;
        }
      for (int __a = 0; __a < 16; __a++)
        if (std::
            abs ((GreaterThan[__a] -
                  GreaterThan_reference[__a]) / (GreaterThan_reference[__a])) >
            1)
        {
          std::cerr << "Mismatch detected in SSE implementation" << std::endl;
          std::cerr << "Variable GreaterThan:" << std::endl;
          std::cerr << "seed=" << seed << ", __a=" << __a << std::endl;
          std::cerr << "GreaterThan SSE=  " << GreaterThan[__a] << std::endl;
          std::
            cerr << "GreaterThan Reference=  " << GreaterThan_reference[__a] <<
            std::endl;
          std::cerr << "GreaterThan Rel Difference=  " << std::
            abs ((GreaterThan[__a] -
                  GreaterThan_reference[__a]) /
                 (GreaterThan_reference[__a])) << std::endl;
          std::cerr << "GreaterThan Abs Difference=  " << std::
            abs (GreaterThan[__a] - GreaterThan_reference[__a]) << std::endl;
          return 1;
        }
      for (int __a = 0; __a < 16; __a++)
        if (std::
            abs ((LessEquals[__a] -
                  LessEquals_reference[__a]) / (LessEquals_reference[__a])) > 1)
        {
          std::cerr << "Mismatch detected in SSE implementation" << std::endl;
          std::cerr << "Variable LessEquals:" << std::endl;
          std::cerr << "seed=" << seed << ", __a=" << __a << std::endl;
          std::cerr << "LessEquals SSE=  " << LessEquals[__a] << std::endl;
          std::
            cerr << "LessEquals Reference=  " << LessEquals_reference[__a] <<
            std::endl;
          std::cerr << "LessEquals Rel Difference=  " << std::
            abs ((LessEquals[__a] -
                  LessEquals_reference[__a]) /
                 (LessEquals_reference[__a])) << std::endl;
          std::cerr << "LessEquals Abs Difference=  " << std::
            abs (LessEquals[__a] - LessEquals_reference[__a]) << std::endl;
          return 1;
        }
      for (int __a = 0; __a < 16; __a++)
        if (std::
            abs ((GreaterEquals[__a] -
                  GreaterEquals_reference[__a]) /
                 (GreaterEquals_reference[__a])) > 1)
        {
          std::cerr << "Mismatch detected in SSE implementation" << std::endl;
          std::cerr << "Variable GreaterEquals:" << std::endl;
          std::cerr << "seed=" << seed << ", __a=" << __a << std::endl;
          std::cerr << "GreaterEquals SSE=  " << GreaterEquals[__a] << std::
            endl;
          std::
            cerr << "GreaterEquals Reference=  " << GreaterEquals_reference[__a]
            << std::endl;
          std::cerr << "GreaterEquals Rel Difference=  " << std::
            abs ((GreaterEquals[__a] -
                  GreaterEquals_reference[__a]) /
                 (GreaterEquals_reference[__a])) << std::endl;
          std::cerr << "GreaterEquals Abs Difference=  " << std::
            abs (GreaterEquals[__a] -
                 GreaterEquals_reference[__a]) << std::endl;
          return 1;
        }
      for (int __a = 0; __a < 16; __a++)
        if (std::
            abs ((Equals[__a] -
                  Equals_reference[__a]) / (Equals_reference[__a])) > 1)
        {
          std::cerr << "Mismatch detected in SSE implementation" << std::endl;
          std::cerr << "Variable Equals:" << std::endl;
          std::cerr << "seed=" << seed << ", __a=" << __a << std::endl;
          std::cerr << "Equals SSE=  " << Equals[__a] << std::endl;
          std::cerr << "Equals Reference=  " << Equals_reference[__a] << std::
            endl;
          std::cerr << "Equals Rel Difference=  " << std::
            abs ((Equals[__a] -
                  Equals_reference[__a]) /
                 (Equals_reference[__a])) << std::endl;
          std::cerr << "Equals Abs Difference=  " << std::abs (Equals[__a] -
                                                               Equals_reference
                                                               [__a]) << std::
            endl;
          return 1;
        }
      for (int __a = 0; __a < 16; __a++)
        if (std::abs ((And[__a] - And_reference[__a]) / (And_reference[__a])) >
            1)
        {
          std::cerr << "Mismatch detected in SSE implementation" << std::endl;
          std::cerr << "Variable And:" << std::endl;
          std::cerr << "seed=" << seed << ", __a=" << __a << std::endl;
          std::cerr << "And SSE=  " << And[__a] << std::endl;
          std::cerr << "And Reference=  " << And_reference[__a] << std::endl;
          std::cerr << "And Rel Difference=  " << std::
            abs ((And[__a] -
                  And_reference[__a]) / (And_reference[__a])) << std::endl;
          std::cerr << "And Abs Difference=  " << std::abs (And[__a] -
                                                            And_reference[__a])
            << std::endl;
          return 1;
        }
      for (int __a = 0; __a < 16; __a++)
        if (std::abs ((Or[__a] - Or_reference[__a]) / (Or_reference[__a])) > 1)
        {
          std::cerr << "Mismatch detected in SSE implementation" << std::endl;
          std::cerr << "Variable Or:" << std::endl;
          std::cerr << "seed=" << seed << ", __a=" << __a << std::endl;
          std::cerr << "Or SSE=  " << Or[__a] << std::endl;
          std::cerr << "Or Reference=  " << Or_reference[__a] << std::endl;
          std::cerr << "Or Rel Difference=  " << std::
            abs ((Or[__a] -
                  Or_reference[__a]) / (Or_reference[__a])) << std::endl;
          std::cerr << "Or Abs Difference=  " << std::abs (Or[__a] -
                                                           Or_reference[__a]) <<
            std::endl;
          return 1;
        }
      for (int __a = 0; __a < 16; __a++)
        if (std::abs ((Xor[__a] - Xor_reference[__a]) / (Xor_reference[__a])) >
            1)
        {
          std::cerr << "Mismatch detected in SSE implementation" << std::endl;
          std::cerr << "Variable Xor:" << std::endl;
          std::cerr << "seed=" << seed << ", __a=" << __a << std::endl;
          std::cerr << "Xor SSE=  " << Xor[__a] << std::endl;
          std::cerr << "Xor Reference=  " << Xor_reference[__a] << std::endl;
          std::cerr << "Xor Rel Difference=  " << std::
            abs ((Xor[__a] -
                  Xor_reference[__a]) / (Xor_reference[__a])) << std::endl;
          std::cerr << "Xor Abs Difference=  " << std::abs (Xor[__a] -
                                                            Xor_reference[__a])
            << std::endl;
          return 1;
        }
      for (int __a = 0; __a < 16; __a++)
        if (std::
            abs ((AndNot[__a] -
                  AndNot_reference[__a]) / (AndNot_reference[__a])) > 1)
        {
          std::cerr << "Mismatch detected in SSE implementation" << std::endl;
          std::cerr << "Variable AndNot:" << std::endl;
          std::cerr << "seed=" << seed << ", __a=" << __a << std::endl;
          std::cerr << "AndNot SSE=  " << AndNot[__a] << std::endl;
          std::cerr << "AndNot Reference=  " << AndNot_reference[__a] << std::
            endl;
          std::cerr << "AndNot Rel Difference=  " << std::
            abs ((AndNot[__a] -
                  AndNot_reference[__a]) /
                 (AndNot_reference[__a])) << std::endl;
          std::cerr << "AndNot Abs Difference=  " << std::abs (AndNot[__a] -
                                                               AndNot_reference
                                                               [__a]) << std::
            endl;
          return 1;
        }
      for (int __a = 0; __a < 16; __a++)
        if (std::abs ((Not[__a] - Not_reference[__a]) / (Not_reference[__a])) >
            1)
        {
          std::cerr << "Mismatch detected in SSE implementation" << std::endl;
          std::cerr << "Variable Not:" << std::endl;
          std::cerr << "seed=" << seed << ", __a=" << __a << std::endl;
          std::cerr << "Not SSE=  " << Not[__a] << std::endl;
          std::cerr << "Not Reference=  " << Not_reference[__a] << std::endl;
          std::cerr << "Not Rel Difference=  " << std::
            abs ((Not[__a] -
                  Not_reference[__a]) / (Not_reference[__a])) << std::endl;
          std::cerr << "Not Abs Difference=  " << std::abs (Not[__a] -
                                                            Not_reference[__a])
            << std::endl;
          return 1;
        }
      for (int __a = 0; __a < 16; __a++)
        if (std::
            abs ((Sqrt[__a] - Sqrt_reference[__a]) / (Sqrt_reference[__a])) > 1)
        {
          std::cerr << "Mismatch detected in SSE implementation" << std::endl;
          std::cerr << "Variable Sqrt:" << std::endl;
          std::cerr << "seed=" << seed << ", __a=" << __a << std::endl;
          std::cerr << "Sqrt SSE=  " << Sqrt[__a] << std::endl;
          std::cerr << "Sqrt Reference=  " << Sqrt_reference[__a] << std::endl;
          std::cerr << "Sqrt Rel Difference=  " << std::
            abs ((Sqrt[__a] -
                  Sqrt_reference[__a]) / (Sqrt_reference[__a])) << std::endl;
          std::cerr << "Sqrt Abs Difference=  " << std::abs (Sqrt[__a] -
                                                             Sqrt_reference
                                                             [__a]) << std::
            endl;
          return 1;
        }
      for (int __a = 0; __a < 16; __a++)
        if (std::
            abs ((RSqrt[__a] - RSqrt_reference[__a]) / (RSqrt_reference[__a])) >
            1)
        {
          std::cerr << "Mismatch detected in SSE implementation" << std::endl;
          std::cerr << "Variable RSqrt:" << std::endl;
          std::cerr << "seed=" << seed << ", __a=" << __a << std::endl;
          std::cerr << "RSqrt SSE=  " << RSqrt[__a] << std::endl;
          std::cerr << "RSqrt Reference=  " << RSqrt_reference[__a] << std::
            endl;
          std::cerr << "RSqrt Rel Difference=  " << std::
            abs ((RSqrt[__a] -
                  RSqrt_reference[__a]) / (RSqrt_reference[__a])) << std::endl;
          std::cerr << "RSqrt Abs Difference=  " << std::abs (RSqrt[__a] -
                                                              RSqrt_reference
                                                              [__a]) << std::
            endl;
          return 1;
        }
      for (int __a = 0; __a < 16; __a++)
        if (std::abs ((Log[__a] - Log_reference[__a]) / (Log_reference[__a])) >
            1)
        {
          std::cerr << "Mismatch detected in SSE implementation" << std::endl;
          std::cerr << "Variable Log:" << std::endl;
          std::cerr << "seed=" << seed << ", __a=" << __a << std::endl;
          std::cerr << "Log SSE=  " << Log[__a] << std::endl;
          std::cerr << "Log Reference=  " << Log_reference[__a] << std::endl;
          std::cerr << "Log Rel Difference=  " << std::
            abs ((Log[__a] -
                  Log_reference[__a]) / (Log_reference[__a])) << std::endl;
          std::cerr << "Log Abs Difference=  " << std::abs (Log[__a] -
                                                            Log_reference[__a])
            << std::endl;
          return 1;
        }
      for (int __a = 0; __a < 16; __a++)
        if (std::abs ((Exp[__a] - Exp_reference[__a]) / (Exp_reference[__a])) >
            1)
        {
          std::cerr << "Mismatch detected in SSE implementation" << std::endl;
          std::cerr << "Variable Exp:" << std::endl;
          std::cerr << "seed=" << seed << ", __a=" << __a << std::endl;
          std::cerr << "Exp SSE=  " << Exp[__a] << std::endl;
          std::cerr << "Exp Reference=  " << Exp_reference[__a] << std::endl;
          std::cerr << "Exp Rel Difference=  " << std::
            abs ((Exp[__a] -
                  Exp_reference[__a]) / (Exp_reference[__a])) << std::endl;
          std::cerr << "Exp Abs Difference=  " << std::abs (Exp[__a] -
                                                            Exp_reference[__a])
            << std::endl;
          return 1;
        }
      for (int __a = 0; __a < 16; __a++)
        if (std::
            abs ((Inverse[__a] -
                  Inverse_reference[__a]) / (Inverse_reference[__a])) > 1)
        {
          std::cerr << "Mismatch detected in SSE implementation" << std::endl;
          std::cerr << "Variable Inverse:" << std::endl;
          std::cerr << "seed=" << seed << ", __a=" << __a << std::endl;
          std::cerr << "Inverse SSE=  " << Inverse[__a] << std::endl;
          std::cerr << "Inverse Reference=  " << Inverse_reference[__a] << std::
            endl;
          std::cerr << "Inverse Rel Difference=  " << std::
            abs ((Inverse[__a] -
                  Inverse_reference[__a]) /
                 (Inverse_reference[__a])) << std::endl;
          std::cerr << "Inverse Abs Difference=  " << std::abs (Inverse[__a] -
                                                                Inverse_reference
                                                                [__a]) << std::
            endl;
          return 1;
        }
      for (int __a = 0; __a < 16; __a++)
        if (std::
            abs ((AbsoluteValue[__a] -
                  AbsoluteValue_reference[__a]) /
                 (AbsoluteValue_reference[__a])) > 1)
        {
          std::cerr << "Mismatch detected in SSE implementation" << std::endl;
          std::cerr << "Variable AbsoluteValue:" << std::endl;
          std::cerr << "seed=" << seed << ", __a=" << __a << std::endl;
          std::cerr << "AbsoluteValue SSE=  " << AbsoluteValue[__a] << std::
            endl;
          std::
            cerr << "AbsoluteValue Reference=  " << AbsoluteValue_reference[__a]
            << std::endl;
          std::cerr << "AbsoluteValue Rel Difference=  " << std::
            abs ((AbsoluteValue[__a] -
                  AbsoluteValue_reference[__a]) /
                 (AbsoluteValue_reference[__a])) << std::endl;
          std::cerr << "AbsoluteValue Abs Difference=  " << std::
            abs (AbsoluteValue[__a] -
                 AbsoluteValue_reference[__a]) << std::endl;
          return 1;
        }
      for (int __a = 0; __a < 16; __a++)
        if (std::
            abs ((Sign[__a] - Sign_reference[__a]) / (Sign_reference[__a])) > 1)
        {
          std::cerr << "Mismatch detected in SSE implementation" << std::endl;
          std::cerr << "Variable Sign:" << std::endl;
          std::cerr << "seed=" << seed << ", __a=" << __a << std::endl;
          std::cerr << "Sign SSE=  " << Sign[__a] << std::endl;
          std::cerr << "Sign Reference=  " << Sign_reference[__a] << std::endl;
          std::cerr << "Sign Rel Difference=  " << std::
            abs ((Sign[__a] -
                  Sign_reference[__a]) / (Sign_reference[__a])) << std::endl;
          std::cerr << "Sign Abs Difference=  " << std::abs (Sign[__a] -
                                                             Sign_reference
                                                             [__a]) << std::
            endl;
          return 1;
        }
      for (int __a = 0; __a < 16; __a++)
        if (std::
            abs ((Minimum[__a] -
                  Minimum_reference[__a]) / (Minimum_reference[__a])) > 1)
        {
          std::cerr << "Mismatch detected in SSE implementation" << std::endl;
          std::cerr << "Variable Minimum:" << std::endl;
          std::cerr << "seed=" << seed << ", __a=" << __a << std::endl;
          std::cerr << "Minimum SSE=  " << Minimum[__a] << std::endl;
          std::cerr << "Minimum Reference=  " << Minimum_reference[__a] << std::
            endl;
          std::cerr << "Minimum Rel Difference=  " << std::
            abs ((Minimum[__a] -
                  Minimum_reference[__a]) /
                 (Minimum_reference[__a])) << std::endl;
          std::cerr << "Minimum Abs Difference=  " << std::abs (Minimum[__a] -
                                                                Minimum_reference
                                                                [__a]) << std::
            endl;
          return 1;
        }
      for (int __a = 0; __a < 16; __a++)
        if (std::
            abs ((Maximum[__a] -
                  Maximum_reference[__a]) / (Maximum_reference[__a])) > 1)
        {
          std::cerr << "Mismatch detected in SSE implementation" << std::endl;
          std::cerr << "Variable Maximum:" << std::endl;
          std::cerr << "seed=" << seed << ", __a=" << __a << std::endl;
          std::cerr << "Maximum SSE=  " << Maximum[__a] << std::endl;
          std::cerr << "Maximum Reference=  " << Maximum_reference[__a] << std::
            endl;
          std::cerr << "Maximum Rel Difference=  " << std::
            abs ((Maximum[__a] -
                  Maximum_reference[__a]) /
                 (Maximum_reference[__a])) << std::endl;
          std::cerr << "Maximum Abs Difference=  " << std::abs (Maximum[__a] -
                                                                Maximum_reference
                                                                [__a]) << std::
            endl;
          return 1;
        }
      for (int __a = 0; __a < 16; __a++)
        if (std::
            abs ((Blend[__a] - Blend_reference[__a]) / (Blend_reference[__a])) >
            1)
        {
          std::cerr << "Mismatch detected in SSE implementation" << std::endl;
          std::cerr << "Variable Blend:" << std::endl;
          std::cerr << "seed=" << seed << ", __a=" << __a << std::endl;
          std::cerr << "Blend SSE=  " << Blend[__a] << std::endl;
          std::cerr << "Blend Reference=  " << Blend_reference[__a] << std::
            endl;
          std::cerr << "Blend Rel Difference=  " << std::
            abs ((Blend[__a] -
                  Blend_reference[__a]) / (Blend_reference[__a])) << std::endl;
          std::cerr << "Blend Abs Difference=  " << std::abs (Blend[__a] -
                                                              Blend_reference
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
      typedef         T (&refArray6)[16];
      typedef         T (&refArray7)[16];
      typedef         T (&refArray8)[16];
      typedef         T (&refArray9)[16];
      typedef         T (&refArray10)[16];
      typedef         T (&refArray11)[16];
      typedef         T (&refArray12)[16];
      typedef         T (&refArray13)[16];
      typedef         T (&refArray14)[16];
      typedef         T (&refArray15)[16];
      typedef         T (&refArray16)[16];
      typedef         T (&refArray17)[16];
      typedef         T (&refArray18)[16];
      typedef         T (&refArray19)[16];
      typedef         T (&refArray20)[16];
      typedef         T (&refArray21)[16];
      typedef         T (&refArray22)[16];
      typedef         T (&refArray23)[16];
      typedef         T (&refArray24)[16];
      typedef         T (&refArray25)[16];
      typedef         T (&refArray26)[16];
      typedef         T (&refArray27)[16];
      typedef         T (&refArray28)[16];
      for (int __a = 0; __a < 16; __a++)
        Add[__a] = Add_original[__a];
      for (int __a = 0; __a < 16; __a++)
        Multiply[__a] = Multiply_original[__a];
      for (int __a = 0; __a < 16; __a++)
        Subtract[__a] = Subtract_original[__a];
      for (int __a = 0; __a < 16; __a++)
        Divide[__a] = Divide_original[__a];
      for (int __a = 0; __a < 16; __a++)
        LessThan[__a] = LessThan_original[__a];
      for (int __a = 0; __a < 16; __a++)
        GreaterThan[__a] = GreaterThan_original[__a];
      for (int __a = 0; __a < 16; __a++)
        LessEquals[__a] = LessEquals_original[__a];
      for (int __a = 0; __a < 16; __a++)
        GreaterEquals[__a] = GreaterEquals_original[__a];
      for (int __a = 0; __a < 16; __a++)
        Equals[__a] = Equals_original[__a];
      for (int __a = 0; __a < 16; __a++)
        And[__a] = And_original[__a];
      for (int __a = 0; __a < 16; __a++)
        Or[__a] = Or_original[__a];
      for (int __a = 0; __a < 16; __a++)
        Xor[__a] = Xor_original[__a];
      for (int __a = 0; __a < 16; __a++)
        AndNot[__a] = AndNot_original[__a];
      for (int __a = 0; __a < 16; __a++)
        Not[__a] = Not_original[__a];
      for (int __a = 0; __a < 16; __a++)
        Sqrt[__a] = Sqrt_original[__a];
      for (int __a = 0; __a < 16; __a++)
        RSqrt[__a] = RSqrt_original[__a];
      for (int __a = 0; __a < 16; __a++)
        Log[__a] = Log_original[__a];
      for (int __a = 0; __a < 16; __a++)
        Exp[__a] = Exp_original[__a];
      for (int __a = 0; __a < 16; __a++)
        Inverse[__a] = Inverse_original[__a];
      for (int __a = 0; __a < 16; __a++)
        AbsoluteValue[__a] = AbsoluteValue_original[__a];
      for (int __a = 0; __a < 16; __a++)
        Sign[__a] = Sign_original[__a];
      for (int __a = 0; __a < 16; __a++)
        Minimum[__a] = Minimum_original[__a];
      for (int __a = 0; __a < 16; __a++)
        Maximum[__a] = Maximum_original[__a];
      for (int __a = 0; __a < 16; __a++)
        Blend[__a] = Blend_original[__a];
      for (int i = 0; i < 16; i += 8)
      {
        refArray1       Ak = reinterpret_cast < refArray1 > (A[i]);
        refArray2       Bk = reinterpret_cast < refArray2 > (B[i]);
        refArray3       Ck = reinterpret_cast < refArray3 > (C[i]);
        refArray4       Dk = reinterpret_cast < refArray4 > (D[i]);
        refArray5       Addk = reinterpret_cast < refArray5 > (Add[i]);
        refArray6       Multiplyk =
          reinterpret_cast < refArray6 > (Multiply[i]);
        refArray7       Subtractk =
          reinterpret_cast < refArray7 > (Subtract[i]);
        refArray8       Dividek = reinterpret_cast < refArray8 > (Divide[i]);
        refArray9       LessThank =
          reinterpret_cast < refArray9 > (LessThan[i]);
        refArray10      GreaterThank =
          reinterpret_cast < refArray10 > (GreaterThan[i]);
        refArray11      LessEqualsk =
          reinterpret_cast < refArray11 > (LessEquals[i]);
        refArray12      GreaterEqualsk =
          reinterpret_cast < refArray12 > (GreaterEquals[i]);
        refArray13      Equalsk = reinterpret_cast < refArray13 > (Equals[i]);
        refArray14      Andk = reinterpret_cast < refArray14 > (And[i]);
        refArray15      Ork = reinterpret_cast < refArray15 > (Or[i]);
        refArray16      Xork = reinterpret_cast < refArray16 > (Xor[i]);
        refArray17      AndNotk = reinterpret_cast < refArray17 > (AndNot[i]);
        refArray18      Notk = reinterpret_cast < refArray18 > (Not[i]);
        refArray19      Sqrtk = reinterpret_cast < refArray19 > (Sqrt[i]);
        refArray20      RSqrtk = reinterpret_cast < refArray20 > (RSqrt[i]);
        refArray21      Logk = reinterpret_cast < refArray21 > (Log[i]);
        refArray22      Expk = reinterpret_cast < refArray22 > (Exp[i]);
        refArray23      Inversek = reinterpret_cast < refArray23 > (Inverse[i]);
        refArray24      AbsoluteValuek =
          reinterpret_cast < refArray24 > (AbsoluteValue[i]);
        refArray25      Signk = reinterpret_cast < refArray25 > (Sign[i]);
        refArray26      Minimumk = reinterpret_cast < refArray26 > (Minimum[i]);
        refArray27      Maximumk = reinterpret_cast < refArray27 > (Maximum[i]);
        refArray28      Blendk = reinterpret_cast < refArray28 > (Blend[i]);
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
      for (int __a = 0; __a < 16; __a++)
        if (std::abs ((Add[__a] - Add_reference[__a]) / (Add_reference[__a])) >
            1)
        {
          std::cerr << "Mismatch detected in AVX implementation" << std::endl;
          std::cerr << "Variable Add:" << std::endl;
          std::cerr << "seed=" << seed << ", __a=" << __a << std::endl;
          std::cerr << "Add AVX=  " << Add[__a] << std::endl;
          std::cerr << "Add Reference=  " << Add_reference[__a] << std::endl;
          std::cerr << "Add Rel Difference=  " << std::
            abs ((Add[__a] -
                  Add_reference[__a]) / (Add_reference[__a])) << std::endl;
          std::cerr << "Add Abs Difference=  " << std::abs (Add[__a] -
                                                            Add_reference[__a])
            << std::endl;
          return 1;
        }
      for (int __a = 0; __a < 16; __a++)
        if (std::
            abs ((Multiply[__a] -
                  Multiply_reference[__a]) / (Multiply_reference[__a])) > 1)
        {
          std::cerr << "Mismatch detected in AVX implementation" << std::endl;
          std::cerr << "Variable Multiply:" << std::endl;
          std::cerr << "seed=" << seed << ", __a=" << __a << std::endl;
          std::cerr << "Multiply AVX=  " << Multiply[__a] << std::endl;
          std::
            cerr << "Multiply Reference=  " << Multiply_reference[__a] << std::
            endl;
          std::cerr << "Multiply Rel Difference=  " << std::
            abs ((Multiply[__a] -
                  Multiply_reference[__a]) /
                 (Multiply_reference[__a])) << std::endl;
          std::cerr << "Multiply Abs Difference=  " << std::abs (Multiply[__a] -
                                                                 Multiply_reference
                                                                 [__a]) << std::
            endl;
          return 1;
        }
      for (int __a = 0; __a < 16; __a++)
        if (std::
            abs ((Subtract[__a] -
                  Subtract_reference[__a]) / (Subtract_reference[__a])) > 1)
        {
          std::cerr << "Mismatch detected in AVX implementation" << std::endl;
          std::cerr << "Variable Subtract:" << std::endl;
          std::cerr << "seed=" << seed << ", __a=" << __a << std::endl;
          std::cerr << "Subtract AVX=  " << Subtract[__a] << std::endl;
          std::
            cerr << "Subtract Reference=  " << Subtract_reference[__a] << std::
            endl;
          std::cerr << "Subtract Rel Difference=  " << std::
            abs ((Subtract[__a] -
                  Subtract_reference[__a]) /
                 (Subtract_reference[__a])) << std::endl;
          std::cerr << "Subtract Abs Difference=  " << std::abs (Subtract[__a] -
                                                                 Subtract_reference
                                                                 [__a]) << std::
            endl;
          return 1;
        }
      for (int __a = 0; __a < 16; __a++)
        if (std::
            abs ((Divide[__a] -
                  Divide_reference[__a]) / (Divide_reference[__a])) > 1)
        {
          std::cerr << "Mismatch detected in AVX implementation" << std::endl;
          std::cerr << "Variable Divide:" << std::endl;
          std::cerr << "seed=" << seed << ", __a=" << __a << std::endl;
          std::cerr << "Divide AVX=  " << Divide[__a] << std::endl;
          std::cerr << "Divide Reference=  " << Divide_reference[__a] << std::
            endl;
          std::cerr << "Divide Rel Difference=  " << std::
            abs ((Divide[__a] -
                  Divide_reference[__a]) /
                 (Divide_reference[__a])) << std::endl;
          std::cerr << "Divide Abs Difference=  " << std::abs (Divide[__a] -
                                                               Divide_reference
                                                               [__a]) << std::
            endl;
          return 1;
        }
      for (int __a = 0; __a < 16; __a++)
        if (std::
            abs ((LessThan[__a] -
                  LessThan_reference[__a]) / (LessThan_reference[__a])) > 1)
        {
          std::cerr << "Mismatch detected in AVX implementation" << std::endl;
          std::cerr << "Variable LessThan:" << std::endl;
          std::cerr << "seed=" << seed << ", __a=" << __a << std::endl;
          std::cerr << "LessThan AVX=  " << LessThan[__a] << std::endl;
          std::
            cerr << "LessThan Reference=  " << LessThan_reference[__a] << std::
            endl;
          std::cerr << "LessThan Rel Difference=  " << std::
            abs ((LessThan[__a] -
                  LessThan_reference[__a]) /
                 (LessThan_reference[__a])) << std::endl;
          std::cerr << "LessThan Abs Difference=  " << std::abs (LessThan[__a] -
                                                                 LessThan_reference
                                                                 [__a]) << std::
            endl;
          return 1;
        }
      for (int __a = 0; __a < 16; __a++)
        if (std::
            abs ((GreaterThan[__a] -
                  GreaterThan_reference[__a]) / (GreaterThan_reference[__a])) >
            1)
        {
          std::cerr << "Mismatch detected in AVX implementation" << std::endl;
          std::cerr << "Variable GreaterThan:" << std::endl;
          std::cerr << "seed=" << seed << ", __a=" << __a << std::endl;
          std::cerr << "GreaterThan AVX=  " << GreaterThan[__a] << std::endl;
          std::
            cerr << "GreaterThan Reference=  " << GreaterThan_reference[__a] <<
            std::endl;
          std::cerr << "GreaterThan Rel Difference=  " << std::
            abs ((GreaterThan[__a] -
                  GreaterThan_reference[__a]) /
                 (GreaterThan_reference[__a])) << std::endl;
          std::cerr << "GreaterThan Abs Difference=  " << std::
            abs (GreaterThan[__a] - GreaterThan_reference[__a]) << std::endl;
          return 1;
        }
      for (int __a = 0; __a < 16; __a++)
        if (std::
            abs ((LessEquals[__a] -
                  LessEquals_reference[__a]) / (LessEquals_reference[__a])) > 1)
        {
          std::cerr << "Mismatch detected in AVX implementation" << std::endl;
          std::cerr << "Variable LessEquals:" << std::endl;
          std::cerr << "seed=" << seed << ", __a=" << __a << std::endl;
          std::cerr << "LessEquals AVX=  " << LessEquals[__a] << std::endl;
          std::
            cerr << "LessEquals Reference=  " << LessEquals_reference[__a] <<
            std::endl;
          std::cerr << "LessEquals Rel Difference=  " << std::
            abs ((LessEquals[__a] -
                  LessEquals_reference[__a]) /
                 (LessEquals_reference[__a])) << std::endl;
          std::cerr << "LessEquals Abs Difference=  " << std::
            abs (LessEquals[__a] - LessEquals_reference[__a]) << std::endl;
          return 1;
        }
      for (int __a = 0; __a < 16; __a++)
        if (std::
            abs ((GreaterEquals[__a] -
                  GreaterEquals_reference[__a]) /
                 (GreaterEquals_reference[__a])) > 1)
        {
          std::cerr << "Mismatch detected in AVX implementation" << std::endl;
          std::cerr << "Variable GreaterEquals:" << std::endl;
          std::cerr << "seed=" << seed << ", __a=" << __a << std::endl;
          std::cerr << "GreaterEquals AVX=  " << GreaterEquals[__a] << std::
            endl;
          std::
            cerr << "GreaterEquals Reference=  " << GreaterEquals_reference[__a]
            << std::endl;
          std::cerr << "GreaterEquals Rel Difference=  " << std::
            abs ((GreaterEquals[__a] -
                  GreaterEquals_reference[__a]) /
                 (GreaterEquals_reference[__a])) << std::endl;
          std::cerr << "GreaterEquals Abs Difference=  " << std::
            abs (GreaterEquals[__a] -
                 GreaterEquals_reference[__a]) << std::endl;
          return 1;
        }
      for (int __a = 0; __a < 16; __a++)
        if (std::
            abs ((Equals[__a] -
                  Equals_reference[__a]) / (Equals_reference[__a])) > 1)
        {
          std::cerr << "Mismatch detected in AVX implementation" << std::endl;
          std::cerr << "Variable Equals:" << std::endl;
          std::cerr << "seed=" << seed << ", __a=" << __a << std::endl;
          std::cerr << "Equals AVX=  " << Equals[__a] << std::endl;
          std::cerr << "Equals Reference=  " << Equals_reference[__a] << std::
            endl;
          std::cerr << "Equals Rel Difference=  " << std::
            abs ((Equals[__a] -
                  Equals_reference[__a]) /
                 (Equals_reference[__a])) << std::endl;
          std::cerr << "Equals Abs Difference=  " << std::abs (Equals[__a] -
                                                               Equals_reference
                                                               [__a]) << std::
            endl;
          return 1;
        }
      for (int __a = 0; __a < 16; __a++)
        if (std::abs ((And[__a] - And_reference[__a]) / (And_reference[__a])) >
            1)
        {
          std::cerr << "Mismatch detected in AVX implementation" << std::endl;
          std::cerr << "Variable And:" << std::endl;
          std::cerr << "seed=" << seed << ", __a=" << __a << std::endl;
          std::cerr << "And AVX=  " << And[__a] << std::endl;
          std::cerr << "And Reference=  " << And_reference[__a] << std::endl;
          std::cerr << "And Rel Difference=  " << std::
            abs ((And[__a] -
                  And_reference[__a]) / (And_reference[__a])) << std::endl;
          std::cerr << "And Abs Difference=  " << std::abs (And[__a] -
                                                            And_reference[__a])
            << std::endl;
          return 1;
        }
      for (int __a = 0; __a < 16; __a++)
        if (std::abs ((Or[__a] - Or_reference[__a]) / (Or_reference[__a])) > 1)
        {
          std::cerr << "Mismatch detected in AVX implementation" << std::endl;
          std::cerr << "Variable Or:" << std::endl;
          std::cerr << "seed=" << seed << ", __a=" << __a << std::endl;
          std::cerr << "Or AVX=  " << Or[__a] << std::endl;
          std::cerr << "Or Reference=  " << Or_reference[__a] << std::endl;
          std::cerr << "Or Rel Difference=  " << std::
            abs ((Or[__a] -
                  Or_reference[__a]) / (Or_reference[__a])) << std::endl;
          std::cerr << "Or Abs Difference=  " << std::abs (Or[__a] -
                                                           Or_reference[__a]) <<
            std::endl;
          return 1;
        }
      for (int __a = 0; __a < 16; __a++)
        if (std::abs ((Xor[__a] - Xor_reference[__a]) / (Xor_reference[__a])) >
            1)
        {
          std::cerr << "Mismatch detected in AVX implementation" << std::endl;
          std::cerr << "Variable Xor:" << std::endl;
          std::cerr << "seed=" << seed << ", __a=" << __a << std::endl;
          std::cerr << "Xor AVX=  " << Xor[__a] << std::endl;
          std::cerr << "Xor Reference=  " << Xor_reference[__a] << std::endl;
          std::cerr << "Xor Rel Difference=  " << std::
            abs ((Xor[__a] -
                  Xor_reference[__a]) / (Xor_reference[__a])) << std::endl;
          std::cerr << "Xor Abs Difference=  " << std::abs (Xor[__a] -
                                                            Xor_reference[__a])
            << std::endl;
          return 1;
        }
      for (int __a = 0; __a < 16; __a++)
        if (std::
            abs ((AndNot[__a] -
                  AndNot_reference[__a]) / (AndNot_reference[__a])) > 1)
        {
          std::cerr << "Mismatch detected in AVX implementation" << std::endl;
          std::cerr << "Variable AndNot:" << std::endl;
          std::cerr << "seed=" << seed << ", __a=" << __a << std::endl;
          std::cerr << "AndNot AVX=  " << AndNot[__a] << std::endl;
          std::cerr << "AndNot Reference=  " << AndNot_reference[__a] << std::
            endl;
          std::cerr << "AndNot Rel Difference=  " << std::
            abs ((AndNot[__a] -
                  AndNot_reference[__a]) /
                 (AndNot_reference[__a])) << std::endl;
          std::cerr << "AndNot Abs Difference=  " << std::abs (AndNot[__a] -
                                                               AndNot_reference
                                                               [__a]) << std::
            endl;
          return 1;
        }
      for (int __a = 0; __a < 16; __a++)
        if (std::abs ((Not[__a] - Not_reference[__a]) / (Not_reference[__a])) >
            1)
        {
          std::cerr << "Mismatch detected in AVX implementation" << std::endl;
          std::cerr << "Variable Not:" << std::endl;
          std::cerr << "seed=" << seed << ", __a=" << __a << std::endl;
          std::cerr << "Not AVX=  " << Not[__a] << std::endl;
          std::cerr << "Not Reference=  " << Not_reference[__a] << std::endl;
          std::cerr << "Not Rel Difference=  " << std::
            abs ((Not[__a] -
                  Not_reference[__a]) / (Not_reference[__a])) << std::endl;
          std::cerr << "Not Abs Difference=  " << std::abs (Not[__a] -
                                                            Not_reference[__a])
            << std::endl;
          return 1;
        }
      for (int __a = 0; __a < 16; __a++)
        if (std::
            abs ((Sqrt[__a] - Sqrt_reference[__a]) / (Sqrt_reference[__a])) > 1)
        {
          std::cerr << "Mismatch detected in AVX implementation" << std::endl;
          std::cerr << "Variable Sqrt:" << std::endl;
          std::cerr << "seed=" << seed << ", __a=" << __a << std::endl;
          std::cerr << "Sqrt AVX=  " << Sqrt[__a] << std::endl;
          std::cerr << "Sqrt Reference=  " << Sqrt_reference[__a] << std::endl;
          std::cerr << "Sqrt Rel Difference=  " << std::
            abs ((Sqrt[__a] -
                  Sqrt_reference[__a]) / (Sqrt_reference[__a])) << std::endl;
          std::cerr << "Sqrt Abs Difference=  " << std::abs (Sqrt[__a] -
                                                             Sqrt_reference
                                                             [__a]) << std::
            endl;
          return 1;
        }
      for (int __a = 0; __a < 16; __a++)
        if (std::
            abs ((RSqrt[__a] - RSqrt_reference[__a]) / (RSqrt_reference[__a])) >
            1)
        {
          std::cerr << "Mismatch detected in AVX implementation" << std::endl;
          std::cerr << "Variable RSqrt:" << std::endl;
          std::cerr << "seed=" << seed << ", __a=" << __a << std::endl;
          std::cerr << "RSqrt AVX=  " << RSqrt[__a] << std::endl;
          std::cerr << "RSqrt Reference=  " << RSqrt_reference[__a] << std::
            endl;
          std::cerr << "RSqrt Rel Difference=  " << std::
            abs ((RSqrt[__a] -
                  RSqrt_reference[__a]) / (RSqrt_reference[__a])) << std::endl;
          std::cerr << "RSqrt Abs Difference=  " << std::abs (RSqrt[__a] -
                                                              RSqrt_reference
                                                              [__a]) << std::
            endl;
          return 1;
        }
      for (int __a = 0; __a < 16; __a++)
        if (std::abs ((Log[__a] - Log_reference[__a]) / (Log_reference[__a])) >
            1)
        {
          std::cerr << "Mismatch detected in AVX implementation" << std::endl;
          std::cerr << "Variable Log:" << std::endl;
          std::cerr << "seed=" << seed << ", __a=" << __a << std::endl;
          std::cerr << "Log AVX=  " << Log[__a] << std::endl;
          std::cerr << "Log Reference=  " << Log_reference[__a] << std::endl;
          std::cerr << "Log Rel Difference=  " << std::
            abs ((Log[__a] -
                  Log_reference[__a]) / (Log_reference[__a])) << std::endl;
          std::cerr << "Log Abs Difference=  " << std::abs (Log[__a] -
                                                            Log_reference[__a])
            << std::endl;
          return 1;
        }
      for (int __a = 0; __a < 16; __a++)
        if (std::abs ((Exp[__a] - Exp_reference[__a]) / (Exp_reference[__a])) >
            1)
        {
          std::cerr << "Mismatch detected in AVX implementation" << std::endl;
          std::cerr << "Variable Exp:" << std::endl;
          std::cerr << "seed=" << seed << ", __a=" << __a << std::endl;
          std::cerr << "Exp AVX=  " << Exp[__a] << std::endl;
          std::cerr << "Exp Reference=  " << Exp_reference[__a] << std::endl;
          std::cerr << "Exp Rel Difference=  " << std::
            abs ((Exp[__a] -
                  Exp_reference[__a]) / (Exp_reference[__a])) << std::endl;
          std::cerr << "Exp Abs Difference=  " << std::abs (Exp[__a] -
                                                            Exp_reference[__a])
            << std::endl;
          return 1;
        }
      for (int __a = 0; __a < 16; __a++)
        if (std::
            abs ((Inverse[__a] -
                  Inverse_reference[__a]) / (Inverse_reference[__a])) > 1)
        {
          std::cerr << "Mismatch detected in AVX implementation" << std::endl;
          std::cerr << "Variable Inverse:" << std::endl;
          std::cerr << "seed=" << seed << ", __a=" << __a << std::endl;
          std::cerr << "Inverse AVX=  " << Inverse[__a] << std::endl;
          std::cerr << "Inverse Reference=  " << Inverse_reference[__a] << std::
            endl;
          std::cerr << "Inverse Rel Difference=  " << std::
            abs ((Inverse[__a] -
                  Inverse_reference[__a]) /
                 (Inverse_reference[__a])) << std::endl;
          std::cerr << "Inverse Abs Difference=  " << std::abs (Inverse[__a] -
                                                                Inverse_reference
                                                                [__a]) << std::
            endl;
          return 1;
        }
      for (int __a = 0; __a < 16; __a++)
        if (std::
            abs ((AbsoluteValue[__a] -
                  AbsoluteValue_reference[__a]) /
                 (AbsoluteValue_reference[__a])) > 1)
        {
          std::cerr << "Mismatch detected in AVX implementation" << std::endl;
          std::cerr << "Variable AbsoluteValue:" << std::endl;
          std::cerr << "seed=" << seed << ", __a=" << __a << std::endl;
          std::cerr << "AbsoluteValue AVX=  " << AbsoluteValue[__a] << std::
            endl;
          std::
            cerr << "AbsoluteValue Reference=  " << AbsoluteValue_reference[__a]
            << std::endl;
          std::cerr << "AbsoluteValue Rel Difference=  " << std::
            abs ((AbsoluteValue[__a] -
                  AbsoluteValue_reference[__a]) /
                 (AbsoluteValue_reference[__a])) << std::endl;
          std::cerr << "AbsoluteValue Abs Difference=  " << std::
            abs (AbsoluteValue[__a] -
                 AbsoluteValue_reference[__a]) << std::endl;
          return 1;
        }
      for (int __a = 0; __a < 16; __a++)
        if (std::
            abs ((Sign[__a] - Sign_reference[__a]) / (Sign_reference[__a])) > 1)
        {
          std::cerr << "Mismatch detected in AVX implementation" << std::endl;
          std::cerr << "Variable Sign:" << std::endl;
          std::cerr << "seed=" << seed << ", __a=" << __a << std::endl;
          std::cerr << "Sign AVX=  " << Sign[__a] << std::endl;
          std::cerr << "Sign Reference=  " << Sign_reference[__a] << std::endl;
          std::cerr << "Sign Rel Difference=  " << std::
            abs ((Sign[__a] -
                  Sign_reference[__a]) / (Sign_reference[__a])) << std::endl;
          std::cerr << "Sign Abs Difference=  " << std::abs (Sign[__a] -
                                                             Sign_reference
                                                             [__a]) << std::
            endl;
          return 1;
        }
      for (int __a = 0; __a < 16; __a++)
        if (std::
            abs ((Minimum[__a] -
                  Minimum_reference[__a]) / (Minimum_reference[__a])) > 1)
        {
          std::cerr << "Mismatch detected in AVX implementation" << std::endl;
          std::cerr << "Variable Minimum:" << std::endl;
          std::cerr << "seed=" << seed << ", __a=" << __a << std::endl;
          std::cerr << "Minimum AVX=  " << Minimum[__a] << std::endl;
          std::cerr << "Minimum Reference=  " << Minimum_reference[__a] << std::
            endl;
          std::cerr << "Minimum Rel Difference=  " << std::
            abs ((Minimum[__a] -
                  Minimum_reference[__a]) /
                 (Minimum_reference[__a])) << std::endl;
          std::cerr << "Minimum Abs Difference=  " << std::abs (Minimum[__a] -
                                                                Minimum_reference
                                                                [__a]) << std::
            endl;
          return 1;
        }
      for (int __a = 0; __a < 16; __a++)
        if (std::
            abs ((Maximum[__a] -
                  Maximum_reference[__a]) / (Maximum_reference[__a])) > 1)
        {
          std::cerr << "Mismatch detected in AVX implementation" << std::endl;
          std::cerr << "Variable Maximum:" << std::endl;
          std::cerr << "seed=" << seed << ", __a=" << __a << std::endl;
          std::cerr << "Maximum AVX=  " << Maximum[__a] << std::endl;
          std::cerr << "Maximum Reference=  " << Maximum_reference[__a] << std::
            endl;
          std::cerr << "Maximum Rel Difference=  " << std::
            abs ((Maximum[__a] -
                  Maximum_reference[__a]) /
                 (Maximum_reference[__a])) << std::endl;
          std::cerr << "Maximum Abs Difference=  " << std::abs (Maximum[__a] -
                                                                Maximum_reference
                                                                [__a]) << std::
            endl;
          return 1;
        }
      for (int __a = 0; __a < 16; __a++)
        if (std::
            abs ((Blend[__a] - Blend_reference[__a]) / (Blend_reference[__a])) >
            1)
        {
          std::cerr << "Mismatch detected in AVX implementation" << std::endl;
          std::cerr << "Variable Blend:" << std::endl;
          std::cerr << "seed=" << seed << ", __a=" << __a << std::endl;
          std::cerr << "Blend AVX=  " << Blend[__a] << std::endl;
          std::cerr << "Blend Reference=  " << Blend_reference[__a] << std::
            endl;
          std::cerr << "Blend Rel Difference=  " << std::
            abs ((Blend[__a] -
                  Blend_reference[__a]) / (Blend_reference[__a])) << std::endl;
          std::cerr << "Blend Abs Difference=  " << std::abs (Blend[__a] -
                                                              Blend_reference
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
      typedef         T (&refArray6)[16];
      typedef         T (&refArray7)[16];
      typedef         T (&refArray8)[16];
      typedef         T (&refArray9)[16];
      typedef         T (&refArray10)[16];
      typedef         T (&refArray11)[16];
      typedef         T (&refArray12)[16];
      typedef         T (&refArray13)[16];
      typedef         T (&refArray14)[16];
      typedef         T (&refArray15)[16];
      typedef         T (&refArray16)[16];
      typedef         T (&refArray17)[16];
      typedef         T (&refArray18)[16];
      typedef         T (&refArray19)[16];
      typedef         T (&refArray20)[16];
      typedef         T (&refArray21)[16];
      typedef         T (&refArray22)[16];
      typedef         T (&refArray23)[16];
      typedef         T (&refArray24)[16];
      typedef         T (&refArray25)[16];
      typedef         T (&refArray26)[16];
      typedef         T (&refArray27)[16];
      typedef         T (&refArray28)[16];
      for (int __a = 0; __a < 16; __a++)
        Add[__a] = Add_original[__a];
      for (int __a = 0; __a < 16; __a++)
        Multiply[__a] = Multiply_original[__a];
      for (int __a = 0; __a < 16; __a++)
        Subtract[__a] = Subtract_original[__a];
      for (int __a = 0; __a < 16; __a++)
        Divide[__a] = Divide_original[__a];
      for (int __a = 0; __a < 16; __a++)
        LessThan[__a] = LessThan_original[__a];
      for (int __a = 0; __a < 16; __a++)
        GreaterThan[__a] = GreaterThan_original[__a];
      for (int __a = 0; __a < 16; __a++)
        LessEquals[__a] = LessEquals_original[__a];
      for (int __a = 0; __a < 16; __a++)
        GreaterEquals[__a] = GreaterEquals_original[__a];
      for (int __a = 0; __a < 16; __a++)
        Equals[__a] = Equals_original[__a];
      for (int __a = 0; __a < 16; __a++)
        And[__a] = And_original[__a];
      for (int __a = 0; __a < 16; __a++)
        Or[__a] = Or_original[__a];
      for (int __a = 0; __a < 16; __a++)
        Xor[__a] = Xor_original[__a];
      for (int __a = 0; __a < 16; __a++)
        AndNot[__a] = AndNot_original[__a];
      for (int __a = 0; __a < 16; __a++)
        Not[__a] = Not_original[__a];
      for (int __a = 0; __a < 16; __a++)
        Sqrt[__a] = Sqrt_original[__a];
      for (int __a = 0; __a < 16; __a++)
        RSqrt[__a] = RSqrt_original[__a];
      for (int __a = 0; __a < 16; __a++)
        Log[__a] = Log_original[__a];
      for (int __a = 0; __a < 16; __a++)
        Exp[__a] = Exp_original[__a];
      for (int __a = 0; __a < 16; __a++)
        Inverse[__a] = Inverse_original[__a];
      for (int __a = 0; __a < 16; __a++)
        AbsoluteValue[__a] = AbsoluteValue_original[__a];
      for (int __a = 0; __a < 16; __a++)
        Sign[__a] = Sign_original[__a];
      for (int __a = 0; __a < 16; __a++)
        Minimum[__a] = Minimum_original[__a];
      for (int __a = 0; __a < 16; __a++)
        Maximum[__a] = Maximum_original[__a];
      for (int __a = 0; __a < 16; __a++)
        Blend[__a] = Blend_original[__a];
      for (int i = 0; i < 16; i += 4)
      {
        refArray1       Ak = reinterpret_cast < refArray1 > (A[i]);
        refArray2       Bk = reinterpret_cast < refArray2 > (B[i]);
        refArray3       Ck = reinterpret_cast < refArray3 > (C[i]);
        refArray4       Dk = reinterpret_cast < refArray4 > (D[i]);
        refArray5       Addk = reinterpret_cast < refArray5 > (Add[i]);
        refArray6       Multiplyk =
          reinterpret_cast < refArray6 > (Multiply[i]);
        refArray7       Subtractk =
          reinterpret_cast < refArray7 > (Subtract[i]);
        refArray8       Dividek = reinterpret_cast < refArray8 > (Divide[i]);
        refArray9       LessThank =
          reinterpret_cast < refArray9 > (LessThan[i]);
        refArray10      GreaterThank =
          reinterpret_cast < refArray10 > (GreaterThan[i]);
        refArray11      LessEqualsk =
          reinterpret_cast < refArray11 > (LessEquals[i]);
        refArray12      GreaterEqualsk =
          reinterpret_cast < refArray12 > (GreaterEquals[i]);
        refArray13      Equalsk = reinterpret_cast < refArray13 > (Equals[i]);
        refArray14      Andk = reinterpret_cast < refArray14 > (And[i]);
        refArray15      Ork = reinterpret_cast < refArray15 > (Or[i]);
        refArray16      Xork = reinterpret_cast < refArray16 > (Xor[i]);
        refArray17      AndNotk = reinterpret_cast < refArray17 > (AndNot[i]);
        refArray18      Notk = reinterpret_cast < refArray18 > (Not[i]);
        refArray19      Sqrtk = reinterpret_cast < refArray19 > (Sqrt[i]);
        refArray20      RSqrtk = reinterpret_cast < refArray20 > (RSqrt[i]);
        refArray21      Logk = reinterpret_cast < refArray21 > (Log[i]);
        refArray22      Expk = reinterpret_cast < refArray22 > (Exp[i]);
        refArray23      Inversek = reinterpret_cast < refArray23 > (Inverse[i]);
        refArray24      AbsoluteValuek =
          reinterpret_cast < refArray24 > (AbsoluteValue[i]);
        refArray25      Signk = reinterpret_cast < refArray25 > (Sign[i]);
        refArray26      Minimumk = reinterpret_cast < refArray26 > (Minimum[i]);
        refArray27      Maximumk = reinterpret_cast < refArray27 > (Maximum[i]);
        refArray28      Blendk = reinterpret_cast < refArray28 > (Blend[i]);
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
      for (int __a = 0; __a < 16; __a++)
        if (std::abs ((Add[__a] - Add_reference[__a]) / (Add_reference[__a])) >
            1)
        {
          std::cerr << "Mismatch detected in NEON implementation" << std::endl;
          std::cerr << "Variable Add:" << std::endl;
          std::cerr << "seed=" << seed << ", __a=" << __a << std::endl;
          std::cerr << "Add NEON=  " << Add[__a] << std::endl;
          std::cerr << "Add Reference=  " << Add_reference[__a] << std::endl;
          std::cerr << "Add Rel Difference=  " << std::
            abs ((Add[__a] -
                  Add_reference[__a]) / (Add_reference[__a])) << std::endl;
          std::cerr << "Add Abs Difference=  " << std::abs (Add[__a] -
                                                            Add_reference[__a])
            << std::endl;
          return 1;
        }
      for (int __a = 0; __a < 16; __a++)
        if (std::
            abs ((Multiply[__a] -
                  Multiply_reference[__a]) / (Multiply_reference[__a])) > 1)
        {
          std::cerr << "Mismatch detected in NEON implementation" << std::endl;
          std::cerr << "Variable Multiply:" << std::endl;
          std::cerr << "seed=" << seed << ", __a=" << __a << std::endl;
          std::cerr << "Multiply NEON=  " << Multiply[__a] << std::endl;
          std::
            cerr << "Multiply Reference=  " << Multiply_reference[__a] << std::
            endl;
          std::cerr << "Multiply Rel Difference=  " << std::
            abs ((Multiply[__a] -
                  Multiply_reference[__a]) /
                 (Multiply_reference[__a])) << std::endl;
          std::cerr << "Multiply Abs Difference=  " << std::abs (Multiply[__a] -
                                                                 Multiply_reference
                                                                 [__a]) << std::
            endl;
          return 1;
        }
      for (int __a = 0; __a < 16; __a++)
        if (std::
            abs ((Subtract[__a] -
                  Subtract_reference[__a]) / (Subtract_reference[__a])) > 1)
        {
          std::cerr << "Mismatch detected in NEON implementation" << std::endl;
          std::cerr << "Variable Subtract:" << std::endl;
          std::cerr << "seed=" << seed << ", __a=" << __a << std::endl;
          std::cerr << "Subtract NEON=  " << Subtract[__a] << std::endl;
          std::
            cerr << "Subtract Reference=  " << Subtract_reference[__a] << std::
            endl;
          std::cerr << "Subtract Rel Difference=  " << std::
            abs ((Subtract[__a] -
                  Subtract_reference[__a]) /
                 (Subtract_reference[__a])) << std::endl;
          std::cerr << "Subtract Abs Difference=  " << std::abs (Subtract[__a] -
                                                                 Subtract_reference
                                                                 [__a]) << std::
            endl;
          return 1;
        }
      for (int __a = 0; __a < 16; __a++)
        if (std::
            abs ((Divide[__a] -
                  Divide_reference[__a]) / (Divide_reference[__a])) > 1)
        {
          std::cerr << "Mismatch detected in NEON implementation" << std::endl;
          std::cerr << "Variable Divide:" << std::endl;
          std::cerr << "seed=" << seed << ", __a=" << __a << std::endl;
          std::cerr << "Divide NEON=  " << Divide[__a] << std::endl;
          std::cerr << "Divide Reference=  " << Divide_reference[__a] << std::
            endl;
          std::cerr << "Divide Rel Difference=  " << std::
            abs ((Divide[__a] -
                  Divide_reference[__a]) /
                 (Divide_reference[__a])) << std::endl;
          std::cerr << "Divide Abs Difference=  " << std::abs (Divide[__a] -
                                                               Divide_reference
                                                               [__a]) << std::
            endl;
          return 1;
        }
      for (int __a = 0; __a < 16; __a++)
        if (std::
            abs ((LessThan[__a] -
                  LessThan_reference[__a]) / (LessThan_reference[__a])) > 1)
        {
          std::cerr << "Mismatch detected in NEON implementation" << std::endl;
          std::cerr << "Variable LessThan:" << std::endl;
          std::cerr << "seed=" << seed << ", __a=" << __a << std::endl;
          std::cerr << "LessThan NEON=  " << LessThan[__a] << std::endl;
          std::
            cerr << "LessThan Reference=  " << LessThan_reference[__a] << std::
            endl;
          std::cerr << "LessThan Rel Difference=  " << std::
            abs ((LessThan[__a] -
                  LessThan_reference[__a]) /
                 (LessThan_reference[__a])) << std::endl;
          std::cerr << "LessThan Abs Difference=  " << std::abs (LessThan[__a] -
                                                                 LessThan_reference
                                                                 [__a]) << std::
            endl;
          return 1;
        }
      for (int __a = 0; __a < 16; __a++)
        if (std::
            abs ((GreaterThan[__a] -
                  GreaterThan_reference[__a]) / (GreaterThan_reference[__a])) >
            1)
        {
          std::cerr << "Mismatch detected in NEON implementation" << std::endl;
          std::cerr << "Variable GreaterThan:" << std::endl;
          std::cerr << "seed=" << seed << ", __a=" << __a << std::endl;
          std::cerr << "GreaterThan NEON=  " << GreaterThan[__a] << std::endl;
          std::
            cerr << "GreaterThan Reference=  " << GreaterThan_reference[__a] <<
            std::endl;
          std::cerr << "GreaterThan Rel Difference=  " << std::
            abs ((GreaterThan[__a] -
                  GreaterThan_reference[__a]) /
                 (GreaterThan_reference[__a])) << std::endl;
          std::cerr << "GreaterThan Abs Difference=  " << std::
            abs (GreaterThan[__a] - GreaterThan_reference[__a]) << std::endl;
          return 1;
        }
      for (int __a = 0; __a < 16; __a++)
        if (std::
            abs ((LessEquals[__a] -
                  LessEquals_reference[__a]) / (LessEquals_reference[__a])) > 1)
        {
          std::cerr << "Mismatch detected in NEON implementation" << std::endl;
          std::cerr << "Variable LessEquals:" << std::endl;
          std::cerr << "seed=" << seed << ", __a=" << __a << std::endl;
          std::cerr << "LessEquals NEON=  " << LessEquals[__a] << std::endl;
          std::
            cerr << "LessEquals Reference=  " << LessEquals_reference[__a] <<
            std::endl;
          std::cerr << "LessEquals Rel Difference=  " << std::
            abs ((LessEquals[__a] -
                  LessEquals_reference[__a]) /
                 (LessEquals_reference[__a])) << std::endl;
          std::cerr << "LessEquals Abs Difference=  " << std::
            abs (LessEquals[__a] - LessEquals_reference[__a]) << std::endl;
          return 1;
        }
      for (int __a = 0; __a < 16; __a++)
        if (std::
            abs ((GreaterEquals[__a] -
                  GreaterEquals_reference[__a]) /
                 (GreaterEquals_reference[__a])) > 1)
        {
          std::cerr << "Mismatch detected in NEON implementation" << std::endl;
          std::cerr << "Variable GreaterEquals:" << std::endl;
          std::cerr << "seed=" << seed << ", __a=" << __a << std::endl;
          std::cerr << "GreaterEquals NEON=  " << GreaterEquals[__a] << std::
            endl;
          std::
            cerr << "GreaterEquals Reference=  " << GreaterEquals_reference[__a]
            << std::endl;
          std::cerr << "GreaterEquals Rel Difference=  " << std::
            abs ((GreaterEquals[__a] -
                  GreaterEquals_reference[__a]) /
                 (GreaterEquals_reference[__a])) << std::endl;
          std::cerr << "GreaterEquals Abs Difference=  " << std::
            abs (GreaterEquals[__a] -
                 GreaterEquals_reference[__a]) << std::endl;
          return 1;
        }
      for (int __a = 0; __a < 16; __a++)
        if (std::
            abs ((Equals[__a] -
                  Equals_reference[__a]) / (Equals_reference[__a])) > 1)
        {
          std::cerr << "Mismatch detected in NEON implementation" << std::endl;
          std::cerr << "Variable Equals:" << std::endl;
          std::cerr << "seed=" << seed << ", __a=" << __a << std::endl;
          std::cerr << "Equals NEON=  " << Equals[__a] << std::endl;
          std::cerr << "Equals Reference=  " << Equals_reference[__a] << std::
            endl;
          std::cerr << "Equals Rel Difference=  " << std::
            abs ((Equals[__a] -
                  Equals_reference[__a]) /
                 (Equals_reference[__a])) << std::endl;
          std::cerr << "Equals Abs Difference=  " << std::abs (Equals[__a] -
                                                               Equals_reference
                                                               [__a]) << std::
            endl;
          return 1;
        }
      for (int __a = 0; __a < 16; __a++)
        if (std::abs ((And[__a] - And_reference[__a]) / (And_reference[__a])) >
            1)
        {
          std::cerr << "Mismatch detected in NEON implementation" << std::endl;
          std::cerr << "Variable And:" << std::endl;
          std::cerr << "seed=" << seed << ", __a=" << __a << std::endl;
          std::cerr << "And NEON=  " << And[__a] << std::endl;
          std::cerr << "And Reference=  " << And_reference[__a] << std::endl;
          std::cerr << "And Rel Difference=  " << std::
            abs ((And[__a] -
                  And_reference[__a]) / (And_reference[__a])) << std::endl;
          std::cerr << "And Abs Difference=  " << std::abs (And[__a] -
                                                            And_reference[__a])
            << std::endl;
          return 1;
        }
      for (int __a = 0; __a < 16; __a++)
        if (std::abs ((Or[__a] - Or_reference[__a]) / (Or_reference[__a])) > 1)
        {
          std::cerr << "Mismatch detected in NEON implementation" << std::endl;
          std::cerr << "Variable Or:" << std::endl;
          std::cerr << "seed=" << seed << ", __a=" << __a << std::endl;
          std::cerr << "Or NEON=  " << Or[__a] << std::endl;
          std::cerr << "Or Reference=  " << Or_reference[__a] << std::endl;
          std::cerr << "Or Rel Difference=  " << std::
            abs ((Or[__a] -
                  Or_reference[__a]) / (Or_reference[__a])) << std::endl;
          std::cerr << "Or Abs Difference=  " << std::abs (Or[__a] -
                                                           Or_reference[__a]) <<
            std::endl;
          return 1;
        }
      for (int __a = 0; __a < 16; __a++)
        if (std::abs ((Xor[__a] - Xor_reference[__a]) / (Xor_reference[__a])) >
            1)
        {
          std::cerr << "Mismatch detected in NEON implementation" << std::endl;
          std::cerr << "Variable Xor:" << std::endl;
          std::cerr << "seed=" << seed << ", __a=" << __a << std::endl;
          std::cerr << "Xor NEON=  " << Xor[__a] << std::endl;
          std::cerr << "Xor Reference=  " << Xor_reference[__a] << std::endl;
          std::cerr << "Xor Rel Difference=  " << std::
            abs ((Xor[__a] -
                  Xor_reference[__a]) / (Xor_reference[__a])) << std::endl;
          std::cerr << "Xor Abs Difference=  " << std::abs (Xor[__a] -
                                                            Xor_reference[__a])
            << std::endl;
          return 1;
        }
      for (int __a = 0; __a < 16; __a++)
        if (std::
            abs ((AndNot[__a] -
                  AndNot_reference[__a]) / (AndNot_reference[__a])) > 1)
        {
          std::cerr << "Mismatch detected in NEON implementation" << std::endl;
          std::cerr << "Variable AndNot:" << std::endl;
          std::cerr << "seed=" << seed << ", __a=" << __a << std::endl;
          std::cerr << "AndNot NEON=  " << AndNot[__a] << std::endl;
          std::cerr << "AndNot Reference=  " << AndNot_reference[__a] << std::
            endl;
          std::cerr << "AndNot Rel Difference=  " << std::
            abs ((AndNot[__a] -
                  AndNot_reference[__a]) /
                 (AndNot_reference[__a])) << std::endl;
          std::cerr << "AndNot Abs Difference=  " << std::abs (AndNot[__a] -
                                                               AndNot_reference
                                                               [__a]) << std::
            endl;
          return 1;
        }
      for (int __a = 0; __a < 16; __a++)
        if (std::abs ((Not[__a] - Not_reference[__a]) / (Not_reference[__a])) >
            1)
        {
          std::cerr << "Mismatch detected in NEON implementation" << std::endl;
          std::cerr << "Variable Not:" << std::endl;
          std::cerr << "seed=" << seed << ", __a=" << __a << std::endl;
          std::cerr << "Not NEON=  " << Not[__a] << std::endl;
          std::cerr << "Not Reference=  " << Not_reference[__a] << std::endl;
          std::cerr << "Not Rel Difference=  " << std::
            abs ((Not[__a] -
                  Not_reference[__a]) / (Not_reference[__a])) << std::endl;
          std::cerr << "Not Abs Difference=  " << std::abs (Not[__a] -
                                                            Not_reference[__a])
            << std::endl;
          return 1;
        }
      for (int __a = 0; __a < 16; __a++)
        if (std::
            abs ((Sqrt[__a] - Sqrt_reference[__a]) / (Sqrt_reference[__a])) > 1)
        {
          std::cerr << "Mismatch detected in NEON implementation" << std::endl;
          std::cerr << "Variable Sqrt:" << std::endl;
          std::cerr << "seed=" << seed << ", __a=" << __a << std::endl;
          std::cerr << "Sqrt NEON=  " << Sqrt[__a] << std::endl;
          std::cerr << "Sqrt Reference=  " << Sqrt_reference[__a] << std::endl;
          std::cerr << "Sqrt Rel Difference=  " << std::
            abs ((Sqrt[__a] -
                  Sqrt_reference[__a]) / (Sqrt_reference[__a])) << std::endl;
          std::cerr << "Sqrt Abs Difference=  " << std::abs (Sqrt[__a] -
                                                             Sqrt_reference
                                                             [__a]) << std::
            endl;
          return 1;
        }
      for (int __a = 0; __a < 16; __a++)
        if (std::
            abs ((RSqrt[__a] - RSqrt_reference[__a]) / (RSqrt_reference[__a])) >
            1)
        {
          std::cerr << "Mismatch detected in NEON implementation" << std::endl;
          std::cerr << "Variable RSqrt:" << std::endl;
          std::cerr << "seed=" << seed << ", __a=" << __a << std::endl;
          std::cerr << "RSqrt NEON=  " << RSqrt[__a] << std::endl;
          std::cerr << "RSqrt Reference=  " << RSqrt_reference[__a] << std::
            endl;
          std::cerr << "RSqrt Rel Difference=  " << std::
            abs ((RSqrt[__a] -
                  RSqrt_reference[__a]) / (RSqrt_reference[__a])) << std::endl;
          std::cerr << "RSqrt Abs Difference=  " << std::abs (RSqrt[__a] -
                                                              RSqrt_reference
                                                              [__a]) << std::
            endl;
          return 1;
        }
      for (int __a = 0; __a < 16; __a++)
        if (std::abs ((Log[__a] - Log_reference[__a]) / (Log_reference[__a])) >
            1)
        {
          std::cerr << "Mismatch detected in NEON implementation" << std::endl;
          std::cerr << "Variable Log:" << std::endl;
          std::cerr << "seed=" << seed << ", __a=" << __a << std::endl;
          std::cerr << "Log NEON=  " << Log[__a] << std::endl;
          std::cerr << "Log Reference=  " << Log_reference[__a] << std::endl;
          std::cerr << "Log Rel Difference=  " << std::
            abs ((Log[__a] -
                  Log_reference[__a]) / (Log_reference[__a])) << std::endl;
          std::cerr << "Log Abs Difference=  " << std::abs (Log[__a] -
                                                            Log_reference[__a])
            << std::endl;
          return 1;
        }
      for (int __a = 0; __a < 16; __a++)
        if (std::abs ((Exp[__a] - Exp_reference[__a]) / (Exp_reference[__a])) >
            1)
        {
          std::cerr << "Mismatch detected in NEON implementation" << std::endl;
          std::cerr << "Variable Exp:" << std::endl;
          std::cerr << "seed=" << seed << ", __a=" << __a << std::endl;
          std::cerr << "Exp NEON=  " << Exp[__a] << std::endl;
          std::cerr << "Exp Reference=  " << Exp_reference[__a] << std::endl;
          std::cerr << "Exp Rel Difference=  " << std::
            abs ((Exp[__a] -
                  Exp_reference[__a]) / (Exp_reference[__a])) << std::endl;
          std::cerr << "Exp Abs Difference=  " << std::abs (Exp[__a] -
                                                            Exp_reference[__a])
            << std::endl;
          return 1;
        }
      for (int __a = 0; __a < 16; __a++)
        if (std::
            abs ((Inverse[__a] -
                  Inverse_reference[__a]) / (Inverse_reference[__a])) > 1)
        {
          std::cerr << "Mismatch detected in NEON implementation" << std::endl;
          std::cerr << "Variable Inverse:" << std::endl;
          std::cerr << "seed=" << seed << ", __a=" << __a << std::endl;
          std::cerr << "Inverse NEON=  " << Inverse[__a] << std::endl;
          std::cerr << "Inverse Reference=  " << Inverse_reference[__a] << std::
            endl;
          std::cerr << "Inverse Rel Difference=  " << std::
            abs ((Inverse[__a] -
                  Inverse_reference[__a]) /
                 (Inverse_reference[__a])) << std::endl;
          std::cerr << "Inverse Abs Difference=  " << std::abs (Inverse[__a] -
                                                                Inverse_reference
                                                                [__a]) << std::
            endl;
          return 1;
        }
      for (int __a = 0; __a < 16; __a++)
        if (std::
            abs ((AbsoluteValue[__a] -
                  AbsoluteValue_reference[__a]) /
                 (AbsoluteValue_reference[__a])) > 1)
        {
          std::cerr << "Mismatch detected in NEON implementation" << std::endl;
          std::cerr << "Variable AbsoluteValue:" << std::endl;
          std::cerr << "seed=" << seed << ", __a=" << __a << std::endl;
          std::cerr << "AbsoluteValue NEON=  " << AbsoluteValue[__a] << std::
            endl;
          std::
            cerr << "AbsoluteValue Reference=  " << AbsoluteValue_reference[__a]
            << std::endl;
          std::cerr << "AbsoluteValue Rel Difference=  " << std::
            abs ((AbsoluteValue[__a] -
                  AbsoluteValue_reference[__a]) /
                 (AbsoluteValue_reference[__a])) << std::endl;
          std::cerr << "AbsoluteValue Abs Difference=  " << std::
            abs (AbsoluteValue[__a] -
                 AbsoluteValue_reference[__a]) << std::endl;
          return 1;
        }
      for (int __a = 0; __a < 16; __a++)
        if (std::
            abs ((Sign[__a] - Sign_reference[__a]) / (Sign_reference[__a])) > 1)
        {
          std::cerr << "Mismatch detected in NEON implementation" << std::endl;
          std::cerr << "Variable Sign:" << std::endl;
          std::cerr << "seed=" << seed << ", __a=" << __a << std::endl;
          std::cerr << "Sign NEON=  " << Sign[__a] << std::endl;
          std::cerr << "Sign Reference=  " << Sign_reference[__a] << std::endl;
          std::cerr << "Sign Rel Difference=  " << std::
            abs ((Sign[__a] -
                  Sign_reference[__a]) / (Sign_reference[__a])) << std::endl;
          std::cerr << "Sign Abs Difference=  " << std::abs (Sign[__a] -
                                                             Sign_reference
                                                             [__a]) << std::
            endl;
          return 1;
        }
      for (int __a = 0; __a < 16; __a++)
        if (std::
            abs ((Minimum[__a] -
                  Minimum_reference[__a]) / (Minimum_reference[__a])) > 1)
        {
          std::cerr << "Mismatch detected in NEON implementation" << std::endl;
          std::cerr << "Variable Minimum:" << std::endl;
          std::cerr << "seed=" << seed << ", __a=" << __a << std::endl;
          std::cerr << "Minimum NEON=  " << Minimum[__a] << std::endl;
          std::cerr << "Minimum Reference=  " << Minimum_reference[__a] << std::
            endl;
          std::cerr << "Minimum Rel Difference=  " << std::
            abs ((Minimum[__a] -
                  Minimum_reference[__a]) /
                 (Minimum_reference[__a])) << std::endl;
          std::cerr << "Minimum Abs Difference=  " << std::abs (Minimum[__a] -
                                                                Minimum_reference
                                                                [__a]) << std::
            endl;
          return 1;
        }
      for (int __a = 0; __a < 16; __a++)
        if (std::
            abs ((Maximum[__a] -
                  Maximum_reference[__a]) / (Maximum_reference[__a])) > 1)
        {
          std::cerr << "Mismatch detected in NEON implementation" << std::endl;
          std::cerr << "Variable Maximum:" << std::endl;
          std::cerr << "seed=" << seed << ", __a=" << __a << std::endl;
          std::cerr << "Maximum NEON=  " << Maximum[__a] << std::endl;
          std::cerr << "Maximum Reference=  " << Maximum_reference[__a] << std::
            endl;
          std::cerr << "Maximum Rel Difference=  " << std::
            abs ((Maximum[__a] -
                  Maximum_reference[__a]) /
                 (Maximum_reference[__a])) << std::endl;
          std::cerr << "Maximum Abs Difference=  " << std::abs (Maximum[__a] -
                                                                Maximum_reference
                                                                [__a]) << std::
            endl;
          return 1;
        }
      for (int __a = 0; __a < 16; __a++)
        if (std::
            abs ((Blend[__a] - Blend_reference[__a]) / (Blend_reference[__a])) >
            1)
        {
          std::cerr << "Mismatch detected in NEON implementation" << std::endl;
          std::cerr << "Variable Blend:" << std::endl;
          std::cerr << "seed=" << seed << ", __a=" << __a << std::endl;
          std::cerr << "Blend NEON=  " << Blend[__a] << std::endl;
          std::cerr << "Blend Reference=  " << Blend_reference[__a] << std::
            endl;
          std::cerr << "Blend Rel Difference=  " << std::
            abs ((Blend[__a] -
                  Blend_reference[__a]) / (Blend_reference[__a])) << std::endl;
          std::cerr << "Blend Abs Difference=  " << std::abs (Blend[__a] -
                                                              Blend_reference
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
      typedef         T (&refArray6)[16];
      typedef         T (&refArray7)[16];
      typedef         T (&refArray8)[16];
      typedef         T (&refArray9)[16];
      typedef         T (&refArray10)[16];
      typedef         T (&refArray11)[16];
      typedef         T (&refArray12)[16];
      typedef         T (&refArray13)[16];
      typedef         T (&refArray14)[16];
      typedef         T (&refArray15)[16];
      typedef         T (&refArray16)[16];
      typedef         T (&refArray17)[16];
      typedef         T (&refArray18)[16];
      typedef         T (&refArray19)[16];
      typedef         T (&refArray20)[16];
      typedef         T (&refArray21)[16];
      typedef         T (&refArray22)[16];
      typedef         T (&refArray23)[16];
      typedef         T (&refArray24)[16];
      typedef         T (&refArray25)[16];
      typedef         T (&refArray26)[16];
      typedef         T (&refArray27)[16];
      typedef         T (&refArray28)[16];
      for (int __a = 0; __a < 16; __a++)
        Add[__a] = Add_original[__a];
      for (int __a = 0; __a < 16; __a++)
        Multiply[__a] = Multiply_original[__a];
      for (int __a = 0; __a < 16; __a++)
        Subtract[__a] = Subtract_original[__a];
      for (int __a = 0; __a < 16; __a++)
        Divide[__a] = Divide_original[__a];
      for (int __a = 0; __a < 16; __a++)
        LessThan[__a] = LessThan_original[__a];
      for (int __a = 0; __a < 16; __a++)
        GreaterThan[__a] = GreaterThan_original[__a];
      for (int __a = 0; __a < 16; __a++)
        LessEquals[__a] = LessEquals_original[__a];
      for (int __a = 0; __a < 16; __a++)
        GreaterEquals[__a] = GreaterEquals_original[__a];
      for (int __a = 0; __a < 16; __a++)
        Equals[__a] = Equals_original[__a];
      for (int __a = 0; __a < 16; __a++)
        And[__a] = And_original[__a];
      for (int __a = 0; __a < 16; __a++)
        Or[__a] = Or_original[__a];
      for (int __a = 0; __a < 16; __a++)
        Xor[__a] = Xor_original[__a];
      for (int __a = 0; __a < 16; __a++)
        AndNot[__a] = AndNot_original[__a];
      for (int __a = 0; __a < 16; __a++)
        Not[__a] = Not_original[__a];
      for (int __a = 0; __a < 16; __a++)
        Sqrt[__a] = Sqrt_original[__a];
      for (int __a = 0; __a < 16; __a++)
        RSqrt[__a] = RSqrt_original[__a];
      for (int __a = 0; __a < 16; __a++)
        Log[__a] = Log_original[__a];
      for (int __a = 0; __a < 16; __a++)
        Exp[__a] = Exp_original[__a];
      for (int __a = 0; __a < 16; __a++)
        Inverse[__a] = Inverse_original[__a];
      for (int __a = 0; __a < 16; __a++)
        AbsoluteValue[__a] = AbsoluteValue_original[__a];
      for (int __a = 0; __a < 16; __a++)
        Sign[__a] = Sign_original[__a];
      for (int __a = 0; __a < 16; __a++)
        Minimum[__a] = Minimum_original[__a];
      for (int __a = 0; __a < 16; __a++)
        Maximum[__a] = Maximum_original[__a];
      for (int __a = 0; __a < 16; __a++)
        Blend[__a] = Blend_original[__a];
      for (int i = 0; i < 16; i += 16)
      {
        refArray1       Ak = reinterpret_cast < refArray1 > (A[i]);
        refArray2       Bk = reinterpret_cast < refArray2 > (B[i]);
        refArray3       Ck = reinterpret_cast < refArray3 > (C[i]);
        refArray4       Dk = reinterpret_cast < refArray4 > (D[i]);
        refArray5       Addk = reinterpret_cast < refArray5 > (Add[i]);
        refArray6       Multiplyk =
          reinterpret_cast < refArray6 > (Multiply[i]);
        refArray7       Subtractk =
          reinterpret_cast < refArray7 > (Subtract[i]);
        refArray8       Dividek = reinterpret_cast < refArray8 > (Divide[i]);
        refArray9       LessThank =
          reinterpret_cast < refArray9 > (LessThan[i]);
        refArray10      GreaterThank =
          reinterpret_cast < refArray10 > (GreaterThan[i]);
        refArray11      LessEqualsk =
          reinterpret_cast < refArray11 > (LessEquals[i]);
        refArray12      GreaterEqualsk =
          reinterpret_cast < refArray12 > (GreaterEquals[i]);
        refArray13      Equalsk = reinterpret_cast < refArray13 > (Equals[i]);
        refArray14      Andk = reinterpret_cast < refArray14 > (And[i]);
        refArray15      Ork = reinterpret_cast < refArray15 > (Or[i]);
        refArray16      Xork = reinterpret_cast < refArray16 > (Xor[i]);
        refArray17      AndNotk = reinterpret_cast < refArray17 > (AndNot[i]);
        refArray18      Notk = reinterpret_cast < refArray18 > (Not[i]);
        refArray19      Sqrtk = reinterpret_cast < refArray19 > (Sqrt[i]);
        refArray20      RSqrtk = reinterpret_cast < refArray20 > (RSqrt[i]);
        refArray21      Logk = reinterpret_cast < refArray21 > (Log[i]);
        refArray22      Expk = reinterpret_cast < refArray22 > (Exp[i]);
        refArray23      Inversek = reinterpret_cast < refArray23 > (Inverse[i]);
        refArray24      AbsoluteValuek =
          reinterpret_cast < refArray24 > (AbsoluteValue[i]);
        refArray25      Signk = reinterpret_cast < refArray25 > (Sign[i]);
        refArray26      Minimumk = reinterpret_cast < refArray26 > (Minimum[i]);
        refArray27      Maximumk = reinterpret_cast < refArray27 > (Maximum[i]);
        refArray28      Blendk = reinterpret_cast < refArray28 > (Blend[i]);
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
      for (int __a = 0; __a < 16; __a++)
        if (std::abs ((Add[__a] - Add_reference[__a]) / (Add_reference[__a])) >
            1)
        {
          std::cerr << "Mismatch detected in MIC implementation" << std::endl;
          std::cerr << "Variable Add:" << std::endl;
          std::cerr << "seed=" << seed << ", __a=" << __a << std::endl;
          std::cerr << "Add MIC=  " << Add[__a] << std::endl;
          std::cerr << "Add Reference=  " << Add_reference[__a] << std::endl;
          std::cerr << "Add Rel Difference=  " << std::
            abs ((Add[__a] -
                  Add_reference[__a]) / (Add_reference[__a])) << std::endl;
          std::cerr << "Add Abs Difference=  " << std::abs (Add[__a] -
                                                            Add_reference[__a])
            << std::endl;
          return 1;
        }
      for (int __a = 0; __a < 16; __a++)
        if (std::
            abs ((Multiply[__a] -
                  Multiply_reference[__a]) / (Multiply_reference[__a])) > 1)
        {
          std::cerr << "Mismatch detected in MIC implementation" << std::endl;
          std::cerr << "Variable Multiply:" << std::endl;
          std::cerr << "seed=" << seed << ", __a=" << __a << std::endl;
          std::cerr << "Multiply MIC=  " << Multiply[__a] << std::endl;
          std::
            cerr << "Multiply Reference=  " << Multiply_reference[__a] << std::
            endl;
          std::cerr << "Multiply Rel Difference=  " << std::
            abs ((Multiply[__a] -
                  Multiply_reference[__a]) /
                 (Multiply_reference[__a])) << std::endl;
          std::cerr << "Multiply Abs Difference=  " << std::abs (Multiply[__a] -
                                                                 Multiply_reference
                                                                 [__a]) << std::
            endl;
          return 1;
        }
      for (int __a = 0; __a < 16; __a++)
        if (std::
            abs ((Subtract[__a] -
                  Subtract_reference[__a]) / (Subtract_reference[__a])) > 1)
        {
          std::cerr << "Mismatch detected in MIC implementation" << std::endl;
          std::cerr << "Variable Subtract:" << std::endl;
          std::cerr << "seed=" << seed << ", __a=" << __a << std::endl;
          std::cerr << "Subtract MIC=  " << Subtract[__a] << std::endl;
          std::
            cerr << "Subtract Reference=  " << Subtract_reference[__a] << std::
            endl;
          std::cerr << "Subtract Rel Difference=  " << std::
            abs ((Subtract[__a] -
                  Subtract_reference[__a]) /
                 (Subtract_reference[__a])) << std::endl;
          std::cerr << "Subtract Abs Difference=  " << std::abs (Subtract[__a] -
                                                                 Subtract_reference
                                                                 [__a]) << std::
            endl;
          return 1;
        }
      for (int __a = 0; __a < 16; __a++)
        if (std::
            abs ((Divide[__a] -
                  Divide_reference[__a]) / (Divide_reference[__a])) > 1)
        {
          std::cerr << "Mismatch detected in MIC implementation" << std::endl;
          std::cerr << "Variable Divide:" << std::endl;
          std::cerr << "seed=" << seed << ", __a=" << __a << std::endl;
          std::cerr << "Divide MIC=  " << Divide[__a] << std::endl;
          std::cerr << "Divide Reference=  " << Divide_reference[__a] << std::
            endl;
          std::cerr << "Divide Rel Difference=  " << std::
            abs ((Divide[__a] -
                  Divide_reference[__a]) /
                 (Divide_reference[__a])) << std::endl;
          std::cerr << "Divide Abs Difference=  " << std::abs (Divide[__a] -
                                                               Divide_reference
                                                               [__a]) << std::
            endl;
          return 1;
        }
      for (int __a = 0; __a < 16; __a++)
        if (std::
            abs ((LessThan[__a] -
                  LessThan_reference[__a]) / (LessThan_reference[__a])) > 1)
        {
          std::cerr << "Mismatch detected in MIC implementation" << std::endl;
          std::cerr << "Variable LessThan:" << std::endl;
          std::cerr << "seed=" << seed << ", __a=" << __a << std::endl;
          std::cerr << "LessThan MIC=  " << LessThan[__a] << std::endl;
          std::
            cerr << "LessThan Reference=  " << LessThan_reference[__a] << std::
            endl;
          std::cerr << "LessThan Rel Difference=  " << std::
            abs ((LessThan[__a] -
                  LessThan_reference[__a]) /
                 (LessThan_reference[__a])) << std::endl;
          std::cerr << "LessThan Abs Difference=  " << std::abs (LessThan[__a] -
                                                                 LessThan_reference
                                                                 [__a]) << std::
            endl;
          return 1;
        }
      for (int __a = 0; __a < 16; __a++)
        if (std::
            abs ((GreaterThan[__a] -
                  GreaterThan_reference[__a]) / (GreaterThan_reference[__a])) >
            1)
        {
          std::cerr << "Mismatch detected in MIC implementation" << std::endl;
          std::cerr << "Variable GreaterThan:" << std::endl;
          std::cerr << "seed=" << seed << ", __a=" << __a << std::endl;
          std::cerr << "GreaterThan MIC=  " << GreaterThan[__a] << std::endl;
          std::
            cerr << "GreaterThan Reference=  " << GreaterThan_reference[__a] <<
            std::endl;
          std::cerr << "GreaterThan Rel Difference=  " << std::
            abs ((GreaterThan[__a] -
                  GreaterThan_reference[__a]) /
                 (GreaterThan_reference[__a])) << std::endl;
          std::cerr << "GreaterThan Abs Difference=  " << std::
            abs (GreaterThan[__a] - GreaterThan_reference[__a]) << std::endl;
          return 1;
        }
      for (int __a = 0; __a < 16; __a++)
        if (std::
            abs ((LessEquals[__a] -
                  LessEquals_reference[__a]) / (LessEquals_reference[__a])) > 1)
        {
          std::cerr << "Mismatch detected in MIC implementation" << std::endl;
          std::cerr << "Variable LessEquals:" << std::endl;
          std::cerr << "seed=" << seed << ", __a=" << __a << std::endl;
          std::cerr << "LessEquals MIC=  " << LessEquals[__a] << std::endl;
          std::
            cerr << "LessEquals Reference=  " << LessEquals_reference[__a] <<
            std::endl;
          std::cerr << "LessEquals Rel Difference=  " << std::
            abs ((LessEquals[__a] -
                  LessEquals_reference[__a]) /
                 (LessEquals_reference[__a])) << std::endl;
          std::cerr << "LessEquals Abs Difference=  " << std::
            abs (LessEquals[__a] - LessEquals_reference[__a]) << std::endl;
          return 1;
        }
      for (int __a = 0; __a < 16; __a++)
        if (std::
            abs ((GreaterEquals[__a] -
                  GreaterEquals_reference[__a]) /
                 (GreaterEquals_reference[__a])) > 1)
        {
          std::cerr << "Mismatch detected in MIC implementation" << std::endl;
          std::cerr << "Variable GreaterEquals:" << std::endl;
          std::cerr << "seed=" << seed << ", __a=" << __a << std::endl;
          std::cerr << "GreaterEquals MIC=  " << GreaterEquals[__a] << std::
            endl;
          std::
            cerr << "GreaterEquals Reference=  " << GreaterEquals_reference[__a]
            << std::endl;
          std::cerr << "GreaterEquals Rel Difference=  " << std::
            abs ((GreaterEquals[__a] -
                  GreaterEquals_reference[__a]) /
                 (GreaterEquals_reference[__a])) << std::endl;
          std::cerr << "GreaterEquals Abs Difference=  " << std::
            abs (GreaterEquals[__a] -
                 GreaterEquals_reference[__a]) << std::endl;
          return 1;
        }
      for (int __a = 0; __a < 16; __a++)
        if (std::
            abs ((Equals[__a] -
                  Equals_reference[__a]) / (Equals_reference[__a])) > 1)
        {
          std::cerr << "Mismatch detected in MIC implementation" << std::endl;
          std::cerr << "Variable Equals:" << std::endl;
          std::cerr << "seed=" << seed << ", __a=" << __a << std::endl;
          std::cerr << "Equals MIC=  " << Equals[__a] << std::endl;
          std::cerr << "Equals Reference=  " << Equals_reference[__a] << std::
            endl;
          std::cerr << "Equals Rel Difference=  " << std::
            abs ((Equals[__a] -
                  Equals_reference[__a]) /
                 (Equals_reference[__a])) << std::endl;
          std::cerr << "Equals Abs Difference=  " << std::abs (Equals[__a] -
                                                               Equals_reference
                                                               [__a]) << std::
            endl;
          return 1;
        }
      for (int __a = 0; __a < 16; __a++)
        if (std::abs ((And[__a] - And_reference[__a]) / (And_reference[__a])) >
            1)
        {
          std::cerr << "Mismatch detected in MIC implementation" << std::endl;
          std::cerr << "Variable And:" << std::endl;
          std::cerr << "seed=" << seed << ", __a=" << __a << std::endl;
          std::cerr << "And MIC=  " << And[__a] << std::endl;
          std::cerr << "And Reference=  " << And_reference[__a] << std::endl;
          std::cerr << "And Rel Difference=  " << std::
            abs ((And[__a] -
                  And_reference[__a]) / (And_reference[__a])) << std::endl;
          std::cerr << "And Abs Difference=  " << std::abs (And[__a] -
                                                            And_reference[__a])
            << std::endl;
          return 1;
        }
      for (int __a = 0; __a < 16; __a++)
        if (std::abs ((Or[__a] - Or_reference[__a]) / (Or_reference[__a])) > 1)
        {
          std::cerr << "Mismatch detected in MIC implementation" << std::endl;
          std::cerr << "Variable Or:" << std::endl;
          std::cerr << "seed=" << seed << ", __a=" << __a << std::endl;
          std::cerr << "Or MIC=  " << Or[__a] << std::endl;
          std::cerr << "Or Reference=  " << Or_reference[__a] << std::endl;
          std::cerr << "Or Rel Difference=  " << std::
            abs ((Or[__a] -
                  Or_reference[__a]) / (Or_reference[__a])) << std::endl;
          std::cerr << "Or Abs Difference=  " << std::abs (Or[__a] -
                                                           Or_reference[__a]) <<
            std::endl;
          return 1;
        }
      for (int __a = 0; __a < 16; __a++)
        if (std::abs ((Xor[__a] - Xor_reference[__a]) / (Xor_reference[__a])) >
            1)
        {
          std::cerr << "Mismatch detected in MIC implementation" << std::endl;
          std::cerr << "Variable Xor:" << std::endl;
          std::cerr << "seed=" << seed << ", __a=" << __a << std::endl;
          std::cerr << "Xor MIC=  " << Xor[__a] << std::endl;
          std::cerr << "Xor Reference=  " << Xor_reference[__a] << std::endl;
          std::cerr << "Xor Rel Difference=  " << std::
            abs ((Xor[__a] -
                  Xor_reference[__a]) / (Xor_reference[__a])) << std::endl;
          std::cerr << "Xor Abs Difference=  " << std::abs (Xor[__a] -
                                                            Xor_reference[__a])
            << std::endl;
          return 1;
        }
      for (int __a = 0; __a < 16; __a++)
        if (std::
            abs ((AndNot[__a] -
                  AndNot_reference[__a]) / (AndNot_reference[__a])) > 1)
        {
          std::cerr << "Mismatch detected in MIC implementation" << std::endl;
          std::cerr << "Variable AndNot:" << std::endl;
          std::cerr << "seed=" << seed << ", __a=" << __a << std::endl;
          std::cerr << "AndNot MIC=  " << AndNot[__a] << std::endl;
          std::cerr << "AndNot Reference=  " << AndNot_reference[__a] << std::
            endl;
          std::cerr << "AndNot Rel Difference=  " << std::
            abs ((AndNot[__a] -
                  AndNot_reference[__a]) /
                 (AndNot_reference[__a])) << std::endl;
          std::cerr << "AndNot Abs Difference=  " << std::abs (AndNot[__a] -
                                                               AndNot_reference
                                                               [__a]) << std::
            endl;
          return 1;
        }
      for (int __a = 0; __a < 16; __a++)
        if (std::abs ((Not[__a] - Not_reference[__a]) / (Not_reference[__a])) >
            1)
        {
          std::cerr << "Mismatch detected in MIC implementation" << std::endl;
          std::cerr << "Variable Not:" << std::endl;
          std::cerr << "seed=" << seed << ", __a=" << __a << std::endl;
          std::cerr << "Not MIC=  " << Not[__a] << std::endl;
          std::cerr << "Not Reference=  " << Not_reference[__a] << std::endl;
          std::cerr << "Not Rel Difference=  " << std::
            abs ((Not[__a] -
                  Not_reference[__a]) / (Not_reference[__a])) << std::endl;
          std::cerr << "Not Abs Difference=  " << std::abs (Not[__a] -
                                                            Not_reference[__a])
            << std::endl;
          return 1;
        }
      for (int __a = 0; __a < 16; __a++)
        if (std::
            abs ((Sqrt[__a] - Sqrt_reference[__a]) / (Sqrt_reference[__a])) > 1)
        {
          std::cerr << "Mismatch detected in MIC implementation" << std::endl;
          std::cerr << "Variable Sqrt:" << std::endl;
          std::cerr << "seed=" << seed << ", __a=" << __a << std::endl;
          std::cerr << "Sqrt MIC=  " << Sqrt[__a] << std::endl;
          std::cerr << "Sqrt Reference=  " << Sqrt_reference[__a] << std::endl;
          std::cerr << "Sqrt Rel Difference=  " << std::
            abs ((Sqrt[__a] -
                  Sqrt_reference[__a]) / (Sqrt_reference[__a])) << std::endl;
          std::cerr << "Sqrt Abs Difference=  " << std::abs (Sqrt[__a] -
                                                             Sqrt_reference
                                                             [__a]) << std::
            endl;
          return 1;
        }
      for (int __a = 0; __a < 16; __a++)
        if (std::
            abs ((RSqrt[__a] - RSqrt_reference[__a]) / (RSqrt_reference[__a])) >
            1)
        {
          std::cerr << "Mismatch detected in MIC implementation" << std::endl;
          std::cerr << "Variable RSqrt:" << std::endl;
          std::cerr << "seed=" << seed << ", __a=" << __a << std::endl;
          std::cerr << "RSqrt MIC=  " << RSqrt[__a] << std::endl;
          std::cerr << "RSqrt Reference=  " << RSqrt_reference[__a] << std::
            endl;
          std::cerr << "RSqrt Rel Difference=  " << std::
            abs ((RSqrt[__a] -
                  RSqrt_reference[__a]) / (RSqrt_reference[__a])) << std::endl;
          std::cerr << "RSqrt Abs Difference=  " << std::abs (RSqrt[__a] -
                                                              RSqrt_reference
                                                              [__a]) << std::
            endl;
          return 1;
        }
      for (int __a = 0; __a < 16; __a++)
        if (std::abs ((Log[__a] - Log_reference[__a]) / (Log_reference[__a])) >
            1)
        {
          std::cerr << "Mismatch detected in MIC implementation" << std::endl;
          std::cerr << "Variable Log:" << std::endl;
          std::cerr << "seed=" << seed << ", __a=" << __a << std::endl;
          std::cerr << "Log MIC=  " << Log[__a] << std::endl;
          std::cerr << "Log Reference=  " << Log_reference[__a] << std::endl;
          std::cerr << "Log Rel Difference=  " << std::
            abs ((Log[__a] -
                  Log_reference[__a]) / (Log_reference[__a])) << std::endl;
          std::cerr << "Log Abs Difference=  " << std::abs (Log[__a] -
                                                            Log_reference[__a])
            << std::endl;
          return 1;
        }
      for (int __a = 0; __a < 16; __a++)
        if (std::abs ((Exp[__a] - Exp_reference[__a]) / (Exp_reference[__a])) >
            1)
        {
          std::cerr << "Mismatch detected in MIC implementation" << std::endl;
          std::cerr << "Variable Exp:" << std::endl;
          std::cerr << "seed=" << seed << ", __a=" << __a << std::endl;
          std::cerr << "Exp MIC=  " << Exp[__a] << std::endl;
          std::cerr << "Exp Reference=  " << Exp_reference[__a] << std::endl;
          std::cerr << "Exp Rel Difference=  " << std::
            abs ((Exp[__a] -
                  Exp_reference[__a]) / (Exp_reference[__a])) << std::endl;
          std::cerr << "Exp Abs Difference=  " << std::abs (Exp[__a] -
                                                            Exp_reference[__a])
            << std::endl;
          return 1;
        }
      for (int __a = 0; __a < 16; __a++)
        if (std::
            abs ((Inverse[__a] -
                  Inverse_reference[__a]) / (Inverse_reference[__a])) > 1)
        {
          std::cerr << "Mismatch detected in MIC implementation" << std::endl;
          std::cerr << "Variable Inverse:" << std::endl;
          std::cerr << "seed=" << seed << ", __a=" << __a << std::endl;
          std::cerr << "Inverse MIC=  " << Inverse[__a] << std::endl;
          std::cerr << "Inverse Reference=  " << Inverse_reference[__a] << std::
            endl;
          std::cerr << "Inverse Rel Difference=  " << std::
            abs ((Inverse[__a] -
                  Inverse_reference[__a]) /
                 (Inverse_reference[__a])) << std::endl;
          std::cerr << "Inverse Abs Difference=  " << std::abs (Inverse[__a] -
                                                                Inverse_reference
                                                                [__a]) << std::
            endl;
          return 1;
        }
      for (int __a = 0; __a < 16; __a++)
        if (std::
            abs ((AbsoluteValue[__a] -
                  AbsoluteValue_reference[__a]) /
                 (AbsoluteValue_reference[__a])) > 1)
        {
          std::cerr << "Mismatch detected in MIC implementation" << std::endl;
          std::cerr << "Variable AbsoluteValue:" << std::endl;
          std::cerr << "seed=" << seed << ", __a=" << __a << std::endl;
          std::cerr << "AbsoluteValue MIC=  " << AbsoluteValue[__a] << std::
            endl;
          std::
            cerr << "AbsoluteValue Reference=  " << AbsoluteValue_reference[__a]
            << std::endl;
          std::cerr << "AbsoluteValue Rel Difference=  " << std::
            abs ((AbsoluteValue[__a] -
                  AbsoluteValue_reference[__a]) /
                 (AbsoluteValue_reference[__a])) << std::endl;
          std::cerr << "AbsoluteValue Abs Difference=  " << std::
            abs (AbsoluteValue[__a] -
                 AbsoluteValue_reference[__a]) << std::endl;
          return 1;
        }
      for (int __a = 0; __a < 16; __a++)
        if (std::
            abs ((Sign[__a] - Sign_reference[__a]) / (Sign_reference[__a])) > 1)
        {
          std::cerr << "Mismatch detected in MIC implementation" << std::endl;
          std::cerr << "Variable Sign:" << std::endl;
          std::cerr << "seed=" << seed << ", __a=" << __a << std::endl;
          std::cerr << "Sign MIC=  " << Sign[__a] << std::endl;
          std::cerr << "Sign Reference=  " << Sign_reference[__a] << std::endl;
          std::cerr << "Sign Rel Difference=  " << std::
            abs ((Sign[__a] -
                  Sign_reference[__a]) / (Sign_reference[__a])) << std::endl;
          std::cerr << "Sign Abs Difference=  " << std::abs (Sign[__a] -
                                                             Sign_reference
                                                             [__a]) << std::
            endl;
          return 1;
        }
      for (int __a = 0; __a < 16; __a++)
        if (std::
            abs ((Minimum[__a] -
                  Minimum_reference[__a]) / (Minimum_reference[__a])) > 1)
        {
          std::cerr << "Mismatch detected in MIC implementation" << std::endl;
          std::cerr << "Variable Minimum:" << std::endl;
          std::cerr << "seed=" << seed << ", __a=" << __a << std::endl;
          std::cerr << "Minimum MIC=  " << Minimum[__a] << std::endl;
          std::cerr << "Minimum Reference=  " << Minimum_reference[__a] << std::
            endl;
          std::cerr << "Minimum Rel Difference=  " << std::
            abs ((Minimum[__a] -
                  Minimum_reference[__a]) /
                 (Minimum_reference[__a])) << std::endl;
          std::cerr << "Minimum Abs Difference=  " << std::abs (Minimum[__a] -
                                                                Minimum_reference
                                                                [__a]) << std::
            endl;
          return 1;
        }
      for (int __a = 0; __a < 16; __a++)
        if (std::
            abs ((Maximum[__a] -
                  Maximum_reference[__a]) / (Maximum_reference[__a])) > 1)
        {
          std::cerr << "Mismatch detected in MIC implementation" << std::endl;
          std::cerr << "Variable Maximum:" << std::endl;
          std::cerr << "seed=" << seed << ", __a=" << __a << std::endl;
          std::cerr << "Maximum MIC=  " << Maximum[__a] << std::endl;
          std::cerr << "Maximum Reference=  " << Maximum_reference[__a] << std::
            endl;
          std::cerr << "Maximum Rel Difference=  " << std::
            abs ((Maximum[__a] -
                  Maximum_reference[__a]) /
                 (Maximum_reference[__a])) << std::endl;
          std::cerr << "Maximum Abs Difference=  " << std::abs (Maximum[__a] -
                                                                Maximum_reference
                                                                [__a]) << std::
            endl;
          return 1;
        }
      for (int __a = 0; __a < 16; __a++)
        if (std::
            abs ((Blend[__a] - Blend_reference[__a]) / (Blend_reference[__a])) >
            1)
        {
          std::cerr << "Mismatch detected in MIC implementation" << std::endl;
          std::cerr << "Variable Blend:" << std::endl;
          std::cerr << "seed=" << seed << ", __a=" << __a << std::endl;
          std::cerr << "Blend MIC=  " << Blend[__a] << std::endl;
          std::cerr << "Blend Reference=  " << Blend_reference[__a] << std::
            endl;
          std::cerr << "Blend Rel Difference=  " << std::
            abs ((Blend[__a] -
                  Blend_reference[__a]) / (Blend_reference[__a])) << std::endl;
          std::cerr << "Blend Abs Difference=  " << std::abs (Blend[__a] -
                                                              Blend_reference
                                                              [__a]) << std::
            endl;
          return 1;
        }

    }
#endif

  }



  std::cout << "SIMD check successful!" << std::endl;

  return 0;

}

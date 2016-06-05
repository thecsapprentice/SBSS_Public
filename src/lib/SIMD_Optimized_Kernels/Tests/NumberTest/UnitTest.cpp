
#include <cstdlib>
#include <iostream>

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
  typedef float T;

  int seed = 1;
  if (argc == 2)
    seed = atoi (argv[1]);
  srand (seed);



  {
    T A __attribute__ ((aligned (4)));
    T B __attribute__ ((aligned (4)));
    T C __attribute__ ((aligned (4)));
    T D __attribute__ ((aligned (4)));
    T Add __attribute__ ((aligned (4)));
    T Add_reference __attribute__ ((aligned (4)));
    T Add_original __attribute__ ((aligned (4)));
    T Multiply __attribute__ ((aligned (4)));
    T Multiply_reference __attribute__ ((aligned (4)));
    T Multiply_original __attribute__ ((aligned (4)));
    T Subtract __attribute__ ((aligned (4)));
    T Subtract_reference __attribute__ ((aligned (4)));
    T Subtract_original __attribute__ ((aligned (4)));
    T Divide __attribute__ ((aligned (4)));
    T Divide_reference __attribute__ ((aligned (4)));
    T Divide_original __attribute__ ((aligned (4)));
    T LessThan __attribute__ ((aligned (4)));
    T LessThan_reference __attribute__ ((aligned (4)));
    T LessThan_original __attribute__ ((aligned (4)));
    T GreaterThan __attribute__ ((aligned (4)));
    T GreaterThan_reference __attribute__ ((aligned (4)));
    T GreaterThan_original __attribute__ ((aligned (4)));
    T LessEquals __attribute__ ((aligned (4)));
    T LessEquals_reference __attribute__ ((aligned (4)));
    T LessEquals_original __attribute__ ((aligned (4)));
    T GreaterEquals __attribute__ ((aligned (4)));
    T GreaterEquals_reference __attribute__ ((aligned (4)));
    T GreaterEquals_original __attribute__ ((aligned (4)));
    T Equals __attribute__ ((aligned (4)));
    T Equals_reference __attribute__ ((aligned (4)));
    T Equals_original __attribute__ ((aligned (4)));
    T And __attribute__ ((aligned (4)));
    T And_reference __attribute__ ((aligned (4)));
    T And_original __attribute__ ((aligned (4)));
    T Or __attribute__ ((aligned (4)));
    T Or_reference __attribute__ ((aligned (4)));
    T Or_original __attribute__ ((aligned (4)));
    T Xor __attribute__ ((aligned (4)));
    T Xor_reference __attribute__ ((aligned (4)));
    T Xor_original __attribute__ ((aligned (4)));
    T AndNot __attribute__ ((aligned (4)));
    T AndNot_reference __attribute__ ((aligned (4)));
    T AndNot_original __attribute__ ((aligned (4)));
    T Not __attribute__ ((aligned (4)));
    T Not_reference __attribute__ ((aligned (4)));
    T Not_original __attribute__ ((aligned (4)));
    T Sqrt __attribute__ ((aligned (4)));
    T Sqrt_reference __attribute__ ((aligned (4)));
    T Sqrt_original __attribute__ ((aligned (4)));
    T RSqrt __attribute__ ((aligned (4)));
    T RSqrt_reference __attribute__ ((aligned (4)));
    T RSqrt_original __attribute__ ((aligned (4)));
    T Log __attribute__ ((aligned (4)));
    T Log_reference __attribute__ ((aligned (4)));
    T Log_original __attribute__ ((aligned (4)));
    T Exp __attribute__ ((aligned (4)));
    T Exp_reference __attribute__ ((aligned (4)));
    T Exp_original __attribute__ ((aligned (4)));
    T Inverse __attribute__ ((aligned (4)));
    T Inverse_reference __attribute__ ((aligned (4)));
    T Inverse_original __attribute__ ((aligned (4)));
    T AbsoluteValue __attribute__ ((aligned (4)));
    T AbsoluteValue_reference __attribute__ ((aligned (4)));
    T AbsoluteValue_original __attribute__ ((aligned (4)));
    T Sign __attribute__ ((aligned (4)));
    T Sign_reference __attribute__ ((aligned (4)));
    T Sign_original __attribute__ ((aligned (4)));
    T Minimum __attribute__ ((aligned (4)));
    T Minimum_reference __attribute__ ((aligned (4)));
    T Minimum_original __attribute__ ((aligned (4)));
    T Maximum __attribute__ ((aligned (4)));
    T Maximum_reference __attribute__ ((aligned (4)));
    T Maximum_original __attribute__ ((aligned (4)));
    T Blend __attribute__ ((aligned (4)));
    T Blend_reference __attribute__ ((aligned (4)));
    T Blend_original __attribute__ ((aligned (4)));



    A = Get_Random < float >();
    B = Get_Random < float >();
    C = Get_Random < float >();
    D = Get_Random < float >();
    {
      Add_original = Get_Random < float >();
      Add = Add_original;
      Add_reference = Add_original;
    }
    {
      Multiply_original = Get_Random < float >();
      Multiply = Multiply_original;
      Multiply_reference = Multiply_original;
    }
    {
      Subtract_original = Get_Random < float >();
      Subtract = Subtract_original;
      Subtract_reference = Subtract_original;
    }
    {
      Divide_original = Get_Random < float >();
      Divide = Divide_original;
      Divide_reference = Divide_original;
    }
    {
      LessThan_original = Get_Random < float >();
      LessThan = LessThan_original;
      LessThan_reference = LessThan_original;
    }
    {
      GreaterThan_original = Get_Random < float >();
      GreaterThan = GreaterThan_original;
      GreaterThan_reference = GreaterThan_original;
    }
    {
      LessEquals_original = Get_Random < float >();
      LessEquals = LessEquals_original;
      LessEquals_reference = LessEquals_original;
    }
    {
      GreaterEquals_original = Get_Random < float >();
      GreaterEquals = GreaterEquals_original;
      GreaterEquals_reference = GreaterEquals_original;
    }
    {
      Equals_original = Get_Random < float >();
      Equals = Equals_original;
      Equals_reference = Equals_original;
    }
    {
      And_original = Get_Random < float >();
      And = And_original;
      And_reference = And_original;
    }
    {
      Or_original = Get_Random < float >();
      Or = Or_original;
      Or_reference = Or_original;
    }
    {
      Xor_original = Get_Random < float >();
      Xor = Xor_original;
      Xor_reference = Xor_original;
    }
    {
      AndNot_original = Get_Random < float >();
      AndNot = AndNot_original;
      AndNot_reference = AndNot_original;
    }
    {
      Not_original = Get_Random < float >();
      Not = Not_original;
      Not_reference = Not_original;
    }
    {
      Sqrt_original = Get_Random < float >();
      Sqrt = Sqrt_original;
      Sqrt_reference = Sqrt_original;
    }
    {
      RSqrt_original = Get_Random < float >();
      RSqrt = RSqrt_original;
      RSqrt_reference = RSqrt_original;
    }
    {
      Log_original = Get_Random < float >();
      Log = Log_original;
      Log_reference = Log_original;
    }
    {
      Exp_original = Get_Random < float >();
      Exp = Exp_original;
      Exp_reference = Exp_original;
    }
    {
      Inverse_original = Get_Random < float >();
      Inverse = Inverse_original;
      Inverse_reference = Inverse_original;
    }
    {
      AbsoluteValue_original = Get_Random < float >();
      AbsoluteValue = AbsoluteValue_original;
      AbsoluteValue_reference = AbsoluteValue_original;
    }
    {
      Sign_original = Get_Random < float >();
      Sign = Sign_original;
      Sign_reference = Sign_original;
    }
    {
      Minimum_original = Get_Random < float >();
      Minimum = Minimum_original;
      Minimum_reference = Minimum_original;
    }
    {
      Maximum_original = Get_Random < float >();
      Maximum = Maximum_original;
      Maximum_reference = Maximum_original;
    }
    {
      Blend_original = Get_Random < float >();
      Blend = Blend_original;
      Blend_reference = Blend_original;
    }


    Add = Add_original;
    Multiply = Multiply_original;
    Subtract = Subtract_original;
    Divide = Divide_original;
    LessThan = LessThan_original;
    GreaterThan = GreaterThan_original;
    LessEquals = LessEquals_original;
    GreaterEquals = GreaterEquals_original;
    Equals = Equals_original;
    And = And_original;
    Or = Or_original;
    Xor = Xor_original;
    AndNot = AndNot_original;
    Not = Not_original;
    Sqrt = Sqrt_original;
    RSqrt = RSqrt_original;
    Log = Log_original;
    Exp = Exp_original;
    Inverse = Inverse_original;
    AbsoluteValue = AbsoluteValue_original;
    Sign = Sign_original;
    Minimum = Minimum_original;
    Maximum = Maximum_original;
    Blend = Blend_original;
    for (int i = 0; i < 1; i += 1)
      {
        NumberTest < float, float, int >(A, B, C, D, Add, Multiply, Subtract,
                                         Divide, LessThan, GreaterThan,
                                         LessEquals, GreaterEquals, Equals, And,
                                         Or, Xor, AndNot, Not, Sqrt, RSqrt, Log,
                                         Exp, Inverse, AbsoluteValue, Sign,
                                         Minimum, Maximum, Blend);
      }

    NumberTest_Reference < float >(A, B, C, D, Add_reference,
                                   Multiply_reference, Subtract_reference,
                                   Divide_reference, LessThan_reference,
                                   GreaterThan_reference, LessEquals_reference,
                                   GreaterEquals_reference, Equals_reference,
                                   And_reference, Or_reference, Xor_reference,
                                   AndNot_reference, Not_reference,
                                   Sqrt_reference, RSqrt_reference,
                                   Log_reference, Exp_reference,
                                   Inverse_reference, AbsoluteValue_reference,
                                   Sign_reference, Minimum_reference,
                                   Maximum_reference, Blend_reference);
    if (!
        (NumberTest_Compare <
         float >(Add, Multiply, Subtract, Divide, LessThan, GreaterThan,
                 LessEquals, GreaterEquals, Equals, And, Or, Xor, AndNot, Not,
                 Sqrt, RSqrt, Log, Exp, Inverse, AbsoluteValue, Sign, Minimum,
                 Maximum, Blend, Add_reference, Multiply_reference,
                 Subtract_reference, Divide_reference, LessThan_reference,
                 GreaterThan_reference, LessEquals_reference,
                 GreaterEquals_reference, Equals_reference, And_reference,
                 Or_reference, Xor_reference, AndNot_reference, Not_reference,
                 Sqrt_reference, RSqrt_reference, Log_reference, Exp_reference,
                 Inverse_reference, AbsoluteValue_reference, Sign_reference,
                 Minimum_reference, Maximum_reference, Blend_reference)))
      {
        std::cout << "Failed to confirm unit test for NumberTest " << std::endl;
        return 1;
      }

  }



  return 0;

}

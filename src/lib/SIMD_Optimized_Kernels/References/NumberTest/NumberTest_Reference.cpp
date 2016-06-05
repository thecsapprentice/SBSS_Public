//#####################################################################
//  Copyright (c) 2011-2013 Nathan Mitchell, Eftychios Sifakis.
//  This file is covered by the FreeBSD license. Please refer to the 
//  license.txt file for more information.
//#####################################################################



#include <cmath>
#include <iostream>

namespace {
    typedef union {
        int i;
        float f;
    } floatConverter;
}

template<class T>
void NumberTest_Reference(const T &A,
                          const T &B,
                          const T &C,
                          const T &D,
                          
                          T &Add,
                          T &Multiply,
                          T &Subtract,
                          T &Divide,
                          
                          T &LessThan,
                          T &GreaterThan,
                          T &LessEquals,
                          T &GreaterEquals,
                          T &Equals,
                          
                          T &And,
                          T &Or,
                          T &Xor,
                          T &AndNot,
                          T &Not,
                          
                          T &Sqrt,
                          T &RSqrt,
                          T &Log,
                          T &Exp,
                          T &Inverse,
                          T &AbsoluteValue,
                          T &Sign,
                          
                          T &Minimum,
                          T &Maximum,
                          T &Blend )
{

    floatConverter ALL_ONES;
    ALL_ONES.i = (0xFFFFFFFF);
    floatConverter iA, iB, iC, iD, temp;
    iA.f = A;     iB.f = B; 
    iC.f = C;     iD.f = D; 

    Add = A + B;
    Multiply = A * B;
    Subtract = A - B;
    Divide = A / B;
    
    LessThan = (A < B) ? ALL_ONES.f : 0;
    GreaterThan = (A > B) ? ALL_ONES.f : 0;
    LessEquals = (A <= B) ? ALL_ONES.f : 0;
    GreaterEquals = (A >= B) ? ALL_ONES.f : 0;
    Equals = (A == B) ? ALL_ONES.f : 0;
       
    temp.i = (iA.i & iB.i);
    And = temp.f;

    temp.i = (iA.i | iB.i);
    Or = temp.f;

    temp.i = (iA.i ^ iB.i);
    Xor = temp.f;

    temp.i = (~iA.i & (iB.i));
    AndNot = temp.f;

    temp.i = (~iA.i);
    Not = temp.f;

    Sqrt = sqrt( A );
    RSqrt = 1.f / sqrt(A);
    Log = log(A);
    Exp = exp(A);
    Inverse = 1.f / A;
    AbsoluteValue = fabs(A);
    Sign = (A >= 0) ? 1.0f : -1.0f;

    Minimum = (A < B) ? A : B;
    Maximum = (A > B) ? A : B;

    Blend = (C < D) ? B : A;
}

template<class T>
bool NumberTest_Compare(const T &Add,
                        const T &Multiply,
                        const T &Subtract,
                        const T &Divide,
                        
                        const T &LessThan,
                        const T &GreaterThan,
                        const T &LessEquals,
                        const T &GreaterEquals,
                        const T &Equals,
                        
                        const T &And,
                        const T &Or,
                        const T &Xor,
                        const T &AndNot,
                        const T &Not,
                        
                        const T &Sqrt,
                        const T &RSqrt,
                        const T &Log,
                        const T &Exp,
                        const T &Inverse,
                        const T &AbsoluteValue,
                        const T &Sign,
                        
                        const T &Minimum,
                        const T &Maximum,
                        const T &Blend, 

                        const T &Add_reference,
                        const T &Multiply_reference,
                        const T &Subtract_reference,
                        const T &Divide_reference,
                        
                        const T &LessThan_reference,
                        const T &GreaterThan_reference,
                        const T &LessEquals_reference,
                        const T &GreaterEquals_reference,
                        const T &Equals_reference,
                        
                        const T &And_reference,
                        const T &Or_reference,
                        const T &Xor_reference,
                        const T &AndNot_reference,
                        const T &Not_reference,
                        
                        const T &Sqrt_reference,
                        const T &RSqrt_reference,
                        const T &Log_reference,
                        const T &Exp_reference,
                        const T &Inverse_reference,
                        const T &AbsoluteValue_reference,
                        const T &Sign_reference,
                        
                        const T &Minimum_reference,
                        const T &Maximum_reference,
                        const T &Blend_reference
                        )
{


    std::cout << "Add Compare   :  " << Add << std::endl;
    std::cout << "Add Reference :  " << Add_reference << std::endl;
    std::cout << "Difference    :  " << Add - Add_reference << std::endl << std::endl;

    std::cout << "Multiply Compare   :  " << Multiply << std::endl;
    std::cout << "Multiply Reference :  " << Multiply_reference << std::endl;
    std::cout << "Difference    :  " << Multiply - Multiply_reference << std::endl << std::endl;

    std::cout << "Subtract Compare   :  " << Subtract << std::endl;
    std::cout << "Subtract Reference :  " << Subtract_reference << std::endl;
    std::cout << "Difference    :  " << Subtract - Subtract_reference << std::endl << std::endl;

    std::cout << "Divide Compare   :  " << Divide << std::endl;
    std::cout << "Divide Reference :  " << Divide_reference << std::endl;
    std::cout << "Difference    :  " << Divide - Divide_reference << std::endl << std::endl;

    std::cout << "LessThan Compare   :  " << LessThan << std::endl;
    std::cout << "LessThan Reference :  " << LessThan_reference << std::endl;
    std::cout << "Difference    :  " << LessThan - LessThan_reference << std::endl << std::endl;

    std::cout << "GreaterThan Compare   :  " << GreaterThan << std::endl;
    std::cout << "GreaterThan Reference :  " << GreaterThan_reference << std::endl;
    std::cout << "Difference    :  " << GreaterThan - GreaterThan_reference << std::endl << std::endl;

    std::cout << "LessEquals Compare   :  " << LessEquals << std::endl;
    std::cout << "LessEquals Reference :  " << LessEquals_reference << std::endl;
    std::cout << "Difference    :  " << LessEquals - LessEquals_reference << std::endl << std::endl;

    std::cout << "GreaterEquals Compare   :  " << GreaterEquals << std::endl;
    std::cout << "GreaterEquals Reference :  " << GreaterEquals_reference << std::endl;
    std::cout << "Difference    :  " << GreaterEquals - GreaterEquals_reference << std::endl << std::endl;

    std::cout << "Equals Compare   :  " << Equals << std::endl;
    std::cout << "Equals Reference :  " << Equals_reference << std::endl;
    std::cout << "Difference    :  " << Equals - Equals_reference << std::endl << std::endl;

    std::cout << "And Compare   :  " << And << std::endl;
    std::cout << "And Reference :  " << And_reference << std::endl;
    std::cout << "Difference    :  " << And - And_reference << std::endl << std::endl;

    std::cout << "Or Compare   :  " << Or << std::endl;
    std::cout << "Or Reference :  " << Or_reference << std::endl;
    std::cout << "Difference    :  " << Or - Or_reference << std::endl << std::endl;

    std::cout << "Xor Compare   :  " << Xor << std::endl;
    std::cout << "Xor Reference :  " << Xor_reference << std::endl;
    std::cout << "Difference    :  " << Xor - Xor_reference << std::endl << std::endl;

    std::cout << "AndNot Compare   :  " << AndNot << std::endl;
    std::cout << "AndNot Reference :  " << AndNot_reference << std::endl;
    std::cout << "Difference    :  " << AndNot - AndNot_reference << std::endl << std::endl;

    std::cout << "Not Compare   :  " << Not << std::endl;
    std::cout << "Not Reference :  " << Not_reference << std::endl;
    std::cout << "Difference    :  " << Not - Not_reference << std::endl << std::endl;

    std::cout << "Sqrt Compare   :  " << Sqrt << std::endl;
    std::cout << "Sqrt Reference :  " << Sqrt_reference << std::endl;
    std::cout << "Difference    :  " << Sqrt - Sqrt_reference << std::endl << std::endl;

    std::cout << "RSqrt Compare   :  " << RSqrt << std::endl;
    std::cout << "RSqrt Reference :  " << RSqrt_reference << std::endl;
    std::cout << "Difference    :  " << RSqrt - RSqrt_reference << std::endl << std::endl;

    std::cout << "Log Compare   :  " << Log << std::endl;
    std::cout << "Log Reference :  " << Log_reference << std::endl;
    std::cout << "Difference    :  " << Log - Log_reference << std::endl << std::endl;

    std::cout << "Exp Compare   :  " << Exp << std::endl;
    std::cout << "Exp Reference :  " << Exp_reference << std::endl;
    std::cout << "Difference    :  " << Exp - Exp_reference << std::endl << std::endl;

    std::cout << "Inverse Compare   :  " << Inverse << std::endl;
    std::cout << "Inverse Reference :  " << Inverse_reference << std::endl;
    std::cout << "Difference    :  " << Inverse - Inverse_reference << std::endl << std::endl;

    std::cout << "AbsoluteValue Compare   :  " << AbsoluteValue << std::endl;
    std::cout << "AbsoluteValue Reference :  " << AbsoluteValue_reference << std::endl;
    std::cout << "Difference    :  " << AbsoluteValue - AbsoluteValue_reference << std::endl << std::endl;

    std::cout << "Sign Compare   :  " << Sign << std::endl;
    std::cout << "Sign Reference :  " << Sign_reference << std::endl;
    std::cout << "Difference    :  " << Sign - Sign_reference << std::endl << std::endl;

    std::cout << "Minimum Compare   :  " << Minimum << std::endl;
    std::cout << "Minimum Reference :  " << Minimum_reference << std::endl;
    std::cout << "Difference    :  " << Minimum - Minimum_reference << std::endl << std::endl;

    std::cout << "Maximum Compare   :  " << Maximum << std::endl;
    std::cout << "Maximum Reference :  " << Maximum_reference << std::endl;
    std::cout << "Difference    :  " << Maximum - Maximum_reference << std::endl << std::endl;

    std::cout << "Blend Compare   :  " << Blend << std::endl;
    std::cout << "Blend Reference :  " << Blend_reference << std::endl;
    std::cout << "Difference    :  " << Blend - Blend_reference << std::endl << std::endl;

    return true;
}



template void NumberTest_Reference(const float &A,
                                   const float &B,
                                   const float &C,
                                   const float &D,
                                   
                                   float &Add,
                                   float &Multiply,
                                   float &Subtract,
                                   float &Divide,
                                   
                                   float &LessThan,
                                   float &GreaterThan,
                                   float &LessEquals,
                                   float &GreaterEquals,
                                   float &Equals,
                                   
                                   float &And,
                                   float &Or,
                                   float &Xor,
                                   float &AndNot,
                                   float &Not,
                                   
                                   float &Sqrt,
                                   float &RSqrt,
                                   float &Log,
                                   float &Exp,
                                   float &Inverse,
                                   float &AbsoluteValue,
                                   float &Sign,
                                   
                                   float &Minimum,
                                   float &Maximum,
                                   float &Blend );

template bool NumberTest_Compare(const float &Add,
                                 const float &Multiply,
                                 const float &Subtract,
                                 const float &Divide,
                                 
                                 const float &LessThan,
                                 const float &GreaterThan,
                                 const float &LessEquals,
                                 const float &GreaterEquals,
                                 const float &Equals,
                                 
                                 const float &And,
                                 const float &Or,
                                 const float &Xor,
                                 const float &AndNot,
                                 const float &Not,
                                 
                                 const float &Sqrt,
                                 const float &RSqrt,
                                 const float &Log,
                                 const float &Exp,
                                 const float &Inverse,
                                 const float &AbsoluteValue,
                                 const float &Sign,
                                 
                                 const float &Minimum,
                                 const float &Maximum,
                                 const float &Blend, 
                                 
                                 const float &Add_reference,
                                 const float &Multiply_reference,
                                 const float &Subtract_reference,
                                 const float &Divide_reference,
                                 
                                 const float &LessThan_reference,
                                 const float &GreaterThan_reference,
                                 const float &LessEquals_reference,
                                 const float &GreaterEquals_reference,
                                 const float &Equals_reference,
                                 
                                 const float &And_reference,
                                 const float &Or_reference,
                                 const float &Xor_reference,
                                 const float &AndNot_reference,
                                 const float &Not_reference,
                                 
                                 const float &Sqrt_reference,
                                 const float &RSqrt_reference,
                                 const float &Log_reference,
                                 const float &Exp_reference,
                                 const float &Inverse_reference,
                                 const float &AbsoluteValue_reference,
                                 const float &Sign_reference,
                                 
                                 const float &Minimum_reference,
                                 const float &Maximum_reference,
                                 const float &Blend_reference
                                 );

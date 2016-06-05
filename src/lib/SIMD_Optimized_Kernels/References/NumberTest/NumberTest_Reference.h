//#####################################################################
//  Copyright (c) 2011-2013 Nathan Mitchell, Eftychios Sifakis.
//  This file is covered by the FreeBSD license. Please refer to the 
//  license.txt file for more information.
//#####################################################################


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
                          T &Blend );

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
                        );

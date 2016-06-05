//#####################################################################
// Copyright 2010, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __Vector_SAXPY_V9_Helper__
#define __Vector_SAXPY_V9_Helper__
namespace MT_Streaming_Kernels{
    
template<class T, bool ACCUM>
    class Vector_SAXPY_V9_Helper
{
    T* const x_f;

    const T* const v_f_1;
    const T* const v_f_2;
    const T* const v_f_3;
    const T* const v_f_4;
    const T* const v_f_5;
    const T* const v_f_6;
    const T* const v_f_7;
    const T* const v_f_8;
    const T* const v_f_9;

    const T* const y_f_1;
    const T* const y_f_2;
    const T* const y_f_3;
    const T* const y_f_4;
    const T* const y_f_5;
    const T* const y_f_6;
    const T* const y_f_7;
    const T* const y_f_8;
    const T* const y_f_9;


    const int size;

public:
    explicit Vector_SAXPY_V9_Helper(T* const x_f_input,
                                    const T* const v_f_1_input,
                                    const T* const v_f_2_input,
                                    const T* const v_f_3_input,
                                    const T* const v_f_4_input,
                                    const T* const v_f_5_input,
                                    const T* const v_f_6_input,
                                    const T* const v_f_7_input,
                                    const T* const v_f_8_input,
                                    const T* const v_f_9_input,
                                   
                                    const T* const y_f_1_input,
                                    const T* const y_f_2_input,
                                    const T* const y_f_3_input,
                                    const T* const y_f_4_input,
                                    const T* const y_f_5_input,
                                    const T* const y_f_6_input,
                                    const T* const y_f_7_input,
                                    const T* const y_f_8_input,
                                    const T* const y_f_9_input,
                               

                                    const int size_input)
        :x_f(x_f_input),
        v_f_1(v_f_1_input),
        v_f_2(v_f_2_input),
        v_f_3(v_f_3_input),
        v_f_4(v_f_4_input),
        v_f_5(v_f_5_input),
        v_f_6(v_f_6_input),
        v_f_7(v_f_7_input),
        v_f_8(v_f_8_input),
        v_f_9(v_f_9_input),

        y_f_1(y_f_1_input),
        y_f_2(y_f_2_input),
        y_f_3(y_f_3_input),
        y_f_4(y_f_4_input),
        y_f_5(y_f_5_input),
        y_f_6(y_f_6_input),
        y_f_7(y_f_7_input),
        y_f_8(y_f_8_input),
        y_f_9(y_f_9_input),

        size(size_input)
    {}

    void Run()
    {
        Run_Index_Range(0,size-1);
    }
  
    // For debugging purposes only

    static void Allocate_Data(T*& x_f,
                              T*& v_f_1,
                              T*& v_f_2,
                              T*& v_f_3,
                              T*& v_f_4,
                              T*& v_f_5,
                              T*& v_f_6,
                              T*& v_f_7,
                              T*& v_f_8,
                              T*& v_f_9,

                              T*& y_f_1,
                              T*& y_f_2,
                              T*& y_f_3,
                              T*& y_f_4,
                              T*& y_f_5,
                              T*& y_f_6,
                              T*& y_f_7,
                              T*& y_f_8,
                              T*& y_f_9,

                              const int size)
    {
        x_f=new T[size];

        v_f_1=new T[size];
        v_f_2=new T[size];
        v_f_3=new T[size];
        v_f_4=new T[size];
        v_f_5=new T[size];
        v_f_6=new T[size];
        v_f_7=new T[size];
        v_f_8=new T[size];
        v_f_9=new T[size];

        y_f_1=new T[size];
        y_f_2=new T[size];
        y_f_3=new T[size];
        y_f_4=new T[size];
        y_f_5=new T[size];
        y_f_6=new T[size];
        y_f_7=new T[size];
        y_f_8=new T[size];
        y_f_9=new T[size];

    }

    static void Initialize_Data(T* x_f,
                                T* v_f_1,
                                T* v_f_2,
                                T* v_f_3,
                                T* v_f_4,
                                T* v_f_5,
                                T* v_f_6,
                                T* v_f_7,
                                T* v_f_8,
                                T* v_f_9,

                                T* y_f_1,
                                T* y_f_2,
                                T* y_f_3,
                                T* y_f_4,
                                T* y_f_5,
                                T* y_f_6,
                                T* y_f_7,
                                T* y_f_8,
                                T* y_f_9,

                              const int size)
    {
        for(int i=0;i<size;i++){
            x_f[i]=T(i);
            v_f_1[i]=v_f_2[i]=v_f_3[i]=T(i); 
            v_f_4[i]=v_f_5[i]=v_f_6[i]=T(i); 
            v_f_7[i]=v_f_8[i]=v_f_9[i]=T(i); 

            y_f_1[i]=y_f_2[i]=y_f_3[i]=T(i); 
            y_f_4[i]=y_f_5[i]=y_f_6[i]=T(i); 
            y_f_7[i]=y_f_8[i]=y_f_9[i]=T(i); 

        }
    }

//#####################################################################
    void Run_Parallel(const int number_of_partitions);
    void Run_Index_Range(const int index_start,const int index_end);
//#####################################################################
};
}
#endif

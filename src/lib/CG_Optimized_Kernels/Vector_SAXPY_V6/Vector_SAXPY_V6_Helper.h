//#####################################################################
// Copyright 2010, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __Vector_SAXPY_V6_Helper__
#define __Vector_SAXPY_V6_Helper__
namespace MT_Streaming_Kernels{
    
template<class T, bool ACCUM>
    class Vector_SAXPY_V6_Helper
{
    T* const x_f;
    T* const x_g;

    const T* const v_f_1;
    const T* const v_f_2;
    const T* const v_f_3;

    const T* const v_g_1;
    const T* const v_g_2;
    const T* const v_g_3;

    const T* const y_f_1;
    const T* const y_f_2;
    const T* const y_f_3;

    const T* const y_g_1;
    const T* const y_g_2;
    const T* const y_g_3;

    const int size;

public:
    explicit Vector_SAXPY_V6_Helper(T* const x_f_input, T* const x_g_input, 
                                    const T* const v_f_1_input,
                                    const T* const v_f_2_input,
                                    const T* const v_f_3_input,

                                    const T* const v_g_1_input,
                                    const T* const v_g_2_input,
                                    const T* const v_g_3_input,
                                    
                                    const T* const y_f_1_input,
                                    const T* const y_f_2_input,
                                    const T* const y_f_3_input,

                                    const T* const y_g_1_input,
                                    const T* const y_g_2_input,
                                    const T* const y_g_3_input,                                  

                                    const int size_input)
        :x_f(x_f_input),x_g(x_g_input),
        v_f_1(v_f_1_input),
        v_f_2(v_f_2_input),
        v_f_3(v_f_3_input),

        v_g_1(v_g_1_input),
        v_g_2(v_g_2_input),
        v_g_3(v_g_3_input),

        y_f_1(y_f_1_input),
        y_f_2(y_f_2_input),
        y_f_3(y_f_3_input),

        y_g_1(y_g_1_input),
        y_g_2(y_g_2_input),
        y_g_3(y_g_3_input),

        size(size_input)
    {}

    void Run()
    {
        Run_Index_Range(0,size-1);
    }
  
    // For debugging purposes only

    static void Allocate_Data(T*& x_f,T*& x_g,
                              T*& v_f_1,
                              T*& v_f_2,
                              T*& v_f_3,
                              T*& v_g_1,
                              T*& v_g_2,
                              T*& v_g_3,
                              T*& y_f_1,
                              T*& y_f_2,
                              T*& y_f_3,
                              T*& y_g_1,
                              T*& y_g_2,
                              T*& y_g_3,
                              const int size)
    {
        x_f=new T[size];
        x_g=new T[size];

        v_f_1=new T[size];
        v_f_2=new T[size];
        v_f_3=new T[size];

        v_g_1=new T[size];
        v_g_2=new T[size];
        v_g_3=new T[size];

        y_f_1=new T[size];
        y_f_2=new T[size];
        y_f_3=new T[size];

        y_g_1=new T[size];
        y_g_2=new T[size];
        y_g_3=new T[size];
    }

    static void Initialize_Data(T* x_f,T* x_g,
                              T* v_f_1,
                              T* v_f_2,
                              T* v_f_3,
                              T* v_g_1,
                              T* v_g_2,
                              T* v_g_3,
                              T* y_f_1,
                              T* y_f_2,
                              T* y_f_3,
                              T* y_g_1,
                              T* y_g_2,
                              T* y_g_3,
                              const int size)
    {
        for(int i=0;i<size;i++){
            x_f[i]=x_g[i]=T(i);
            v_f_1[i]=v_f_2[i]=v_f_3[i]=T(i); 
            v_g_1[i]=v_g_2[i]=v_g_3[i]=T(i); 
            y_f_1[i]=y_f_2[i]=y_f_3[i]=T(i); 
            y_g_1[i]=y_g_2[i]=y_g_3[i]=T(i); 
        }
    }

//#####################################################################
    void Run_Parallel(const int number_of_partitions);
    void Run_Index_Range(const int index_start,const int index_end);
//#####################################################################
};
}
#endif

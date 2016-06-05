//#####################################################################
// Copyright 2010, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __Vector_SAXPY_V_Helper__
#define __Vector_SAXPY_V_Helper__
namespace MT_Streaming_Kernels{

    template<class T, bool ACCUM>
class Vector_SAXPY_V_Helper
{
    T* const x;
    const T* const v;
    const T* const y;
    const int size;

public:
    explicit Vector_SAXPY_V_Helper(T* const x_input,const T* const y_input,const T* const v_input, const int size_input)
        :x(x_input),v(v_input),y(y_input),size(size_input)
    {}

    void Run()
    {
        Run_Index_Range(0,size-1);
    }
  
    // For debugging purposes only

    static void Allocate_Data(T*& x,T*& v,T*& y,const int size)
    {
        x=new T[size];
        v=new T[size];
        y=new T[size];
    }

    static void Initialize_Data(T* const x,T* const v,T* const y,const int size)
    {
        for(int i=0;i<size;i++)
            x[i]=v[i]=y[i]=(T)i;
    }

//#####################################################################
    void Run_Parallel(const int number_of_partitions);
    void Run_Index_Range(const int index_start,const int index_end);
//#####################################################################
};
}
#endif

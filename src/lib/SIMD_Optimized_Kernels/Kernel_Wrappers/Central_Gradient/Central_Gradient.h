//#####################################################################
//  Copyright (c) 2011-2013 Nathan Mitchell, Eftychios Sifakis.
//  This file is covered by the FreeBSD license. Please refer to the 
//  license.txt file for more information.
//#####################################################################

//#####################################################################
// Function Central_Gradient
//#####################################################################
#ifndef __Central_Gradient_Wrapper__
#define __Central_Gradient_Wrapper__

#include "../../Kernels/Unweighted_Gradient/Unweighted_Gradient.h"

namespace PhysBAM{

namespace {

    template<class T,int d> void
        Gradient_Matrix(MATRIX_MXN<T>& G,const VECTOR<T,d>& weights,const T h)
    {
        typedef VECTOR<int,d> T_INDEX;
        enum {vertices_per_cell=POWER<2,d>::value};
        
        G.Resize(d,vertices_per_cell);
        PHYSBAM_ASSERT(weights.Min()>=0. && weights.Max()<=1.);
        
        int vertex=1;
        for(RANGE_ITERATOR<d> iterator(RANGE<T_INDEX>(T_INDEX(),T_INDEX::All_Ones_Vector()));iterator.Valid();iterator.Next(),vertex++){
            const T_INDEX& index=iterator.Index();
            
            for(int v=1;v<=d;v++){
                T& coefficient=G(v,vertex);
                coefficient=1.;
                
                for(int w=1;w<=d;w++)
                    if(w==v)
                        if(index(w)==0) coefficient*=-1./h; else coefficient*=1./h;
                    else
                        if(index(w)==0) coefficient*=(1.-weights(w)); else coefficient*=weights(w);}}
    }

template<class T,int d> inline void
Central_Gradient(const MATRIX_MXN<T>& Du, MATRIX<T,d>& F,const T one_over_h)
{
    MATRIX_MXN<T> G;
    Gradient_Matrix(G,VECTOR<T,d>::All_Ones_Vector()*(T).5,1.f/one_over_h);
    F=MATRIX<T,d>(Du*G.Transposed());
}

template<> inline void
Central_Gradient(const MATRIX_MXN<float>& Du, MATRIX<float,3>& F,const float one_over_h)
{
    typedef const float (&refConstMat3x8)[3][8];
    typedef float (&refMat3)[9];

    ::Unweighted_Gradient<float,float,int>(
        reinterpret_cast<refConstMat3x8>(Du.Transposed().x[0]),
        reinterpret_cast<refMat3>(F),
        one_over_h
    );
}

 }
}

#endif

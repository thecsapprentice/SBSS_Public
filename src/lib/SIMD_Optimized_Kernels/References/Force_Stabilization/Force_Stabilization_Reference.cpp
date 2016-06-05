//#####################################################################
//  Copyright (c) 2011-2013 Nathan Mitchell, Eftychios Sifakis.
//  This file is covered by the FreeBSD license. Please refer to the 
//  license.txt file for more information.
//#####################################################################


#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Arrays_Computations/MAGNITUDE.h>
#include <PhysBAM_Tools/Matrices/DIAGONAL_MATRIX_3X3.h>
#include <PhysBAM_Tools/Matrices/MATRIX_3X3.h>
#include <PhysBAM_Tools/Matrices/SYMMETRIC_MATRIX_3X3.h>
#include <PhysBAM_Tools/Matrices/MATRIX_MXN.h>
#include <PhysBAM_Tools/Math_Tools/RANGE.h>

#include "ROTATED_STRESS_DERIVATIVE.h"
#include "STENCIL.h"
#include "RANGE_ITERATOR.h"

using namespace PhysBAM;


namespace{


typedef float T;
const int d = 3;
typedef VECTOR<int,d> T_INDEX;
typedef STENCIL<T,d> T_STENCIL;

 float h;
 float stabilization_factor;

ARRAY<MATRIX_MXN<T> > G_HQ;
MATRIX_MXN<T> G_One;
MATRIX_MXN<T> K_stab;
VECTOR<T_STENCIL,d> cell_centered_derivative_operator;
ARRAY<T> cell_volume_weights;

template<class T_MATRIX>
void Print_Formatted(const T_MATRIX& A,std::ostream& output)
{
    output << "[" << std::endl;

    for(int i=1;i<=A.Rows();i++){
        for(int j=1;j<=A.Columns();j++){
            if(A.Valid_Index(i,j))
                output<<std::setw(12)<<A(i,j);
            else
                output<<"            ";
            if(j<A.Columns()) output<<", ";}
        output<<";"<<std::endl;}
    output << "]" << std::endl;
}

template<class T>
void Initialize_Cell_Centered_Quadrature()
{
    // Construct G_One (gradient) matrix for one-point quadrature

    G_One.Resize(d,8);

    int vertex=1;
    for(RANGE_ITERATOR<d> iterator(RANGE<T_INDEX>(T_INDEX(),T_INDEX::All_Ones_Vector()));iterator.Valid();iterator.Next(),vertex++){
        const T_INDEX& index=iterator.Index();

	for(int v=1;v<=d;v++)
	    if(d==2)
	        if(index(v)==0) G_One(v,vertex)=-1./(2.*h); else G_One(v,vertex)=1./(2.*h);
	    else
	        if(index(v)==0) G_One(v,vertex)=-1./(4.*h); else G_One(v,vertex)=1./(4.*h);
    }
}




template<class T>
void Initialize_Gauss_Quadrature()
{
    // Set up locations of quadrature points and cell volume corresponding to each quadrature point

    const int number_of_quadrature_points=8;
    MATRIX_MXN<T> quadrature_weights(d,number_of_quadrature_points);
    cell_volume_weights.Resize(number_of_quadrature_points);

    int quadrature_point=1;
    for(RANGE_ITERATOR<d> iterator(RANGE<T_INDEX>(T_INDEX(),T_INDEX::All_Ones_Vector()));iterator.Valid();iterator.Next(),quadrature_point++){
        cell_volume_weights(quadrature_point)=(d==2)?sqr(h)/((T)number_of_quadrature_points):(h*h*h)/((T)number_of_quadrature_points);
        const T_INDEX& index=iterator.Index();
        VECTOR<T,d> weights;
        for(int i=1;i<=d;i++)
	    if(index(i)==0) quadrature_weights(i,quadrature_point)=(1-one_over_root_three)*.5;
	    else quadrature_weights(i,quadrature_point)=(1+one_over_root_three)*.5;}

    // Construct G_HQ (gradient) matrix from Gauss quadrature points

    G_HQ.Resize(number_of_quadrature_points);

    for(int q=1;q<=number_of_quadrature_points;q++){
        G_HQ(q).Resize(d,8);
        VECTOR<T,d> weights;
        for(int v=1;v<=d;v++)
            weights(v)=quadrature_weights(v,q);

        int vertex=1;
	for(RANGE_ITERATOR<d> iterator(RANGE<T_INDEX>(T_INDEX(),T_INDEX::All_Ones_Vector()));iterator.Valid();iterator.Next(),vertex++){
            const T_INDEX& index=iterator.Index();

            for(int v=1;v<=d;v++){
                T& coefficient=G_HQ(q)(v,vertex);
                coefficient=1.;

                for(int w=1;w<=d;w++)
                    if(w==v)
                        if(index(w)==0) coefficient*=-1./h; else coefficient*=1./h;
                    else
		        if(index(w)==0) coefficient*=(1.-weights(w)); else coefficient*=weights(w);}}
    }
}




template<class T>
void Initialize_Stabilization_Kernel()
{
    const T cell_volume=(d==2)?sqr(h):(h*h*h);

    K_stab=G_One.Normal_Equations_Matrix()*cell_volume;
    for(int q=1;q<=G_HQ.m;q++)
        K_stab-=G_HQ(q).Normal_Equations_Matrix()*cell_volume_weights(q);

    K_stab*=(2.*stabilization_factor);
}

template void Initialize_Cell_Centered_Quadrature<float>();
template void Initialize_Gauss_Quadrature<float>();
template void Initialize_Stabilization_Kernel<float>();


}

template<class T>
void Force_Stabilization_Reference(const T Du[3][8], const T constant, /*const T h, const T stabilization_factor,*/
                                   /*const T weight,*/ T dH[3][8] )
{

    ::h = h;
    ::stabilization_factor = stabilization_factor;

    Initialize_Cell_Centered_Quadrature<float>();
    Initialize_Gauss_Quadrature<float>();
    Initialize_Stabilization_Kernel<float>();

    //    Print_Formatted(G_One,std::cout);

    //    std::cout << std::endl;
    //    std::cout << std::endl;

    //    for( int i = 0; i < 8 ; i++){
    //        Print_Formatted(G_HQ(i+1),std::cout);     std::cout << std::endl;}

    //    std::cout << std::endl;
    //    std::cout << std::endl;

    //    Print_Formatted(K_stab,std::cout);

    //    std::cout << std::endl;
    //    std::cout << std::endl;

    //    Print_Formatted(K_stab*24/h,std::cout);

    //    std::cout << std::endl;
    //    std::cout << std::endl;

    MATRIX_MXN<T> dDs;
    dDs.Resize(3,8);
    for(int i=0;i<3;i++) for(int j=0;j<8;j++) dDs(i+1,j+1)=Du[i][j];

    MATRIX_MXN<T> _dH =dDs*K_stab*constant;
    
    Print_Formatted(_dH,std::cout);


    for(int i=0;i<3;i++) for(int j=0;j<8;j++) dH[i][j] += _dH(i+1,j+1);
}



template<class T>
bool Force_Stabilization_Compare(const T dH[3][8],  const T dH_reference[3][8])
{

    MATRIX_MXN<T> _dH(3,8);
    for(int i=0;i<3;i++) for(int j=0;j<8;j++) _dH(i+1,j+1)=dH[i][j];
    MATRIX_MXN<T> _dH_reference(3,8);
    for(int i=0;i<3;i++) for(int j=0;j<8;j++) _dH_reference(i+1,j+1)=dH_reference[i][j];

    std::cout<<"Computed matrix dH :"<<std::endl;Print_Formatted(_dH,std::cout);
    std::cout<<"Reference matrix dH :"<<std::endl;Print_Formatted(_dH_reference,std::cout);
    std::cout<<"Difference = "<<(_dH-_dH_reference).Frobenius_Norm()<<std::endl;

    if( (_dH-_dH_reference).Frobenius_Norm() < 0.00001 )
        return true;
    else
        return false;



}


template void Force_Stabilization_Reference(const float Du[3][8], const float constant, /* const float h, const float stabilization_factor,*/ /*const float weight,*/ float dH[3][8] );
template bool Force_Stabilization_Compare(const float dH[3][8],  const float dH_reference[3][8]);
 

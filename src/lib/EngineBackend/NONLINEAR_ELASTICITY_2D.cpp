//#####################################################################
// Copyright 2011, Taylor Patterson, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class NONLINEAR_ELASTICITY (2D specializations)
//#####################################################################
#include "NONLINEAR_ELASTICITY.h"
#include <Common_Tools/Math_Tools/RANGE_ITERATOR.h>

#ifdef LOG_DETAILED_PERFORMANCE
#define LOG_DETAILED_PERFORMANCE_NE
#endif

using namespace PhysBAM;

template<class T,int d>
VECTOR<typename NONLINEAR_ELASTICITY<T,d>::T_SCALAR_VARIABLE_VIEW_CONST,d> NONLINEAR_ELASTICITY<T,d>::
View_Convert(const VECTOR<T_SCALAR_VARIABLE_VIEW,d>& in)
{
    return VECTOR<T_SCALAR_VARIABLE_VIEW_CONST,d>(in.x,in.y);
};

//#####################################################################
// Function Update_Position_Based_State
//#####################################################################
template<class T,int d> template<bool enable_constraints,bool enable_muscles> void NONLINEAR_ELASTICITY<T,d>::
Update_Position_Based_State_Specialized(T_VECTOR_VARIABLE_VIEW_CONST u, T_SCALAR_VARIABLE_VIEW_CONST p, T_VECTOR_VARIABLE_VIEW diag)
{PHYSBAM_NOT_IMPLEMENTED("Specialized Kernels only available in 3D.");}
//#####################################################################
// Function CompactData_Specialized
//#####################################################################
template<class T,int d> void NONLINEAR_ELASTICITY<T,d>::
CompactData_Specialized(T_VECTOR_VARIABLE_VIEW_CONST x, T_SCALAR_VARIABLE_VIEW_CONST p) const 
{PHYSBAM_NOT_IMPLEMENTED("Specialized Kernels only available in 3D.");}
//#####################################################################
// Function UnCompactData_Specialized
//#####################################################################
template<class T,int d> void NONLINEAR_ELASTICITY<T,d>::
UnCompactData_Specialized(T_VECTOR_VARIABLE_VIEW x,T_SCALAR_VARIABLE_VIEW p) const 
{PHYSBAM_NOT_IMPLEMENTED("Specialized Kernels only available in 3D.");}
//#####################################################################
// Function CompactData_Specialized
//#####################################################################
template<class T,int d> void NONLINEAR_ELASTICITY<T,d>::
CompactData_Specialized(T_VECTOR_VARIABLE_VIEW_CONST x) const 
{PHYSBAM_NOT_IMPLEMENTED("Specialized Kernels only available in 3D.");}
//#####################################################################
// Function UnCompactData_Specialized
//#####################################################################
template<class T,int d> void NONLINEAR_ELASTICITY<T,d>::
UnCompactData_Specialized(T_VECTOR_VARIABLE_VIEW x) const 
{PHYSBAM_NOT_IMPLEMENTED("Specialized Kernels only available in 3D.");}
//#####################################################################
// Function Update_Position_Based_State
//#####################################################################
template<class T,int d>  template<bool enable_constraints,bool enable_muscles> void NONLINEAR_ELASTICITY<T,d>::
Add_Force_First_Order_Elasticity_Specialized(T_VECTOR_VARIABLE_VIEW_CONST u,T_SCALAR_VARIABLE_VIEW_CONST p, T_VECTOR_VARIABLE_VIEW f,T_SCALAR_VARIABLE_VIEW q) 
{PHYSBAM_NOT_IMPLEMENTED("Specialized Kernels only available in 3D.");}
//#####################################################################
// Function Update_Position_Based_State
//#####################################################################
template<class T,int d>  template<bool enable_constraints,bool enable_muscles> void NONLINEAR_ELASTICITY<T,d>::
Add_Force_Differential_Elasticity_Specialized(T_VECTOR_VARIABLE_VIEW_CONST du,T_SCALAR_VARIABLE_VIEW_CONST dp,T_VECTOR_VARIABLE_VIEW df,T_SCALAR_VARIABLE_VIEW dq) const
{PHYSBAM_NOT_IMPLEMENTED("Specialized Kernels only available in 3D.");}
//#####################################################################
// Function Update_Position_Based_State
//#####################################################################
template<class T,int d>  template<bool enable_constraints,bool enable_muscles> void NONLINEAR_ELASTICITY<T,d>::
Add_Force_Differential_Elasticity_Specialized(T_VECTOR_VARIABLE_VIEW_CONST du,T_SCALAR_VARIABLE_VIEW_CONST dp,T_VECTOR_VARIABLE_VIEW df,T_SCALAR_VARIABLE_VIEW dq, const int subdomain) const
{PHYSBAM_NOT_IMPLEMENTED("Specialized Kernels only available in 3D.");}
//#####################################################################
// Function Initialize_Blocks
//#####################################################################
template<class T,int d> void NONLINEAR_ELASTICITY<T,d>::
Initialize_Blocks_Specialized(const int number_of_partitions)
{/*NOP*/}

//#####################################################################
template class NONLINEAR_ELASTICITY<float,2>;
template void NONLINEAR_ELASTICITY<float,2>::Update_Position_Based_State_Specialized<true,true>(T_VECTOR_VARIABLE_VIEW_CONST u,T_SCALAR_VARIABLE_VIEW_CONST p, T_VECTOR_VARIABLE_VIEW diag);
template void NONLINEAR_ELASTICITY<float,2>::Update_Position_Based_State_Specialized<true,false>(T_VECTOR_VARIABLE_VIEW_CONST u,T_SCALAR_VARIABLE_VIEW_CONST p, T_VECTOR_VARIABLE_VIEW diag);
template void NONLINEAR_ELASTICITY<float,2>::Update_Position_Based_State_Specialized<false,true>(T_VECTOR_VARIABLE_VIEW_CONST u,T_SCALAR_VARIABLE_VIEW_CONST p, T_VECTOR_VARIABLE_VIEW diag);
template void NONLINEAR_ELASTICITY<float,2>::Update_Position_Based_State_Specialized<false,false>(T_VECTOR_VARIABLE_VIEW_CONST u,T_SCALAR_VARIABLE_VIEW_CONST p, T_VECTOR_VARIABLE_VIEW diag);

template void NONLINEAR_ELASTICITY<float,2>::Add_Force_First_Order_Elasticity_Specialized<true,true>(T_VECTOR_VARIABLE_VIEW_CONST u,T_SCALAR_VARIABLE_VIEW_CONST p,T_VECTOR_VARIABLE_VIEW f,T_SCALAR_VARIABLE_VIEW q) ;
template void NONLINEAR_ELASTICITY<float,2>::Add_Force_First_Order_Elasticity_Specialized<true,false>(T_VECTOR_VARIABLE_VIEW_CONST u,T_SCALAR_VARIABLE_VIEW_CONST p,T_VECTOR_VARIABLE_VIEW f,T_SCALAR_VARIABLE_VIEW q) ;
template void NONLINEAR_ELASTICITY<float,2>::Add_Force_First_Order_Elasticity_Specialized<false,true>(T_VECTOR_VARIABLE_VIEW_CONST u,T_SCALAR_VARIABLE_VIEW_CONST p,T_VECTOR_VARIABLE_VIEW f,T_SCALAR_VARIABLE_VIEW q) ;
template void NONLINEAR_ELASTICITY<float,2>::Add_Force_First_Order_Elasticity_Specialized<false,false>(T_VECTOR_VARIABLE_VIEW_CONST u,T_SCALAR_VARIABLE_VIEW_CONST p,T_VECTOR_VARIABLE_VIEW f,T_SCALAR_VARIABLE_VIEW q) ;

template void NONLINEAR_ELASTICITY<float,2>::Add_Force_Differential_Elasticity_Specialized<true,true>(T_VECTOR_VARIABLE_VIEW_CONST du,T_SCALAR_VARIABLE_VIEW_CONST dp,T_VECTOR_VARIABLE_VIEW df,T_SCALAR_VARIABLE_VIEW dq) const;
template void NONLINEAR_ELASTICITY<float,2>::Add_Force_Differential_Elasticity_Specialized<true,false>(T_VECTOR_VARIABLE_VIEW_CONST du,T_SCALAR_VARIABLE_VIEW_CONST dp,T_VECTOR_VARIABLE_VIEW df,T_SCALAR_VARIABLE_VIEW dq) const;
template void NONLINEAR_ELASTICITY<float,2>::Add_Force_Differential_Elasticity_Specialized<false,true>(T_VECTOR_VARIABLE_VIEW_CONST du,T_SCALAR_VARIABLE_VIEW_CONST dp,T_VECTOR_VARIABLE_VIEW df,T_SCALAR_VARIABLE_VIEW dq) const;
template void NONLINEAR_ELASTICITY<float,2>::Add_Force_Differential_Elasticity_Specialized<false,false>(T_VECTOR_VARIABLE_VIEW_CONST du,T_SCALAR_VARIABLE_VIEW_CONST dp,T_VECTOR_VARIABLE_VIEW df,T_SCALAR_VARIABLE_VIEW dq) const;

template void NONLINEAR_ELASTICITY<float,2>::Add_Force_Differential_Elasticity_Specialized<true,true>(T_VECTOR_VARIABLE_VIEW_CONST du,T_SCALAR_VARIABLE_VIEW_CONST dp,T_VECTOR_VARIABLE_VIEW df,T_SCALAR_VARIABLE_VIEW dq, const int subdomain) const;
template void NONLINEAR_ELASTICITY<float,2>::Add_Force_Differential_Elasticity_Specialized<true,false>(T_VECTOR_VARIABLE_VIEW_CONST du,T_SCALAR_VARIABLE_VIEW_CONST dp,T_VECTOR_VARIABLE_VIEW df,T_SCALAR_VARIABLE_VIEW dq, const int subdomain) const;
template void NONLINEAR_ELASTICITY<float,2>::Add_Force_Differential_Elasticity_Specialized<false,true>(T_VECTOR_VARIABLE_VIEW_CONST du,T_SCALAR_VARIABLE_VIEW_CONST dp,T_VECTOR_VARIABLE_VIEW df,T_SCALAR_VARIABLE_VIEW dq, const int subdomain) const;
template void NONLINEAR_ELASTICITY<float,2>::Add_Force_Differential_Elasticity_Specialized<false,false>(T_VECTOR_VARIABLE_VIEW_CONST du,T_SCALAR_VARIABLE_VIEW_CONST dp,T_VECTOR_VARIABLE_VIEW df,T_SCALAR_VARIABLE_VIEW dq, const int subdomain) const;

#ifndef USE_SPECIALIZED_KERNELS
template class NONLINEAR_ELASTICITY<double,2>;
template void NONLINEAR_ELASTICITY<double,2>::Update_Position_Based_State_Specialized<true,true>(T_VECTOR_VARIABLE_VIEW_CONST u,T_SCALAR_VARIABLE_VIEW_CONST p, T_VECTOR_VARIABLE_VIEW diag);
template void NONLINEAR_ELASTICITY<double,2>::Update_Position_Based_State_Specialized<true,false>(T_VECTOR_VARIABLE_VIEW_CONST u,T_SCALAR_VARIABLE_VIEW_CONST p, T_VECTOR_VARIABLE_VIEW diag);
template void NONLINEAR_ELASTICITY<double,2>::Update_Position_Based_State_Specialized<false,true>(T_VECTOR_VARIABLE_VIEW_CONST u,T_SCALAR_VARIABLE_VIEW_CONST p, T_VECTOR_VARIABLE_VIEW diag);
template void NONLINEAR_ELASTICITY<double,2>::Update_Position_Based_State_Specialized<false,false>(T_VECTOR_VARIABLE_VIEW_CONST u,T_SCALAR_VARIABLE_VIEW_CONST p, T_VECTOR_VARIABLE_VIEW diag);

template void NONLINEAR_ELASTICITY<double,2>::Add_Force_First_Order_Elasticity_Specialized<true,true>(T_VECTOR_VARIABLE_VIEW_CONST u,T_SCALAR_VARIABLE_VIEW_CONST p,T_VECTOR_VARIABLE_VIEW f,T_SCALAR_VARIABLE_VIEW q) ;
template void NONLINEAR_ELASTICITY<double,2>::Add_Force_First_Order_Elasticity_Specialized<true,false>(T_VECTOR_VARIABLE_VIEW_CONST u,T_SCALAR_VARIABLE_VIEW_CONST p,T_VECTOR_VARIABLE_VIEW f,T_SCALAR_VARIABLE_VIEW q) ;
template void NONLINEAR_ELASTICITY<double,2>::Add_Force_First_Order_Elasticity_Specialized<false,true>(T_VECTOR_VARIABLE_VIEW_CONST u,T_SCALAR_VARIABLE_VIEW_CONST p,T_VECTOR_VARIABLE_VIEW f,T_SCALAR_VARIABLE_VIEW q) ;
template void NONLINEAR_ELASTICITY<double,2>::Add_Force_First_Order_Elasticity_Specialized<false,false>(T_VECTOR_VARIABLE_VIEW_CONST u,T_SCALAR_VARIABLE_VIEW_CONST p,T_VECTOR_VARIABLE_VIEW f,T_SCALAR_VARIABLE_VIEW q) ;

template void NONLINEAR_ELASTICITY<double,2>::Add_Force_Differential_Elasticity_Specialized<true,true>(T_VECTOR_VARIABLE_VIEW_CONST du,T_SCALAR_VARIABLE_VIEW_CONST dp,T_VECTOR_VARIABLE_VIEW df,T_SCALAR_VARIABLE_VIEW dq) const;
template void NONLINEAR_ELASTICITY<double,2>::Add_Force_Differential_Elasticity_Specialized<true,false>(T_VECTOR_VARIABLE_VIEW_CONST du,T_SCALAR_VARIABLE_VIEW_CONST dp,T_VECTOR_VARIABLE_VIEW df,T_SCALAR_VARIABLE_VIEW dq) const;
template void NONLINEAR_ELASTICITY<double,2>::Add_Force_Differential_Elasticity_Specialized<false,true>(T_VECTOR_VARIABLE_VIEW_CONST du,T_SCALAR_VARIABLE_VIEW_CONST dp,T_VECTOR_VARIABLE_VIEW df,T_SCALAR_VARIABLE_VIEW dq) const;
template void NONLINEAR_ELASTICITY<double,2>::Add_Force_Differential_Elasticity_Specialized<false,false>(T_VECTOR_VARIABLE_VIEW_CONST du,T_SCALAR_VARIABLE_VIEW_CONST dp,T_VECTOR_VARIABLE_VIEW df,T_SCALAR_VARIABLE_VIEW dq) const;

template void NONLINEAR_ELASTICITY<double,2>::Add_Force_Differential_Elasticity_Specialized<true,true>(T_VECTOR_VARIABLE_VIEW_CONST du,T_SCALAR_VARIABLE_VIEW_CONST dp,T_VECTOR_VARIABLE_VIEW df,T_SCALAR_VARIABLE_VIEW dq, const int subdomain) const;
template void NONLINEAR_ELASTICITY<double,2>::Add_Force_Differential_Elasticity_Specialized<true,false>(T_VECTOR_VARIABLE_VIEW_CONST du,T_SCALAR_VARIABLE_VIEW_CONST dp,T_VECTOR_VARIABLE_VIEW df,T_SCALAR_VARIABLE_VIEW dq, const int subdomain) const;
template void NONLINEAR_ELASTICITY<double,2>::Add_Force_Differential_Elasticity_Specialized<false,true>(T_VECTOR_VARIABLE_VIEW_CONST du,T_SCALAR_VARIABLE_VIEW_CONST dp,T_VECTOR_VARIABLE_VIEW df,T_SCALAR_VARIABLE_VIEW dq, const int subdomain) const;
template void NONLINEAR_ELASTICITY<double,2>::Add_Force_Differential_Elasticity_Specialized<false,false>(T_VECTOR_VARIABLE_VIEW_CONST du,T_SCALAR_VARIABLE_VIEW_CONST dp,T_VECTOR_VARIABLE_VIEW df,T_SCALAR_VARIABLE_VIEW dq, const int subdomain) const;
#endif

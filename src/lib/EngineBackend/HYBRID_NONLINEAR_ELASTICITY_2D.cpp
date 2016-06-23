#include "HYBRID_NONLINEAR_ELASTICITY.h"
#include <Common_Tools/Math_Tools/RANGE_ITERATOR.h>
using namespace PhysBAM;

#ifdef LOG_DETAILED_PERFORMANCE
#define LOG_DETAILED_PERFORMANCE_NE
#endif

template<class T_NONLINEAR_ELASTICITY>
HYBRID_NONLINEAR_ELASTICITY_STATE<T_NONLINEAR_ELASTICITY>::HYBRID_NONLINEAR_ELASTICITY_STATE() :
    x_mesh(T_SCALAR_VARIABLE_MESH_VIEW(0,NULL),T_SCALAR_VARIABLE_MESH_VIEW(0,NULL)), p_mesh(0,NULL)
{}

template<class T_NONLINEAR_ELASTICITY>
typename HYBRID_NONLINEAR_ELASTICITY<T_NONLINEAR_ELASTICITY>::T_VECTOR_VARIABLE_MESH_VIEW_CONST 
HYBRID_NONLINEAR_ELASTICITY<T_NONLINEAR_ELASTICITY>::
View_Convert(const VECTOR<T_SCALAR_VARIABLE_MESH_VIEW,d>& in)
{
    return VECTOR<T_SCALAR_VARIABLE_MESH_VIEW_CONST,d>(in.x,in.y);
};

//#####################################################################
// Function Initialize_Blocks_Elasticity
//#####################################################################
template<class T_NONLINEAR_ELASTICITY> void HYBRID_NONLINEAR_ELASTICITY<T_NONLINEAR_ELASTICITY>::
Initialize_Blocks_Elasticity(int number_of_partitions)
{/*NOP*/}

//#####################################################################
// Function CompactData_Specialized
//#####################################################################
template<class T_NONLINEAR_ELASTICITY> void HYBRID_NONLINEAR_ELASTICITY<T_NONLINEAR_ELASTICITY>::
CompactData_Specialized(T_VECTOR_VARIABLE_VIEW_CONST x, T_SCALAR_VARIABLE_VIEW_CONST p, T_VECTOR_VARIABLE_MESH_VIEW_CONST x_mesh, T_SCALAR_VARIABLE_MESH_VIEW_CONST p_mesh) const 
{PHYSBAM_NOT_IMPLEMENTED("Specialized Kernels only available in 3D.");}
//#####################################################################
// Function UnCompactData_Specialized
//#####################################################################
template<class T_NONLINEAR_ELASTICITY> void HYBRID_NONLINEAR_ELASTICITY<T_NONLINEAR_ELASTICITY>::
UnCompactData_Specialized(T_VECTOR_VARIABLE_VIEW x,T_SCALAR_VARIABLE_VIEW p, T_VECTOR_VARIABLE_MESH_VIEW x_mesh,T_SCALAR_VARIABLE_MESH_VIEW p_mesh) const 
{PHYSBAM_NOT_IMPLEMENTED("Specialized Kernels only available in 3D.");}
//#####################################################################
// Function CompactData_Specialized
//#####################################################################
template<class T_NONLINEAR_ELASTICITY> void HYBRID_NONLINEAR_ELASTICITY<T_NONLINEAR_ELASTICITY>::
CompactData_Specialized(T_VECTOR_VARIABLE_VIEW_CONST x, T_VECTOR_VARIABLE_MESH_VIEW_CONST x_mesh) const 
{PHYSBAM_NOT_IMPLEMENTED("Specialized Kernels only available in 3D.");}
//#####################################################################
// Function UnCompactData_Specialized
//#####################################################################
template<class T_NONLINEAR_ELASTICITY> void HYBRID_NONLINEAR_ELASTICITY<T_NONLINEAR_ELASTICITY>::
UnCompactData_Specialized(T_VECTOR_VARIABLE_VIEW x, T_VECTOR_VARIABLE_MESH_VIEW x_mesh) const 
{PHYSBAM_NOT_IMPLEMENTED("Specialized Kernels only available in 3D.");}
//#####################################################################
// Function Update_Position_Based_State
//#####################################################################
template<class T_NONLINEAR_ELASTICITY> void HYBRID_NONLINEAR_ELASTICITY<T_NONLINEAR_ELASTICITY>::
Update_Position_Based_State_Specialized(T_VECTOR_VARIABLE_VIEW_CONST u,T_SCALAR_VARIABLE_VIEW_CONST p,T_VECTOR_VARIABLE_VIEW diag,T_VECTOR_VARIABLE_MESH_VIEW_CONST u_mesh,T_SCALAR_VARIABLE_MESH_VIEW_CONST p_mesh, T_VECTOR_VARIABLE_MESH_VIEW diag_mesh)
{PHYSBAM_NOT_IMPLEMENTED("Specialized Kernels only available in 3D.");}
//#####################################################################
// Function Update_Position_Based_State
//#####################################################################
template<class T_NONLINEAR_ELASTICITY> void HYBRID_NONLINEAR_ELASTICITY<T_NONLINEAR_ELASTICITY>::
Add_Force_First_Order_Elasticity_Specialized(T_VECTOR_VARIABLE_VIEW_CONST u,T_SCALAR_VARIABLE_VIEW_CONST p,T_VECTOR_VARIABLE_MESH_VIEW_CONST u_mesh,T_SCALAR_VARIABLE_MESH_VIEW_CONST p_mesh, T_VECTOR_VARIABLE_VIEW f,T_SCALAR_VARIABLE_VIEW q,T_VECTOR_VARIABLE_MESH_VIEW f_mesh,T_SCALAR_VARIABLE_MESH_VIEW q_mesh) 
{PHYSBAM_NOT_IMPLEMENTED("Specialized Kernels only available in 3D.");}
//#####################################################################
// Function Update_Position_Based_State
//#####################################################################
template<class T_NONLINEAR_ELASTICITY> void HYBRID_NONLINEAR_ELASTICITY<T_NONLINEAR_ELASTICITY>::
Add_Force_Differential_Elasticity_Specialized(T_VECTOR_VARIABLE_VIEW_CONST du,T_SCALAR_VARIABLE_VIEW_CONST dp,T_VECTOR_VARIABLE_MESH_VIEW_CONST du_mesh,T_SCALAR_VARIABLE_MESH_VIEW_CONST dp_mesh,
                                              T_VECTOR_VARIABLE_VIEW df,T_SCALAR_VARIABLE_VIEW dq, T_VECTOR_VARIABLE_MESH_VIEW df_mesh,T_SCALAR_VARIABLE_MESH_VIEW dq_mesh) const
{PHYSBAM_NOT_IMPLEMENTED("Specialized Kernels only available in 3D.");}
//#####################################################################
// Function Update_Position_Based_State
//#####################################################################
template<class T_NONLINEAR_ELASTICITY> void HYBRID_NONLINEAR_ELASTICITY<T_NONLINEAR_ELASTICITY>::
Add_Force_Differential_Elasticity_Specialized(T_VECTOR_VARIABLE_VIEW_CONST du,T_SCALAR_VARIABLE_VIEW_CONST dp,T_VECTOR_VARIABLE_MESH_VIEW_CONST du_mesh,T_SCALAR_VARIABLE_MESH_VIEW_CONST dp_mesh,
                                              T_VECTOR_VARIABLE_VIEW df,T_SCALAR_VARIABLE_VIEW dq, T_VECTOR_VARIABLE_MESH_VIEW df_mesh,T_SCALAR_VARIABLE_MESH_VIEW dq_mesh, const int subdomain) const
{PHYSBAM_NOT_IMPLEMENTED("Specialized Kernels only available in 3D.");}

//#####################################################################
// Function Stress_Mesh
//#####################################################################
template<class T_NONLINEAR_ELASTICITY>
typename HYBRID_NONLINEAR_ELASTICITY<T_NONLINEAR_ELASTICITY>::TV
HYBRID_NONLINEAR_ELASTICITY<T_NONLINEAR_ELASTICITY>::
Stress_Mesh( const int& mesh_element, const TV& multilinear_coordinates )
{PHYSBAM_NOT_IMPLEMENTED("Stress/Strain access available only in 3D."); }
//#####################################################################
// Function Stress_Grid
//#####################################################################
template<class T_NONLINEAR_ELASTICITY>
typename HYBRID_NONLINEAR_ELASTICITY<T_NONLINEAR_ELASTICITY>::TV
HYBRID_NONLINEAR_ELASTICITY<T_NONLINEAR_ELASTICITY>::
Stress_Grid( const T_INDEX& cell_index, const TV& multilinear_coordinates )
{PHYSBAM_NOT_IMPLEMENTED("Stress/Strain access available only in 3D."); }
//#####################################################################
// Function Strain_Mesh
//#####################################################################
template<class T_NONLINEAR_ELASTICITY>
typename HYBRID_NONLINEAR_ELASTICITY<T_NONLINEAR_ELASTICITY>::TV
HYBRID_NONLINEAR_ELASTICITY<T_NONLINEAR_ELASTICITY>::
Strain_Mesh( const int& mesh_element, const TV& multilinear_coordinates )
{PHYSBAM_NOT_IMPLEMENTED("Stress/Strain access available only in 3D."); }
//#####################################################################
// Function Strain_Grid
//#####################################################################
template<class T_NONLINEAR_ELASTICITY>
typename HYBRID_NONLINEAR_ELASTICITY<T_NONLINEAR_ELASTICITY>::TV
HYBRID_NONLINEAR_ELASTICITY<T_NONLINEAR_ELASTICITY>::
Strain_Grid( const T_INDEX& cell_index, const TV& multilinear_coordinates )
{PHYSBAM_NOT_IMPLEMENTED("Stress/Strain access available only in 3D."); }

template struct HYBRID_NONLINEAR_ELASTICITY_STATE<NONLINEAR_ELASTICITY<float,2> >;
template struct HYBRID_NONLINEAR_ELASTICITY_STATE<SKINNING_NONLINEAR_ELASTICITY<float,2,true,true> >;
template struct HYBRID_NONLINEAR_ELASTICITY_STATE<SKINNING_NONLINEAR_ELASTICITY<float,2,true,false> >;
template class HYBRID_NONLINEAR_ELASTICITY<NONLINEAR_ELASTICITY<float,2> >;
template class HYBRID_NONLINEAR_ELASTICITY<SKINNING_NONLINEAR_ELASTICITY<float,2,true,true> >;
template class HYBRID_NONLINEAR_ELASTICITY<SKINNING_NONLINEAR_ELASTICITY<float,2,true,false> >;
#ifndef USE_SPECIALIZED_KERNELS
template struct HYBRID_NONLINEAR_ELASTICITY_STATE<NONLINEAR_ELASTICITY<double,2> >;
template struct HYBRID_NONLINEAR_ELASTICITY_STATE<SKINNING_NONLINEAR_ELASTICITY<double,2,true,true> >;
template struct HYBRID_NONLINEAR_ELASTICITY_STATE<SKINNING_NONLINEAR_ELASTICITY<double,2,true,false> >;
template class HYBRID_NONLINEAR_ELASTICITY<NONLINEAR_ELASTICITY<double,2> >;
template class HYBRID_NONLINEAR_ELASTICITY<SKINNING_NONLINEAR_ELASTICITY<double,2,true,true> >;
template class HYBRID_NONLINEAR_ELASTICITY<SKINNING_NONLINEAR_ELASTICITY<double,2,true,false> >;
#endif

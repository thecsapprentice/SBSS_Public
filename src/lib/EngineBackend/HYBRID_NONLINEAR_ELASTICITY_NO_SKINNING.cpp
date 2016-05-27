#include "HYBRID_NONLINEAR_ELASTICITY.h"
#include <Common/RANGE_ITERATOR.h>
using namespace PhysBAM;

#ifdef LOG_DETAILED_PERFORMANCE
#define LOG_DETAILED_PERFORMANCE_NE
#endif

//#####################################################################
// Function Initialize_Blocks_Constraints
//#####################################################################
template<class T_NONLINEAR_ELASTICITY> void HYBRID_NONLINEAR_ELASTICITY<T_NONLINEAR_ELASTICITY>::
Initialize_Blocks_Constraints(int number_of_partitions)
{/*NOP*/}
//#####################################################################
// Function Initialize_Blocks_Muscles
//#####################################################################
template<class T_NONLINEAR_ELASTICITY> void HYBRID_NONLINEAR_ELASTICITY<T_NONLINEAR_ELASTICITY>::
Initialize_Blocks_Muscles(int number_of_partitions)
{/*NOP*/}
//#####################################################################
// Function Initialize_Muscles
//#####################################################################
template<class T_NONLINEAR_ELASTICITY> void HYBRID_NONLINEAR_ELASTICITY<T_NONLINEAR_ELASTICITY>::
Initialize_Muscles()
{ /*NOP*/ }
//#####################################################################
// Function Update_Position_Based_State_Muscles_And_Constraints
//#####################################################################
template<class T_NONLINEAR_ELASTICITY> void HYBRID_NONLINEAR_ELASTICITY<T_NONLINEAR_ELASTICITY>::
Update_Position_Based_State_Muscles_And_Constraints(const T_STATE& state_u)
{/*NOP*/}
//#####################################################################
// Function Add_Force_Muscles_And_Constraints
//#####################################################################
template<class T_NONLINEAR_ELASTICITY> void HYBRID_NONLINEAR_ELASTICITY<T_NONLINEAR_ELASTICITY>::
Add_Force_Muscles_And_Constraints(const T_STATE& state_u, T_STATE& state_f) const
{/*NOP*/}    
//#####################################################################
// Function Add_Force_Differential_Muscles_And_Constraints
//#####################################################################
template<class T_NONLINEAR_ELASTICITY> void HYBRID_NONLINEAR_ELASTICITY<T_NONLINEAR_ELASTICITY>::
Add_Force_Differential_Muscles_And_Constraints(const T_STATE& state_du, T_STATE& state_df) const
{/*NOP*/}
//#####################################################################
// Function Cell_Has_Constraint
//#####################################################################
template<class T_NONLINEAR_ELASTICITY> bool HYBRID_NONLINEAR_ELASTICITY<T_NONLINEAR_ELASTICITY>::
Cell_Has_Constraint( const T_INDEX& index, int mesh_index ){
    return false;
}

template class HYBRID_NONLINEAR_ELASTICITY<NONLINEAR_ELASTICITY<float,3> >;
template class HYBRID_NONLINEAR_ELASTICITY<NONLINEAR_ELASTICITY<float,2> >;
#ifndef USE_SPECIALIZED_KERNELS
template class HYBRID_NONLINEAR_ELASTICITY<NONLINEAR_ELASTICITY<double,3> >;
template class HYBRID_NONLINEAR_ELASTICITY<NONLINEAR_ELASTICITY<double,2> >;
#endif

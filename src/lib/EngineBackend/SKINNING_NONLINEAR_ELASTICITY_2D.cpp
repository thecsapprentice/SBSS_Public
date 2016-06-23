#include "SKINNING_NONLINEAR_ELASTICITY.h"
#include <Common_Tools/Math_Tools/RANGE_ITERATOR.h>

using namespace PhysBAM;

//#####################################################################
// Function Initialize_Blocks_Constraints
//#####################################################################
template<class T,int d,bool enable_constraints,bool enable_muscles> void SKINNING_NONLINEAR_ELASTICITY<T,d,enable_constraints,enable_muscles>::
Initialize_Blocks_Constraints(const int number_of_partitions)
{/*NOP*/}
//#####################################################################
// Function Initialize_Blocks_Muscles
//#####################################################################
template<class T,int d,bool enable_constraints,bool enable_muscles> void SKINNING_NONLINEAR_ELASTICITY<T,d,enable_constraints,enable_muscles>::
Initialize_Blocks_Muscles(const int number_of_partitions)
{/*NOP*/}

template class SKINNING_NONLINEAR_ELASTICITY<float,2,true,true>;
template class SKINNING_NONLINEAR_ELASTICITY<float,2,true,false>;
#ifndef USE_SPECIALIZED_KERNELS
template class SKINNING_NONLINEAR_ELASTICITY<double,2,true,true>;
template class SKINNING_NONLINEAR_ELASTICITY<double,2,true,false>;
#endif

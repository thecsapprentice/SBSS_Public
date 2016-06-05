#include "OPENVDB_LEVELSET_COLLISION.h"
#include <PhysBAM_Tools/Log/DEBUG_UTILITIES.h>
using namespace PhysBAM;

template<class T, int d>
OPENVDB_LEVELSET_COLLISION<T,d>::OPENVDB_LEVELSET_COLLISION(const std::vector<std::vector<float> > vertices, const std::vector<std::vector<int> > triangles, const float refinement){}
template<class T, int d> T 
OPENVDB_LEVELSET_COLLISION<T,d>::Phi_Implementation( const T world_location[d]) const
{PHYSBAM_NOT_IMPLEMENTED();};
template<class T, int d> void
OPENVDB_LEVELSET_COLLISION<T,d>::Normal_Implementation( const T world_location[d], T normal[d]) const
{PHYSBAM_NOT_IMPLEMENTED();};

template class OPENVDB_LEVELSET_COLLISION<float, 2>;
template class OPENVDB_LEVELSET_COLLISION<double, 2>;

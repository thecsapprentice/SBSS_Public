#include "OPENVDB_LEVELSET_COLLISION.h"
#include <PhysBAM_Tools/Log/DEBUG_UTILITIES.h>
#include <vector>
//#include <openvdb/tools/GridOperators.h>
//#include <openvdb/tools/MeshToVolume.h>
//#include <openvdb/util/NullInterrupter.h>
using namespace PhysBAM;

template<class T, int d>
OPENVDB_LEVELSET_COLLISION<T,d>::OPENVDB_LEVELSET_COLLISION(const std::vector<std::vector<float> > vertices, const std::vector<std::vector<int> > triangles, const float refinement)
{
/*
    openvdb::initialize();

    std::vector< openvdb::Vec3s > _points;
    std::vector< openvdb::Vec3I > _triangles;
    std::vector< openvdb::Vec4I > _quads;

    for( std::vector<std::vector<float> >::const_iterator iter = vertices.begin();
         iter != vertices.end(); 
         iter++){
        openvdb::Vec3s vertex;
        for(int i=0; i<3; i++)
            vertex[i] = (*iter)[i];
        _points.push_back( vertex );
    }

    for( std::vector<std::vector<int> >::const_iterator iter = triangles.begin();
         iter != triangles.end(); 
         iter++){
        openvdb::Vec3I triangle;
        for(int i=0; i<3; i++)
            triangle[i] = (*iter)[i];
        _triangles.push_back( triangle );
    }

    float exBandWidth = 3;
    float inBandWidth = 3;
    
    xtrans = openvdb::math::Transform::createLinearTransform(refinement);

    phi_grid = openvdb::tools::meshToSignedDistanceField<PhiGrid>( *xtrans,
                                          _points,
                                          _triangles,
                                          _quads,
                                          exBandWidth,
                                          inBandWidth);

    phi_grid->setGridClass(openvdb::GRID_LEVEL_SET);
    phi_grid->setName("CollisionLevelSet");
    normal_grid = openvdb::tools::gradient<PhiGrid>( *phi_grid, true );
*/
};

template<class T, int d> T 
OPENVDB_LEVELSET_COLLISION<T,d>::Phi_Implementation( const T world_location[d]) const
{
    PHYSBAM_NOT_IMPLEMENTED();
/*
    typename PhiGrid::Accessor accessor = phi_grid->getAccessor();
    return accessor.getValue(xtrans->worldToIndexCellCentered(WorldCoord(world_location)));
*/
};

template<class T, int d> void
OPENVDB_LEVELSET_COLLISION<T,d>::Normal_Implementation( const T world_location[d], T normal[d]) const
{
    PHYSBAM_NOT_IMPLEMENTED();
/*
    typename NormalGrid::Accessor accessor = normal_grid->getAccessor();
    NormalType normal_vdb;
    normal_vdb = accessor.getValue(xtrans->worldToIndexCellCentered(WorldCoord(world_location)));
    normal[0] = normal_vdb.x();  normal[1] = normal_vdb.y();  normal[2] = normal_vdb.z();     
*/
};

template class OPENVDB_LEVELSET_COLLISION<float, 3>;
template class OPENVDB_LEVELSET_COLLISION<double, 3>;

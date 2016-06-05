#ifndef __WRITE_OUTPUT_EMBEDDED_H__
#define __WRITE_OUTPUT_EMBEDDED_H__

#include <EngineInterface/CELL_TYPE.h>
#include <PhysBAM_Geometry/Geometry_Particles/REGISTER_GEOMETRY_READ_WRITE.h>
#include <PhysBAM_Geometry/Solids_Geometry/DEFORMABLE_GEOMETRY_COLLECTION.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/FREE_PARTICLES.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/SEGMENTED_CURVE.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TETRAHEDRALIZED_VOLUME.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <PhysBAM_Tools/Read_Write/Arrays/READ_WRITE_ARRAY.h>
#include <PhysBAM_Tools/Read_Write/Grids_Uniform_Arrays/READ_WRITE_ARRAYS.h>

#include <Common/Write_Output.h>

#include "ELASTIC_LATTICE_DEFORMER.h"

class CLElib;

namespace PhysBAM{
    void Create_Output_Data_Embedded(const ELASTIC_LATTICE_DEFORMER& deformer,
                                     GEOMETRY_PARTICLES< VECTOR<float,3> >& particles,
                                     ARRAY<STRUCTURE<VECTOR<float,3> >* >& collection_structures);

    void Create_Output_Data_FineGrid(const ELASTIC_LATTICE_DEFORMER& deformer,
                                     GEOMETRY_PARTICLES< VECTOR<float,3> >& particles,
                                     ARRAY<STRUCTURE<VECTOR<float,3> >* >& collection_structures);

    void Write_Output_Embedded(STREAM_TYPE stream_type, const ELASTIC_LATTICE_DEFORMER& deformer, const std::string directory,const int frame);


}

#endif

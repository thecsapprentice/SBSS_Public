#ifndef __WRITE_OUTPUT_COMMON_H__
#define __WRITE_OUTPUT_COMMON_H__

#include <PhysBAM_Geometry/Geometry_Particles/REGISTER_GEOMETRY_READ_WRITE.h>
#include <PhysBAM_Geometry/Solids_Geometry/DEFORMABLE_GEOMETRY_COLLECTION.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/FREE_PARTICLES.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/SEGMENTED_CURVE.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TETRAHEDRALIZED_VOLUME.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <PhysBAM_Tools/Read_Write/Arrays/READ_WRITE_ARRAY.h>
#include <PhysBAM_Tools/Read_Write/Grids_Uniform_Arrays/READ_WRITE_ARRAYS.h>

class CLElib;

namespace PhysBAM{

    template< class T, int d >
        void Write_Output_Structures( STREAM_TYPE stream_type,
                                      GEOMETRY_PARTICLES< VECTOR<T,d> >& particles,
                                      ARRAY<STRUCTURE<VECTOR<T,d> >* >& collection_structures,
                                      const std::string directory,const int frame ){

        LOG::SCOPE scope( "Writing structures to disk..." );

        FILE_UTILITIES::Create_Directory(directory+"/"+STRING_UTILITIES::Value_To_String(frame));
        typedef VECTOR<T,d> TV;

        DEFORMABLE_GEOMETRY_COLLECTION<TV> collection(particles);
        for(int i = 1; i <= collection_structures.m; i++)
            collection.Add_Structure( collection_structures(i) );
        
        collection.Write(stream_type,directory,frame,frame,true);
    }


}

#endif

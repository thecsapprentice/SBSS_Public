#ifndef __WRITE_OUTPUT_LATTICE_H__
#define __WRITE_OUTPUT_LATTICE_H__

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

class CLElib;

namespace PhysBAM{

//  template<typename T_ELASTICITY>
//      void Write_Output(STREAM_TYPE stream_type,CLElib& deformer,const std::string directory,const int frame, const int impl);
    

    template<typename T_ELASTICITY>
        void Create_Output_Data_Lattice(const T_ELASTICITY& elasticity,
                                        const typename T_ELASTICITY::T_STATE& state,
                                        GEOMETRY_PARTICLES< VECTOR<typename T_ELASTICITY::SCALAR, T_ELASTICITY::DIM> >& particles,
                                        ARRAY<STRUCTURE<VECTOR<typename T_ELASTICITY::SCALAR, T_ELASTICITY::DIM> >* >& collection_structures);

//#####################################################################
//
//                      Write_Output_Minimal
//
//#####################################################################


    template<typename T_ELASTICITY>
        void Write_Output_Minimal(STREAM_TYPE stream_type, const T_ELASTICITY& elasticity, const typename T_ELASTICITY::T_STATE& state,const std::string directory,const int frame, const int impl){
  
        LOG::SCOPE scope( "Write Output Lattice");

        typedef typename T_ELASTICITY::SCALAR T;
        static const int d=T_ELASTICITY::DIM;
        
        typedef VECTOR<T,d> TV;
        typedef VECTOR<int,d> T_INDEX;
        
        Initialize_Geometry_Particle();Initialize_Read_Write_Structures();    
        GEOMETRY_PARTICLES<TV> particles;
        ARRAY<STRUCTURE<TV>* > collection_structures;
        
        Create_Output_Data_Lattice(elasticity, state, particles, collection_structures);

        for(int s=1;s<=collection_structures.m;s++)
            collection_structures(s)->Update_Number_Nodes();

        Write_Output_Structures(stream_type, particles, collection_structures, directory, frame);
    }


}

#endif

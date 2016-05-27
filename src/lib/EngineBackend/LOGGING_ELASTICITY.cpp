#include "LOGGING_ELASTICITY.h"

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>
#include <PhysBAM_Tools/Data_Structures/HASHTABLE.h>

#include <PhysBAM_Tools/Read_Write/Grids_Uniform_Arrays/READ_WRITE_ARRAYS.h>
#include <PhysBAM_Tools/Read_Write/Arrays/READ_WRITE_ARRAY.h>
#include <PhysBAM_Tools/Read_Write/Vectors/READ_WRITE_VECTOR_3D.h>

#include "Archive/ARCHIVE.h"
#include <stdio.h>
#include <string.h>

using namespace PhysBAM;


template<class TYPE> void
Write_Record( std::string name, const TYPE& data)
{
    LOG::SCOPE scope("Writing " + name + " to archive.");

    DATA_BUFFER temp;
    std::ostringstream stream;
    Write_Binary<float, TYPE >(stream, data);
    
    // Fast Copy from stringstream to PhysBAM Array.
    temp.Resize(stream.str().length());
    const char* ss_buf = stream.str().c_str();
    char* d_buf = temp.Get_Array_Pointer();
    memcpy( d_buf, ss_buf, stream.str().length() );
    // End fast copy.
    
    ARCHIVE::WriteDataBuffer(name,temp);
}


template<class TYPE> void
Read_Record( std::string name, TYPE& data)
{
    LOG::SCOPE scope("Reading " + name + " from archive.");

    DATA_BUFFER temp;
    ARCHIVE::ReadDataBuffer(name,temp);
    std::string d;
    d.resize(temp.m);
    for(int i=0;i<temp.m;i++){
        d[i] = temp(i+1);
    }
    std::istringstream stream(d);
    Read_Binary<float, TYPE >(stream, data);    

}


template<class T_ELASTICITY> void 
LOGGING_ELASTICITY<T_ELASTICITY>::
ExportElasticityData(T_DISCRETIZATION& discretization, const std::string& filename)
{
    LOG::SCOPE scope("ExportElasticityData::Exporting Sim Data");
    DATA_BUFFER temp;

    ARCHIVE::CreateArchive(filename);

    Write_Record("h", discretization.h);
    Write_Record("n", discretization.n);
    Write_Record("origin", discretization.grid.domain.min_corner);
    Write_Record("mu", discretization.constant_mu);
    Write_Record("kappa", discretization.constant_kappa);
    Write_Record("alpha", discretization.constant_alpha);
    Write_Record("cutoff", discretization.cutoff_value);
    Write_Record("stabilization", discretization.stabilization_factor);
    Write_Record("cell_type", discretization.cell_type);
    Write_Record("number_of_mesh_cells", discretization.Number_Of_Mesh_Cells());
    Write_Record("number_of_mesh_nodes", discretization.Number_Of_Mesh_Nodes());
    Write_Record("cell_type_mesh", discretization.cell_type_mesh);
    Write_Record("cell_indices_mesh", discretization.cell_indices_mesh);
    Write_Record("cells_mesh", discretization.cells_mesh);

    // Close
    ARCHIVE::CloseArchive();

    // Write to Disk
    ARCHIVE::WriteArchive();

    //Cleanup and Release resources
    ARCHIVE::CleanUpArchive();   
}

template<class T_ELASTICITY> void 
LOGGING_ELASTICITY<T_ELASTICITY>::
ImportElasticityData(const std::string& filename,
                     float& h, VECTOR<int,3>& n,
                     VECTOR<float,3>& origin,
                     float& mu,
                     float& alpha,
                     float& kappa,
                     float& cutoff_value,
                     float& stabilization_factor,
                     ARRAY<CELL_TYPE, VECTOR<int, 3> >& cell_type,
                     int& number_of_mesh_cells, 
                     int& number_of_mesh_nodes,
                     ARRAY<CELL_TYPE>& cell_type_mesh,
                     ARRAY<VECTOR<int,3> >& cell_indices_mesh,
                     ARRAY<VECTOR<int,8> >& cells_mesh )
{

    LOG::SCOPE scope("ImportElasticityData::Importing Sim Data");
    ARCHIVE::OpenArchive(filename);

    Read_Record("h", h);
    Read_Record("n", n);
    Read_Record("origin", origin);
    Read_Record("mu", mu);
    Read_Record("kappa", kappa);
    Read_Record("alpha", alpha);
    Read_Record("cutoff", cutoff_value);
    Read_Record("stabilization", stabilization_factor);
    Read_Record("cell_type", cell_type);
    Read_Record("number_of_mesh_cells", number_of_mesh_cells);
    Read_Record("number_of_mesh_nodes", number_of_mesh_nodes);
    Read_Record("cell_type_mesh", cell_type_mesh);
    Read_Record("cell_indices_mesh", cell_indices_mesh);
    Read_Record("cells_mesh", cells_mesh);

    // Close
    ARCHIVE::CloseArchive();

    //Cleanup and Release resources
    ARCHIVE::CleanUpArchive();   
}


template class LOGGING_ELASTICITY<HYBRID_NONLINEAR_ELASTICITY<SKINNING_NONLINEAR_ELASTICITY<float,3,true,false> > >;

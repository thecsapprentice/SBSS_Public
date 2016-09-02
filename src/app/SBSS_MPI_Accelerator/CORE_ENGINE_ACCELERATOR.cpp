#include <mpi.h>
#include <cstring>
#include "CORE_ENGINE_ACCELERATOR.h"

#include <iostream>
#include <PhysBAM_Tools/Read_Write/Grids_Uniform_Arrays/READ_WRITE_ARRAYS.h>
#include <PhysBAM_Tools/Read_Write/Arrays/READ_WRITE_ARRAY.h>
#include <PhysBAM_Tools/Read_Write/Vectors/READ_WRITE_VECTOR_3D.h>
#include <PhysBAM_Tools/Read_Write/Data_Structures/READ_WRITE_TRIPLE.h>
#include <PhysBAM_Tools/Data_Structures/TRIPLE.h>


#define COMMAND_TAG 1
#define DATA_TAG 2

namespace PhysBAM{

template<class TYPE> void
Marshal( const TYPE& data, ARRAY<unsigned char>& bytes )
{
    //LOG::SCOPE scope("Writing to byte stream.");

    std::ostringstream stream;
    Write_Binary<float, TYPE >(stream, data);
    
    // Fast Copy from stringstream to PhysBAM Array.
    bytes.Resize(stream.str().length());
    std::string s = stream.str();
    const char* ss_buf = s.c_str();
    unsigned char* d_buf = bytes.Get_Array_Pointer();
    memcpy( d_buf, ss_buf, stream.str().length() );
    // End fast copy.
}


template<class TYPE> void
UnMarshal( TYPE& data, const ARRAY<unsigned char>& bytes)
{
    //LOG::SCOPE scope("Reading from byte stream.");

    std::string d;
    d.resize(bytes.m);
    for(int i=0;i<bytes.m;i++){
        d[i] = bytes(i+1);
    }
    std::istringstream stream(d);
    Read_Binary<float, TYPE >(stream, data);    
}


void Read_Data( ARRAY<unsigned char>& bytes ){
    int buffer_size;
    int m;
    MPI_Status status;
    MPI_Recv(&(buffer_size),1,MPI_INT,0,DATA_TAG,MPI_COMM_WORLD,&status);
    MPI_Recv(&(m),1,MPI_INT,0,DATA_TAG,MPI_COMM_WORLD,&status);
    bytes.Exact_Resize(buffer_size, false);
    MPI_Recv(bytes.base_pointer,buffer_size,MPI_BYTE,0,DATA_TAG,MPI_COMM_WORLD,&status);
    bytes.m = m;
}

void Write_Data( ARRAY<unsigned char>& bytes ){
    MPI_Send(&(bytes.buffer_size),1,MPI_INT,0,DATA_TAG,MPI_COMM_WORLD);
    MPI_Send(&(bytes.m),1,MPI_INT,0,DATA_TAG,MPI_COMM_WORLD);
    MPI_Send(bytes.base_pointer,bytes.buffer_size,MPI_BYTE,0,DATA_TAG,MPI_COMM_WORLD);
}
        



    CORE_ENGINE_ACCELERATOR::CORE_ENGINE_ACCELERATOR(int partition_count_input)
    {
        int flag;
        MPI_Initialized(&flag);
        if( !flag ){
            LOG::cout << "CORE ENGINE FAILURE! No MPI Comm is available. " << std::endl;
            exit(1);
        }
    }

    CORE_ENGINE_ACCELERATOR::~CORE_ENGINE_ACCELERATOR()
    {
    }


void CORE_ENGINE_ACCELERATOR::
CreateEngine(){
    // Get data from MPI connection
    T_INDEX bounds;
    T dx;
    TV min_corner;
    T mu, lambda, alpha, cutoff_value, stabilization;
    MPI_Status status;

    MPI_Recv(&(bounds.x),3,MPI_INT,0,DATA_TAG,MPI_COMM_WORLD,&status);
    MPI_Recv(&(dx),1,MPI_FLOAT,0,DATA_TAG,MPI_COMM_WORLD,&status);
    MPI_Recv(&(min_corner.x),3,MPI_FLOAT,0,DATA_TAG,MPI_COMM_WORLD,&status);
    MPI_Recv(&(mu),1,MPI_FLOAT,0,DATA_TAG,MPI_COMM_WORLD,&status);
    MPI_Recv(&(lambda),1,MPI_FLOAT,0,DATA_TAG,MPI_COMM_WORLD,&status);
    MPI_Recv(&(alpha),1,MPI_FLOAT,0,DATA_TAG,MPI_COMM_WORLD,&status);
    MPI_Recv(&(cutoff_value),1,MPI_FLOAT,0,DATA_TAG,MPI_COMM_WORLD,&status);
    MPI_Recv(&(stabilization),1,MPI_FLOAT,0,DATA_TAG,MPI_COMM_WORLD,&status);
    
    engineInterface.CreateEngine( bounds, dx, min_corner, mu, lambda, alpha, cutoff_value, stabilization);
}

void  CORE_ENGINE_ACCELERATOR::
Initialize_Muscles(){
    engineInterface.Initialize_Muscles();
}

void  CORE_ENGINE_ACCELERATOR::
Initialize_Mesh(){
    // Get data from MPI connection
    int mesh_cell_count;
    int mesh_node_count;
    MPI_Status status;

    MPI_Recv(&(mesh_cell_count),1,MPI_INT,0,DATA_TAG,MPI_COMM_WORLD,&status);
    MPI_Recv(&(mesh_node_count),1,MPI_INT,0,DATA_TAG,MPI_COMM_WORLD,&status);
    
    engineInterface.Initialize_Mesh(mesh_cell_count, mesh_node_count);
}

void  CORE_ENGINE_ACCELERATOR::
InitializeEngine(){
    engineInterface.InitializeEngine();
}

void  CORE_ENGINE_ACCELERATOR::
Exact_Solve(){
    static int frame = 0;

    // Get data from MPI connection
    int  krylov_iterations;
    int  newton_iterations;
    T    krylov_tolerance; 
    T    newton_tolerance; 
    bool no_cut_cells;     
    MPI_Status status;              
    
    MPI_Recv(&(krylov_iterations),1,MPI_INT,0,DATA_TAG,MPI_COMM_WORLD,&status);
    MPI_Recv(&(newton_iterations),1,MPI_INT,0,DATA_TAG,MPI_COMM_WORLD,&status);
    MPI_Recv(&(krylov_tolerance),1,MPI_FLOAT,0,DATA_TAG,MPI_COMM_WORLD,&status);
    MPI_Recv(&(newton_tolerance),1,MPI_FLOAT,0,DATA_TAG,MPI_COMM_WORLD,&status);
    MPI_Recv(&(no_cut_cells),1,MPI_C_BOOL,0,DATA_TAG,MPI_COMM_WORLD,&status);

    T norm;
    engineInterface.Exact_Solve( krylov_iterations, newton_iterations, krylov_tolerance,
                                 newton_tolerance, no_cut_cells, norm );

    MPI_Send(&(norm),1,MPI_FLOAT,0,DATA_TAG,MPI_COMM_WORLD);
}

void  CORE_ENGINE_ACCELERATOR::
DestroyEngine(){
    engineInterface.DestroyEngine();
}

void  CORE_ENGINE_ACCELERATOR::
Displacement_Grid() const {
    // Get data from MPI connection
    T_INDEX cell;
    TV weights;
    TV result;
    MPI_Status status;              

    MPI_Recv(&(cell.x),3,MPI_INT,0,DATA_TAG,MPI_COMM_WORLD,&status);
    MPI_Recv(&(weights.x),3,MPI_FLOAT,0,DATA_TAG,MPI_COMM_WORLD,&status);

    result = engineInterface.Displacement_Grid( cell, weights);
    
    // Send data back
    MPI_Send(&(result.x),3,MPI_FLOAT,0,DATA_TAG,MPI_COMM_WORLD);
}

void  CORE_ENGINE_ACCELERATOR::
Displacement_Mesh() const {
    // Get data from MPI connection
    int cell;
    TV weights;
    TV result;
    MPI_Status status;              

    MPI_Recv(&(cell),1,MPI_INT,0,DATA_TAG,MPI_COMM_WORLD,&status);
    MPI_Recv(&(weights.x),3,MPI_FLOAT,0,DATA_TAG,MPI_COMM_WORLD,&status);

    result = engineInterface.Displacement_Mesh(cell, weights );

    // Send data back
    MPI_Send(&(result.x),3,MPI_FLOAT,0,DATA_TAG,MPI_COMM_WORLD);
}

void  CORE_ENGINE_ACCELERATOR::
Displacement() const {
    ARRAY<TRIPLE<int, T_INDEX, TV> > queue;
    ARRAY<TV> updates;

    ARRAY<unsigned char> bytes;
    Read_Data( bytes );
    UnMarshal< ARRAY<TRIPLE<int, T_INDEX, TV> > >(queue,bytes);
    updates.Resize( queue.m );

    engineInterface.Displacement( queue, updates );

    Marshal< ARRAY<TV> >(updates,bytes);
    Write_Data(bytes);    
}

void  CORE_ENGINE_ACCELERATOR::
Deformation() const {
    ARRAY<TRIPLE<int, T_INDEX, TV> > queue;
    ARRAY<TV> updates;

    ARRAY<unsigned char> bytes;
    Read_Data( bytes );
    UnMarshal< ARRAY<TRIPLE<int, T_INDEX, TV> > >(queue,bytes);
    updates.Resize( queue.m );

    engineInterface.Deformation( queue, updates );

    Marshal< ARRAY<TV> >(updates,bytes);
    Write_Data(bytes);    
}

void  CORE_ENGINE_ACCELERATOR::
Stress() const {
    ARRAY<TRIPLE<int, T_INDEX, TV> > queue;
    ARRAY<TV> updates;
    ARRAY<unsigned char> bytes;
    Read_Data( bytes );
    UnMarshal< ARRAY<TRIPLE<int, T_INDEX, TV> > >(queue,bytes);
    updates.Resize( queue.m );
    engineInterface.Stress( queue, updates );
    Marshal< ARRAY<TV> >(updates,bytes);
    Write_Data(bytes);    
}

void  CORE_ENGINE_ACCELERATOR::
Strain() const {
    ARRAY<TRIPLE<int, T_INDEX, TV> > queue;
    ARRAY<TV> updates;
    ARRAY<unsigned char> bytes;
    Read_Data( bytes );
    UnMarshal< ARRAY<TRIPLE<int, T_INDEX, TV> > >(queue,bytes);
    updates.Resize( queue.m );
    engineInterface.Strain( queue, updates );
    Marshal< ARRAY<TV> >(updates,bytes);
    Write_Data(bytes);    
}

void   CORE_ENGINE_ACCELERATOR::
h() const {
    T result = engineInterface.h();    
    // Send data back
    MPI_Send(&(result),1,MPI_FLOAT,0,DATA_TAG,MPI_COMM_WORLD);
}

void  CORE_ENGINE_ACCELERATOR::
Cell() const {
    // Get data from MPI connection
    TV location;
    T_INDEX result;
    MPI_Status status;      

    MPI_Recv(&(location.x),3,MPI_FLOAT,0,DATA_TAG,MPI_COMM_WORLD,&status);
    
    result = engineInterface.Cell( location );

    // Send data back
    MPI_Send(&(result.x),3,MPI_INT,0,DATA_TAG,MPI_COMM_WORLD);
}

void CORE_ENGINE_ACCELERATOR::
Node() const {
    // Get data from MPI connection
    T_INDEX location;
    TV result;
    MPI_Status status;      

    MPI_Recv(&(location.x),3,MPI_INT,0,DATA_TAG,MPI_COMM_WORLD,&status);
    
    result = engineInterface.Node(location);
    
    // Send data back
    MPI_Send(&(result.x),3,MPI_FLOAT,0,DATA_TAG,MPI_COMM_WORLD);
}

void  CORE_ENGINE_ACCELERATOR::
Padded_Node_Domain() const {
    RANGE<T_INDEX> result = engineInterface.Padded_Node_Domain();
    // Send data back
    MPI_Send(&(result.min_corner.x),6,MPI_INT,0,DATA_TAG,MPI_COMM_WORLD);           
}

void  CORE_ENGINE_ACCELERATOR::
Unpadded_Node_Domain() const {
    RANGE<T_INDEX> result = engineInterface.Unpadded_Node_Domain();
    // Send data back
    MPI_Send(&(result.min_corner.x),6,MPI_INT,0,DATA_TAG,MPI_COMM_WORLD);   
}

void  CORE_ENGINE_ACCELERATOR::
Padded_Cell_Domain() const {
    RANGE<T_INDEX> result = engineInterface.Padded_Cell_Domain();
    // Send data back
    MPI_Send(&(result.min_corner.x),6,MPI_INT,0,DATA_TAG,MPI_COMM_WORLD);   
}

void  CORE_ENGINE_ACCELERATOR::
Unpadded_Cell_Domain() const {
    RANGE<T_INDEX> result = engineInterface.Unpadded_Cell_Domain();
    // Send data back
    MPI_Send(&(result.min_corner.x),6,MPI_INT,0,DATA_TAG,MPI_COMM_WORLD);       
}

void  CORE_ENGINE_ACCELERATOR::
Grid() const {
    GRID<TV> result = engineInterface.Grid();
    // Send data back
    MPI_Send(&(result.counts.x),80,MPI_BYTE,0,DATA_TAG,MPI_COMM_WORLD); 
}

void  CORE_ENGINE_ACCELERATOR::
GetCoarseGrid() const {
    ARRAY<unsigned char> bytes;
    VECTOR< ARRAY<T, T_INDEX>, d> coarse_grid;
    engineInterface.GetCoarseGrid( coarse_grid );
    Marshal< VECTOR< ARRAY<T, T_INDEX>, d> >( coarse_grid, bytes );
    Write_Data( bytes );
}


void  CORE_ENGINE_ACCELERATOR::
GetCellType() const {
    ARRAY<unsigned char> bytes;
    ARRAY<CELL_TYPE, T_INDEX> cell_type;
    engineInterface.GetCellType( cell_type );
    Marshal< ARRAY<CELL_TYPE, T_INDEX> >( cell_type, bytes );
    Write_Data( bytes );
}

void  CORE_ENGINE_ACCELERATOR::
SetCellType() {
    ARRAY<unsigned char> bytes;
    Read_Data( bytes );
    ARRAY<CELL_TYPE, T_INDEX> cell_type;
    UnMarshal< ARRAY<CELL_TYPE, T_INDEX> >( cell_type, bytes );
    engineInterface.SetCellType(cell_type);   
}

void  CORE_ENGINE_ACCELERATOR::
GetCellTypeMesh() const {
    ARRAY<unsigned char> bytes;
    ARRAY<CELL_TYPE, int> cell_type_mesh;
    engineInterface.GetCellTypeMesh( cell_type_mesh );
    Marshal< ARRAY<CELL_TYPE, int> >( cell_type_mesh, bytes );
    Write_Data( bytes );
}

void  CORE_ENGINE_ACCELERATOR::
SetCellTypeMesh() {
    ARRAY<unsigned char> bytes;
    Read_Data( bytes );
    ARRAY<CELL_TYPE, int> cell_type_mesh;
    UnMarshal< ARRAY<CELL_TYPE, int> >( cell_type_mesh, bytes );
    engineInterface.SetCellTypeMesh( cell_type_mesh );
}

void  CORE_ENGINE_ACCELERATOR::
GetCellIndicesMesh() const {
    ARRAY<unsigned char> bytes;
    ARRAY<T_INDEX, int> cell_indices_mesh;
    engineInterface.GetCellIndicesMesh( cell_indices_mesh );
    Marshal< ARRAY<T_INDEX, int> >( cell_indices_mesh, bytes );
    Write_Data( bytes );
}

void  CORE_ENGINE_ACCELERATOR::
SetCellIndicesMesh() {
    ARRAY<unsigned char> bytes;
    ARRAY<T_INDEX, int> cell_indices_mesh;
    Read_Data( bytes );
    UnMarshal< ARRAY<T_INDEX, int> >( cell_indices_mesh, bytes );
    engineInterface.SetCellIndicesMesh( cell_indices_mesh );
}

void  CORE_ENGINE_ACCELERATOR::
GetCellsMesh() const {
    ARRAY<unsigned char> bytes;
    ARRAY<VECTOR<int,8>, int> cells_mesh;
    engineInterface.GetCellsMesh(cells_mesh);
    Marshal< ARRAY<VECTOR<int,8>, int> >( cells_mesh, bytes );
    Write_Data( bytes ); 
}

void  CORE_ENGINE_ACCELERATOR::
SetCellsMesh() {
    ARRAY<unsigned char> bytes;
    ARRAY<VECTOR<int,8>, int> cells_mesh;
    Read_Data( bytes );
    UnMarshal< ARRAY<VECTOR<int,8>, int> >( cells_mesh, bytes );
    engineInterface.SetCellsMesh(cells_mesh);
}

void  CORE_ENGINE_ACCELERATOR::
Number_Of_Mesh_Nodes() const {
    int result = engineInterface.Number_Of_Mesh_Nodes();
    // Send data back
    MPI_Send(&(result),1,MPI_INT,0,DATA_TAG,MPI_COMM_WORLD);
}

void  CORE_ENGINE_ACCELERATOR::
Number_Of_Mesh_Cells() const {
    int result = engineInterface.Number_Of_Mesh_Cells();
    // Send data back
    MPI_Send(&(result),1,MPI_INT,0,DATA_TAG,MPI_COMM_WORLD);
}

void CORE_ENGINE_ACCELERATOR::
CreateNewConstraint() {
    ENGINE_INTERFACE::CONSTRAINT_TYPE ctype;
    int cid;
    MPI_Status status; 
    MPI_Recv(&(ctype),1,MPI_INT,0,DATA_TAG,MPI_COMM_WORLD,&status);
    cid = engineInterface.CreateNewConstraint(ctype);
    // Send data back
    MPI_Send(&(cid),1,MPI_INT,0,DATA_TAG,MPI_COMM_WORLD);
}

void  CORE_ENGINE_ACCELERATOR::
GetConstraint() const {
    ENGINE_INTERFACE::CONSTRAINT_TYPE ctype;
    int cid;
    MPI_Status status; 
    MPI_Recv(&(ctype),1,MPI_INT,0,DATA_TAG,MPI_COMM_WORLD,&status);
    MPI_Recv(&(cid),1,MPI_INT,0,DATA_TAG,MPI_COMM_WORLD,&status);

    CONSTRAINT_SEGMENT<T,d> cs;
    engineInterface.GetConstraint(ctype, cid, cs );
    // Send data back
    MPI_Send(&(cs),sizeof(CONSTRAINT_SEGMENT<T,d>),MPI_BYTE,0,DATA_TAG,MPI_COMM_WORLD);
}

void  CORE_ENGINE_ACCELERATOR::
SetConstraint() {
    ENGINE_INTERFACE::CONSTRAINT_TYPE ctype;
    int cid;
    CONSTRAINT_SEGMENT<T,d> cs;
    MPI_Status status; 
    MPI_Recv(&(ctype),1,MPI_INT,0,DATA_TAG,MPI_COMM_WORLD,&status);
    MPI_Recv(&(cid),1,MPI_INT,0,DATA_TAG,MPI_COMM_WORLD,&status);
    MPI_Recv(&(cs),sizeof(CONSTRAINT_SEGMENT<T,d>),MPI_BYTE,0,DATA_TAG,MPI_COMM_WORLD,&status);
    engineInterface.SetConstraint(ctype, cid, cs );
}

void  CORE_ENGINE_ACCELERATOR::
RemoveConstraint() {
    ENGINE_INTERFACE::CONSTRAINT_TYPE ctype;
    int cid;
    MPI_Status status; 
    MPI_Recv(&(ctype),1,MPI_INT,0,DATA_TAG,MPI_COMM_WORLD,&status);
    MPI_Recv(&(cid),1,MPI_INT,0,DATA_TAG,MPI_COMM_WORLD,&status);
    engineInterface.RemoveConstraint( ctype, cid );
    // Send data back
    MPI_Send(&(cid),1,MPI_INT,0,DATA_TAG,MPI_COMM_WORLD);
}

void  CORE_ENGINE_ACCELERATOR::
NumberOfConstraints() const {
    ENGINE_INTERFACE::CONSTRAINT_TYPE ctype;
    int result;
    MPI_Status status; 
    MPI_Recv(&(ctype),1,MPI_INT,0,DATA_TAG,MPI_COMM_WORLD,&status);
    result=engineInterface.NumberOfConstraints(ctype);
    // Send data back
    MPI_Send(&(result),1,MPI_INT,0,DATA_TAG,MPI_COMM_WORLD);
}

void CORE_ENGINE_ACCELERATOR::
GetCollisionConstants() const {
    ARRAY<unsigned char> bytes;
    ARRAY<T, int> collision_spring_constants;
    engineInterface.GetCollisionConstants(collision_spring_constants);
    Marshal< ARRAY<T, int> >( collision_spring_constants, bytes );
    Write_Data( bytes );
}

void  CORE_ENGINE_ACCELERATOR::
SetCollisionConstants() {
    ARRAY<unsigned char> bytes;
    ARRAY<T, int> collision_spring_constants;
    Read_Data( bytes );
    UnMarshal< ARRAY<T, int> >( collision_spring_constants, bytes );
    engineInterface.SetCollisionConstants(collision_spring_constants);
}

void  CORE_ENGINE_ACCELERATOR::
GetCollisionLocations() const {
    ARRAY<unsigned char> bytes;
    ARRAY<TV, int> collision_spring_locations;
    engineInterface.GetCollisionLocations(collision_spring_locations);
    Marshal< ARRAY<TV, int> >(collision_spring_locations, bytes );
    Write_Data( bytes );
}

void  CORE_ENGINE_ACCELERATOR::
SetCollisionLocations() {
    ARRAY<unsigned char> bytes;
    ARRAY<TV, int> collision_spring_locations;
    Read_Data( bytes );
    UnMarshal< ARRAY<TV, int> >( collision_spring_locations, bytes );
    engineInterface.SetCollisionLocations(collision_spring_locations );
}

void CORE_ENGINE_ACCELERATOR::
GetMuscleData() const {
    MPI_Status status; 
    ARRAY<unsigned char> bytes;
    int result;
 
    int num_muscles;
    ARRAY< ARRAY<int>, T_INDEX> cell_muscles;
    ARRAY< ARRAY<int>, int> cell_muscles_mesh;
    ARRAY< ARRAY<T>, T_INDEX> cell_densities;
    ARRAY< ARRAY<T>, int> cell_densities_mesh;
    ARRAY< ARRAY<TV>, T_INDEX> cell_fibers;
    ARRAY< ARRAY<TV>, int>  cell_fibers_mesh;

    engineInterface.GetMuscleData(num_muscles,
                                  cell_muscles,
                                  cell_muscles_mesh,
                                  cell_densities,
                                  cell_densities_mesh,
                                  cell_fibers,
                                  cell_fibers_mesh );
    
    MPI_Send(&(num_muscles),1,MPI_INT,0,DATA_TAG,MPI_COMM_WORLD);
    Marshal< ARRAY< ARRAY<int>, T_INDEX> >( cell_muscles, bytes );
    Write_Data( bytes );

    Marshal< ARRAY< ARRAY<int>, int> >( cell_muscles_mesh, bytes );
    Write_Data( bytes );

    Marshal< ARRAY< ARRAY<T>, T_INDEX> >( cell_densities, bytes );
    Write_Data( bytes );

    Marshal< ARRAY< ARRAY<T>, int> >( cell_densities_mesh, bytes );
    Write_Data( bytes );

    Marshal< ARRAY< ARRAY<TV>, T_INDEX> >( cell_fibers, bytes );
    Write_Data( bytes );
    
    Marshal< ARRAY< ARRAY<TV>, int> >( cell_fibers_mesh, bytes );
    Write_Data( bytes );

}

void CORE_ENGINE_ACCELERATOR::
SetMuscleData() {

    MPI_Status status; 
    ARRAY<unsigned char> bytes;

    int num_muscles;
    ARRAY< ARRAY<int>, T_INDEX> cell_muscles;
    ARRAY< ARRAY<int>, int> cell_muscles_mesh;
    ARRAY< ARRAY<T>, T_INDEX> cell_densities;
    ARRAY< ARRAY<T>, int> cell_densities_mesh;
    ARRAY< ARRAY<TV>, T_INDEX> cell_fibers;
    ARRAY< ARRAY<TV>, int>  cell_fibers_mesh;

    MPI_Recv(&(num_muscles),1,MPI_INT,0,DATA_TAG,MPI_COMM_WORLD,&status);
    Read_Data( bytes );
    UnMarshal< ARRAY< ARRAY<int>, T_INDEX> >( cell_muscles, bytes );

    Read_Data( bytes );
    UnMarshal< ARRAY< ARRAY<int>, int> >( cell_muscles_mesh, bytes );

    Read_Data( bytes );
    UnMarshal< ARRAY< ARRAY<T>, T_INDEX> >( cell_densities, bytes );

    Read_Data( bytes );
    UnMarshal< ARRAY< ARRAY<T>, int> >( cell_densities_mesh, bytes );

    Read_Data( bytes );
    UnMarshal< ARRAY< ARRAY<TV>, T_INDEX> >( cell_fibers, bytes );

    Read_Data( bytes );
    UnMarshal< ARRAY< ARRAY<TV>, int> >( cell_fibers_mesh, bytes );

    engineInterface.SetMuscleData(num_muscles,
                                  cell_muscles,
                                  cell_muscles_mesh,
                                  cell_densities,
                                  cell_densities_mesh,
                                  cell_fibers,
                                  cell_fibers_mesh );
}

void CORE_ENGINE_ACCELERATOR::
GetMuscleActivations() const {
    ARRAY<unsigned char> bytes;
    ARRAY<T, int> muscle_activations;
    engineInterface.GetMuscleActivations( muscle_activations );
    Marshal< ARRAY<T, int> >(muscle_activations, bytes );
    Write_Data( bytes );
}

void CORE_ENGINE_ACCELERATOR::
SetMuscleActivations(){
    ARRAY<unsigned char> bytes;
    ARRAY<T, int> muscle_activations;
    Read_Data( bytes );
    UnMarshal< ARRAY<T, int> >(muscle_activations, bytes );
    engineInterface.SetMuscleActivations( muscle_activations );
}

void CORE_ENGINE_ACCELERATOR::
GetMuscleMaxStress() const {
    ARRAY<unsigned char> bytes;
    ARRAY<T, int> fiber_max_stress;
    engineInterface.GetMuscleFiberMaxStress(fiber_max_stress);
    Marshal< ARRAY<T, int> >( fiber_max_stress, bytes );
    Write_Data( bytes );
}

void CORE_ENGINE_ACCELERATOR::
SetMuscleMaxStress() {
    ARRAY<unsigned char> bytes;
    ARRAY<T, int> fiber_max_stress;
    Read_Data( bytes );
    UnMarshal< ARRAY<T, int> >( fiber_max_stress, bytes );
    engineInterface.SetMuscleFiberMaxStress(fiber_max_stress);
}


void CORE_ENGINE_ACCELERATOR::
EngineReady() const {
    bool result = engineInterface.EngineReady();
    // Send data back
    MPI_Send(&(result),sizeof(bool),MPI_BYTE,0,DATA_TAG,MPI_COMM_WORLD);        
}
    
void CORE_ENGINE_ACCELERATOR::
EngineInitialized() {
    bool result = engineInterface.EngineInitialized();
    // Send data back
    MPI_Send(&(result),sizeof(bool),MPI_BYTE,0,DATA_TAG,MPI_COMM_WORLD);        
}


}

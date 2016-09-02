#ifdef ENABLE_MPI_OFFLOAD
#include <mpi.h>
#include <cstring>
#include <iostream>
#include <PhysBAM_Tools/Read_Write/Grids_Uniform_Arrays/READ_WRITE_ARRAYS.h>
#include <PhysBAM_Tools/Read_Write/Arrays/READ_WRITE_ARRAY.h>
#include <PhysBAM_Tools/Read_Write/Vectors/READ_WRITE_VECTOR_3D.h>
#include <PhysBAM_Tools/Read_Write/Data_Structures/READ_WRITE_TRIPLE.h>

#include <EngineInterfaceMPI/ENGINE_INTERFACE_MPI.h>
#include <EngineInterfaceMPI/COMMAND_ID.h>

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
        MPI_Recv(&(buffer_size),1,MPI_INT,1,DATA_TAG,MPI_COMM_WORLD,&status);
        MPI_Recv(&(m),1,MPI_INT,1,DATA_TAG,MPI_COMM_WORLD,&status);
        bytes.Exact_Resize(buffer_size, false);
        MPI_Recv(bytes.base_pointer,buffer_size,MPI_BYTE,1,DATA_TAG,MPI_COMM_WORLD,&status);
        bytes.m = m;
    }
    
    void Write_Data( ARRAY<unsigned char>& bytes ){
        MPI_Send(&(bytes.buffer_size),1,MPI_INT,1,DATA_TAG,MPI_COMM_WORLD);
        MPI_Send(&(bytes.m),1,MPI_INT,1,DATA_TAG,MPI_COMM_WORLD);
        MPI_Send(bytes.base_pointer,bytes.buffer_size,MPI_BYTE,1,DATA_TAG,MPI_COMM_WORLD);
    }
    
    void Start_Command( MPI_COMMAND_ID command ){
        MPI_Send(&command, 1, MPI_INT, 1, COMMAND_TAG, MPI_COMM_WORLD);
    }
    
#if 0
#define START_COMMAND(C) LOG::SCOPE scope(std::string("Executing Method ") + std::string(MPI_METHOD_NAMES[int(C)]) ); Start_Command(C);
    #else
#define START_COMMAND(C) Start_Command(C);
#endif
#define END_COMMAND() MPI_Barrier( MPI_COMM_WORLD );
        

    ENGINE_INTERFACE_MPI::ENGINE_INTERFACE_MPI()    {
        int flag;
        MPI_Initialized(&flag);
        if( !flag ){
            LOG::cout << "CORE ENGINE FAILURE! No MPI Comm is available. " << std::endl;
            exit(1);
        }
        else
            LOG::cout << "MPI Initialized and Ready." << std::endl;
    }

    ENGINE_INTERFACE_MPI::~ENGINE_INTERFACE_MPI()
    {
    }


    void ENGINE_INTERFACE_MPI::
CreateEngine(T_INDEX bounds, T dx, TV min_corner, T mu, T lambda, T alpha, T cutoff_value, T stabilization){
    START_COMMAND(CREATE_ENGINE);

    MPI_Send(&(bounds.x),3,MPI_INT,1,DATA_TAG,MPI_COMM_WORLD);
    MPI_Send(&(dx),1,MPI_FLOAT,1,DATA_TAG,MPI_COMM_WORLD);
    MPI_Send(&(min_corner.x),3,MPI_FLOAT,1,DATA_TAG,MPI_COMM_WORLD);
    MPI_Send(&(mu),1,MPI_FLOAT,1,DATA_TAG,MPI_COMM_WORLD);
    MPI_Send(&(lambda),1,MPI_FLOAT,1,DATA_TAG,MPI_COMM_WORLD);
    MPI_Send(&(alpha),1,MPI_FLOAT,1,DATA_TAG,MPI_COMM_WORLD);
    MPI_Send(&(cutoff_value),1,MPI_FLOAT,1,DATA_TAG,MPI_COMM_WORLD);
    MPI_Send(&(stabilization),1,MPI_FLOAT,1,DATA_TAG,MPI_COMM_WORLD);

    END_COMMAND();
}

void  ENGINE_INTERFACE_MPI::
Initialize_Muscles(){
    START_COMMAND(INITIALIZE_MUSCLES);
    END_COMMAND();
}

void  ENGINE_INTERFACE_MPI::
Initialize_Mesh(int mesh_cell_count, int mesh_node_count){
    START_COMMAND(INITIALIZE_MESH);
    MPI_Send(&(mesh_cell_count),1,MPI_INT,1,DATA_TAG,MPI_COMM_WORLD);
    MPI_Send(&(mesh_node_count),1,MPI_INT,1,DATA_TAG,MPI_COMM_WORLD);

    END_COMMAND();
}

void  ENGINE_INTERFACE_MPI::
InitializeEngine(){
    START_COMMAND(INITIALIZE_ENGINE);
    END_COMMAND();
}

void  ENGINE_INTERFACE_MPI::
Exact_Solve(int krylov_iterations, int newton_iterations, T krylov_tolerance, T newton_tolerance, bool no_cut_cells, float& _result){
    START_COMMAND(EXACT_SOLVE);
    MPI_Send(&(krylov_iterations),1,MPI_INT,1,DATA_TAG,MPI_COMM_WORLD);
    MPI_Send(&(newton_iterations),1,MPI_INT,1,DATA_TAG,MPI_COMM_WORLD);
    MPI_Send(&(krylov_tolerance),1,MPI_FLOAT,1,DATA_TAG,MPI_COMM_WORLD);
    MPI_Send(&(newton_tolerance),1,MPI_FLOAT,1,DATA_TAG,MPI_COMM_WORLD);
    MPI_Send(&(no_cut_cells),1,MPI_C_BOOL,1,DATA_TAG,MPI_COMM_WORLD);

    //TODO - BLOCK HERE FOR EXACT SOLVE TO FINISH
    T result;
    MPI_Status status; 
    MPI_Recv(&(result),1,MPI_FLOAT,1,DATA_TAG,MPI_COMM_WORLD,&status);
    _result = result;
    END_COMMAND();
}

void  ENGINE_INTERFACE_MPI::
DestroyEngine(){
    START_COMMAND(DESTROY_ENGINE);
    END_COMMAND();
}

ENGINE_INTERFACE_MPI::TV  ENGINE_INTERFACE_MPI::
Displacement_Grid( T_INDEX cell, TV weights ) const {
    START_COMMAND(DISPLACEMENT_GRID);
    MPI_Send(&(cell.x),3,MPI_INT,1,DATA_TAG,MPI_COMM_WORLD);
    MPI_Send(&(weights.x),3,MPI_FLOAT,1,DATA_TAG,MPI_COMM_WORLD);

    TV result;
    MPI_Status status; 
    MPI_Recv(&(result.x),3,MPI_FLOAT,1,DATA_TAG,MPI_COMM_WORLD,&status);
    END_COMMAND();
    return result;
}

ENGINE_INTERFACE_MPI::TV  ENGINE_INTERFACE_MPI::
Displacement_Mesh( int cell, TV weights ) const {
    START_COMMAND(DISPLACEMENT_MESH);
    MPI_Send(&(cell),1,MPI_INT,1,DATA_TAG,MPI_COMM_WORLD);
    MPI_Send(&(weights.x),3,MPI_FLOAT,1,DATA_TAG,MPI_COMM_WORLD);

    TV result;
    MPI_Status status; 
    MPI_Recv(&(result.x),3,MPI_FLOAT,1,DATA_TAG,MPI_COMM_WORLD,&status);
    END_COMMAND();
    return result;
}

void  ENGINE_INTERFACE_MPI::
Displacement( ARRAY<TRIPLE<int, T_INDEX, TV> >& queue, ARRAY<TV>& updates) const {
    START_COMMAND(DISPLACEMENT);
    ARRAY<unsigned char> bytes;
    Marshal< ARRAY<TRIPLE<int, T_INDEX, TV> > >( queue, bytes );
    Write_Data( bytes );

    Read_Data(bytes);
    UnMarshal< ARRAY<TV> >( updates, bytes );    
    END_COMMAND();
}

void  ENGINE_INTERFACE_MPI::
Deformation( ARRAY<TRIPLE<int, T_INDEX, TV> >& queue, ARRAY<TV>& updates) const {
    START_COMMAND(DEFORMATION);
    ARRAY<unsigned char> bytes;
    Marshal< ARRAY<TRIPLE<int, T_INDEX, TV> > >( queue, bytes );
    Write_Data( bytes );

    Read_Data(bytes);
    UnMarshal< ARRAY<TV> >( updates, bytes );   
    END_COMMAND();
}

void  ENGINE_INTERFACE_MPI::
Stress( ARRAY<TRIPLE<int, T_INDEX, TV> >& queue, ARRAY<TV>& updates) const {
    START_COMMAND(STRESS);
    ARRAY<unsigned char> bytes;
    Marshal< ARRAY<TRIPLE<int, T_INDEX, TV> > >( queue, bytes );
    Write_Data( bytes );

    Read_Data(bytes);
    UnMarshal< ARRAY<TV> >( updates, bytes );   
    END_COMMAND();
}

void  ENGINE_INTERFACE_MPI::
Strain( ARRAY<TRIPLE<int, T_INDEX, TV> >& queue, ARRAY<TV>& updates) const {
    START_COMMAND(STRAIN);
    ARRAY<unsigned char> bytes;
    Marshal< ARRAY<TRIPLE<int, T_INDEX, TV> > >( queue, bytes );
    Write_Data( bytes );

    Read_Data(bytes);
    UnMarshal< ARRAY<TV> >( updates, bytes );   
    END_COMMAND();
}

ENGINE_INTERFACE_MPI::T   ENGINE_INTERFACE_MPI::
h() const {
    START_COMMAND(H);
    T result;
    MPI_Status status; 
    MPI_Recv(&(result),1,MPI_FLOAT,1,DATA_TAG,MPI_COMM_WORLD,&status);
    END_COMMAND();
    return result;
}

ENGINE_INTERFACE_MPI::T_INDEX  ENGINE_INTERFACE_MPI::
Cell( TV& location) const {
    START_COMMAND(CELL);
    MPI_Send(&(location.x),3,MPI_FLOAT,1,DATA_TAG,MPI_COMM_WORLD);
    
    T_INDEX result;
    MPI_Status status; 
    MPI_Recv(&(result.x),3,MPI_INT,1,DATA_TAG,MPI_COMM_WORLD,&status);
    END_COMMAND();
    return result;
}

ENGINE_INTERFACE_MPI::TV  ENGINE_INTERFACE_MPI::
Node( T_INDEX& location) const {
    START_COMMAND(NODE);
    MPI_Send(&(location.x),3,MPI_INT,1,DATA_TAG,MPI_COMM_WORLD);
    
    TV result;
    MPI_Status status; 
    MPI_Recv(&(result.x),3,MPI_FLOAT,1,DATA_TAG,MPI_COMM_WORLD,&status);
    END_COMMAND();
    return result;  
}

RANGE<ENGINE_INTERFACE_MPI::T_INDEX>  ENGINE_INTERFACE_MPI::
Padded_Node_Domain() const {
    START_COMMAND(PADDED_NODE_DOMAIN);
    RANGE<T_INDEX> result;
    MPI_Status status; 
    MPI_Recv(&(result.min_corner.x),6,MPI_INT,1,DATA_TAG,MPI_COMM_WORLD,&status);
    END_COMMAND();
    return result;  
}

RANGE<ENGINE_INTERFACE_MPI::T_INDEX>  ENGINE_INTERFACE_MPI::
Unpadded_Node_Domain() const {
    START_COMMAND(UNPADDED_NODE_DOMAIN);
    RANGE<T_INDEX> result;
    MPI_Status status; 
    MPI_Recv(&(result.min_corner.x),6,MPI_INT,1,DATA_TAG,MPI_COMM_WORLD,&status);
    END_COMMAND();
    return result;     
}

RANGE<ENGINE_INTERFACE_MPI::T_INDEX>  ENGINE_INTERFACE_MPI::
Padded_Cell_Domain() const {
    START_COMMAND(PADDED_CELL_DOMAIN);
    RANGE<T_INDEX> result;
    MPI_Status status; 
    MPI_Recv(&(result.min_corner.x),6,MPI_INT,1,DATA_TAG,MPI_COMM_WORLD,&status);
    END_COMMAND();
    return result;  
}

RANGE<ENGINE_INTERFACE_MPI::T_INDEX>  ENGINE_INTERFACE_MPI::
Unpadded_Cell_Domain() const {
    START_COMMAND(UNPADDED_CELL_DOMAIN);
    RANGE<T_INDEX> result;
    MPI_Status status; 
    MPI_Recv(&(result.min_corner.x),6,MPI_INT,1,DATA_TAG,MPI_COMM_WORLD,&status);
    END_COMMAND();
    return result;  
}

GRID<ENGINE_INTERFACE_MPI::TV>  ENGINE_INTERFACE_MPI::
Grid() const {
    START_COMMAND(GET_GRID);
    GRID<TV> result;
    MPI_Status status; 
    MPI_Recv(&(result.counts.x),80,MPI_BYTE,1,DATA_TAG,MPI_COMM_WORLD,&status);
    END_COMMAND();
    return result;     
}

void ENGINE_INTERFACE_MPI::
GetCoarseGrid( VECTOR< ARRAY< T, T_INDEX >, d>& grid_u ) const {
    START_COMMAND(GET_COARSEGRID);
    MPI_Status status; 
    ARRAY<unsigned char> bytes;
    Read_Data( bytes );
    UnMarshal< VECTOR< ARRAY< T, T_INDEX >, d> >( grid_u, bytes );
    END_COMMAND();
}


void  ENGINE_INTERFACE_MPI::
GetCellType(ARRAY<CELL_TYPE, T_INDEX>& cell_type) const {
    START_COMMAND(GET_CELL_TYPE);
    ARRAY<unsigned char> bytes;
    Read_Data( bytes );
    UnMarshal< ARRAY<CELL_TYPE, T_INDEX> >( cell_type, bytes );
    END_COMMAND();
}

void  ENGINE_INTERFACE_MPI::
SetCellType( ARRAY<CELL_TYPE, T_INDEX>& cell_type) {
    START_COMMAND(SET_CELL_TYPE);
    ARRAY<unsigned char> bytes;
    Marshal< ARRAY<CELL_TYPE, T_INDEX> >( cell_type, bytes );
    Write_Data( bytes );
    END_COMMAND();
}

void  ENGINE_INTERFACE_MPI::
GetCellTypeMesh(ARRAY<CELL_TYPE, int>& cell_type_mesh) const {
    START_COMMAND(GET_CELL_TYPE_MESH);
    ARRAY<unsigned char> bytes;
    Read_Data( bytes );
    UnMarshal< ARRAY<CELL_TYPE, int> >( cell_type_mesh, bytes );
    END_COMMAND();
}

void  ENGINE_INTERFACE_MPI::
SetCellTypeMesh( ARRAY<CELL_TYPE, int>& cell_type_mesh) {
    START_COMMAND(SET_CELL_TYPE_MESH);
    ARRAY<unsigned char> bytes;
    Marshal< ARRAY<CELL_TYPE, int> >( cell_type_mesh, bytes );
    Write_Data( bytes );
    END_COMMAND();
}

void  ENGINE_INTERFACE_MPI::
GetCellIndicesMesh(ARRAY<T_INDEX,  int>& cell_indices_mesh) const {
    START_COMMAND(GET_CELL_INDICES_MESH);
    ARRAY<unsigned char> bytes;
    Read_Data( bytes );
    UnMarshal< ARRAY<T_INDEX, int> >(cell_indices_mesh, bytes );
    END_COMMAND();
}

void  ENGINE_INTERFACE_MPI::
SetCellIndicesMesh( ARRAY<T_INDEX, int>& cell_indices_mesh) {
    START_COMMAND(SET_CELL_INDICES_MESH);
    ARRAY<unsigned char> bytes;
    Marshal< ARRAY<T_INDEX, int> >( cell_indices_mesh, bytes );
    Write_Data( bytes );
    END_COMMAND();
}

void  ENGINE_INTERFACE_MPI::
GetCellsMesh(ARRAY<VECTOR<int, 8>, int>& cells_mesh) const {
    START_COMMAND(GET_CELLS_MESH);
    ARRAY<unsigned char> bytes;
    Read_Data( bytes );
    UnMarshal< ARRAY<VECTOR<int, 8>, int> >(cells_mesh, bytes );
    END_COMMAND();
}

void  ENGINE_INTERFACE_MPI::
SetCellsMesh( ARRAY<VECTOR<int, 8>, int>& cells_mesh) {
    START_COMMAND(SET_CELLS_MESH);
    ARRAY<unsigned char> bytes;
    Marshal< ARRAY<VECTOR<int, 8>, int> >( cells_mesh, bytes );
    Write_Data( bytes );
    END_COMMAND();
}

int  ENGINE_INTERFACE_MPI::
Number_Of_Mesh_Nodes() const {
    START_COMMAND(NUMBER_OF_MESH_NODES);
    int result;
    MPI_Status status; 
    MPI_Recv(&(result),1,MPI_INT,1,DATA_TAG,MPI_COMM_WORLD,&status);
    END_COMMAND();
    return result;     
}

int  ENGINE_INTERFACE_MPI::
Number_Of_Mesh_Cells() const {
    START_COMMAND(NUMBER_OF_MESH_CELLS);
    int result;
    MPI_Status status; 
    MPI_Recv(&(result),1,MPI_INT,1,DATA_TAG,MPI_COMM_WORLD,&status);
    END_COMMAND();
    return result;     
}

void ENGINE_INTERFACE_MPI::
Output_Structures(GEOMETRY_PARTICLES<TV>& particles,ARRAY<STRUCTURE<TV>*>& collection_structures) const{
}

int  ENGINE_INTERFACE_MPI::
CreateNewConstraint(CONSTRAINT_TYPE ctype) {
    START_COMMAND(CREATE_NEW_CONSTRAINT);
    MPI_Send(&(ctype),1,MPI_INT,1,DATA_TAG,MPI_COMM_WORLD);
    int result;
    MPI_Status status; 
    MPI_Recv(&(result),1,MPI_INT,1,DATA_TAG,MPI_COMM_WORLD,&status);
    END_COMMAND();
    return result;   
}

void  ENGINE_INTERFACE_MPI::
GetConstraint(CONSTRAINT_TYPE ctype, int cid, CONSTRAINT_SEGMENT<T,d>& cs) const {
    START_COMMAND(GET_CONSTRAINT);
    MPI_Send(&(ctype),1,MPI_INT,1,DATA_TAG,MPI_COMM_WORLD);
    MPI_Send(&(cid),1,MPI_INT,1,DATA_TAG,MPI_COMM_WORLD);
    
    MPI_Status status; 
    MPI_Recv(&(cs),sizeof(CONSTRAINT_SEGMENT<T,d>),MPI_BYTE,1,DATA_TAG,MPI_COMM_WORLD,&status);
    END_COMMAND();
}

void  ENGINE_INTERFACE_MPI::
SetConstraint(CONSTRAINT_TYPE ctype, int cid,  CONSTRAINT_SEGMENT<T,d>& cs) {
    START_COMMAND(SET_CONSTRAINT);
    MPI_Send(&(ctype),1,MPI_INT,1,DATA_TAG,MPI_COMM_WORLD);
    MPI_Send(&(cid),1,MPI_INT,1,DATA_TAG,MPI_COMM_WORLD);
    MPI_Send(&(cs),sizeof(CONSTRAINT_SEGMENT<T,d>),MPI_BYTE,1,DATA_TAG,MPI_COMM_WORLD);
    END_COMMAND();
}

void  ENGINE_INTERFACE_MPI::
RemoveConstraint(CONSTRAINT_TYPE ctype, int& cid) {
    START_COMMAND(REMOVE_CONSTRAINT);
    MPI_Send(&(ctype),1,MPI_INT,1,DATA_TAG,MPI_COMM_WORLD);
    MPI_Send(&(cid),1,MPI_INT,1,DATA_TAG,MPI_COMM_WORLD);
    MPI_Status status; 
    MPI_Recv(&(cid),1,MPI_INT,1,DATA_TAG,MPI_COMM_WORLD,&status);
    END_COMMAND();
}

int  ENGINE_INTERFACE_MPI::
NumberOfConstraints(CONSTRAINT_TYPE ctype) const {
    START_COMMAND(NUMBER_OF_CONSTRAINTS);
    MPI_Send(&(ctype),1,MPI_INT,1,DATA_TAG,MPI_COMM_WORLD);

    int result;
    MPI_Status status; 
    MPI_Recv(&(result),1,MPI_INT,1,DATA_TAG,MPI_COMM_WORLD,&status);
    END_COMMAND();
    return result;   
}

void  ENGINE_INTERFACE_MPI::
GetCollisionConstants(ARRAY<T>& collision_constants) const {
    START_COMMAND(GET_COLLISION_CONSTANTS);
    ARRAY<unsigned char> bytes;
    Read_Data( bytes );
    UnMarshal< ARRAY<T, int> >(collision_constants, bytes );
    END_COMMAND();
}

void  ENGINE_INTERFACE_MPI::
SetCollisionConstants( ARRAY<T>& collision_constants) {
    START_COMMAND(SET_COLLISION_CONSTANTS);
    ARRAY<unsigned char> bytes;
    Marshal< ARRAY<T, int> >( collision_constants, bytes );
    Write_Data( bytes );
    END_COMMAND();
}

void  ENGINE_INTERFACE_MPI::
GetCollisionLocations(ARRAY<TV>& collision_locations) const {
    START_COMMAND(GET_COLLISION_LOCATIONS);
    ARRAY<unsigned char> bytes;
    Read_Data( bytes );
    UnMarshal< ARRAY<TV, int> >(collision_locations, bytes );
    END_COMMAND();
}

void  ENGINE_INTERFACE_MPI::
SetCollisionLocations( ARRAY<TV>& collision_locations) {
    START_COMMAND(SET_COLLISION_LOCATIONS);
    ARRAY<unsigned char> bytes;
    Marshal< ARRAY<TV, int> >( collision_locations, bytes );
    Write_Data( bytes );
    END_COMMAND();
}

void ENGINE_INTERFACE_MPI::
GetMuscleData(int& max_muscle,
              ARRAY< ARRAY<int>, T_INDEX>& muscle_ids,
              ARRAY< ARRAY<int>, int>& muscle_ids_mesh, 
              ARRAY< ARRAY<T>, T_INDEX>& muscle_density,
              ARRAY< ARRAY<T>, int>& muscle_density_mesh, 
              ARRAY< ARRAY<TV>, T_INDEX>& muscle_fiber,
              ARRAY< ARRAY<TV>, int>& muscle_fiber_mesh) const {
        START_COMMAND(GET_MUSCLE_DATA);
        MPI_Status status; 
        ARRAY<unsigned char> bytes;

        MPI_Recv(&(max_muscle),1,MPI_INT,1,DATA_TAG,MPI_COMM_WORLD,&status);
        Read_Data( bytes );
        UnMarshal< ARRAY< ARRAY<int>, T_INDEX> >(muscle_ids, bytes );
        Read_Data( bytes );
        UnMarshal< ARRAY< ARRAY<int>, int> >(muscle_ids_mesh, bytes );
        Read_Data( bytes );
        UnMarshal< ARRAY< ARRAY<T>, T_INDEX> >(muscle_density, bytes );
        Read_Data( bytes );
        UnMarshal< ARRAY< ARRAY<T>, int> >(muscle_density_mesh, bytes );
        Read_Data( bytes );
        UnMarshal< ARRAY< ARRAY<TV>, T_INDEX> >(muscle_fiber, bytes );
        Read_Data( bytes );
        UnMarshal< ARRAY< ARRAY<TV>, int> >(muscle_fiber_mesh, bytes );

        END_COMMAND();
}

void ENGINE_INTERFACE_MPI::
SetMuscleData(int& max_muscle,
              ARRAY< ARRAY<int>, T_INDEX>& muscle_ids,
              ARRAY< ARRAY<int>, int>& muscle_ids_mesh, 
              ARRAY< ARRAY<T>, T_INDEX>& muscle_density,
              ARRAY< ARRAY<T>, int>& muscle_density_mesh, 
              ARRAY< ARRAY<TV>, T_INDEX>& muscle_fiber,
              ARRAY< ARRAY<TV>, int>& muscle_fiber_mesh){
        START_COMMAND(SET_MUSCLE_DATA);
        ARRAY<unsigned char> bytes;

        MPI_Send(&(max_muscle),1,MPI_INT,1,DATA_TAG,MPI_COMM_WORLD);
        Marshal< ARRAY< ARRAY<int>, T_INDEX> >(muscle_ids, bytes );
        Write_Data( bytes );
        Marshal< ARRAY< ARRAY<int>, int> >(muscle_ids_mesh, bytes );
        Write_Data( bytes );
        Marshal< ARRAY< ARRAY<T>, T_INDEX> >(muscle_density, bytes );
        Write_Data( bytes );
        Marshal< ARRAY< ARRAY<T>, int> >(muscle_density_mesh, bytes );
        Write_Data( bytes );
        Marshal< ARRAY< ARRAY<TV>, T_INDEX> >(muscle_fiber, bytes );
        Write_Data( bytes );
        Marshal< ARRAY< ARRAY<TV>, int> >(muscle_fiber_mesh, bytes );
        Write_Data( bytes );

        END_COMMAND();
}
void ENGINE_INTERFACE_MPI::
GetMuscleActivations( ARRAY<T, int>& muscle_activations ) const {
        START_COMMAND(GET_MUSCLE_ACTIVATIONS);
        ARRAY<unsigned char> bytes;
        Read_Data( bytes );
        UnMarshal< ARRAY<T, int> >(muscle_activations, bytes );
        END_COMMAND();
}
void ENGINE_INTERFACE_MPI::
SetMuscleActivations( ARRAY<T, int>& muscle_activations ){
        START_COMMAND(SET_MUSCLE_ACTIVATIONS);
        ARRAY<unsigned char> bytes;
        Marshal< ARRAY<T, int> >(muscle_activations, bytes );
        Write_Data( bytes );
        END_COMMAND();
}

void ENGINE_INTERFACE_MPI::
GetMuscleFiberMaxStress( ARRAY<T, int>& fiber_max_stress ) const {
        START_COMMAND(GET_MUSCLE_MAX_STRESS);
        ARRAY<unsigned char> bytes;
        Read_Data( bytes );
        UnMarshal< ARRAY<T, int> >(fiber_max_stress, bytes );
        END_COMMAND();
}

void ENGINE_INTERFACE_MPI::
SetMuscleFiberMaxStress( ARRAY<T, int>& fiber_max_stress ){
        START_COMMAND(SET_MUSCLE_MAX_STRESS);
        ARRAY<unsigned char> bytes;
        Marshal< ARRAY<T, int> >(fiber_max_stress, bytes );
        Write_Data( bytes );
        END_COMMAND();
}
    

bool ENGINE_INTERFACE_MPI::
EngineReady() const {
    START_COMMAND(ENGINE_READY);
    bool result;
    MPI_Status status; 
    MPI_Recv(&(result),sizeof(bool),MPI_BYTE,1,DATA_TAG,MPI_COMM_WORLD,&status);
    END_COMMAND();
    return result;  
}
    
bool ENGINE_INTERFACE_MPI::
EngineInitialized() const {
    START_COMMAND(ENGINE_INITIALIZED);
    bool result;
    MPI_Status status; 
    MPI_Recv(&(result),sizeof(bool),MPI_BYTE,1,DATA_TAG,MPI_COMM_WORLD,&status);
    END_COMMAND();
    return result;
}


}
#endif

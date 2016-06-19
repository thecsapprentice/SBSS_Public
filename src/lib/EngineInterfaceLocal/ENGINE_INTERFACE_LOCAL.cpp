#include <PhysBAM_Tools/Krylov_Solvers/SYMMQMR.h>
#include <PhysBAM_Tools/Krylov_Solvers/CONJUGATE_GRADIENT.h>
#include <PhysBAM_Tools/Data_Structures/HASHTABLE.h>

#include <Common/ALIGNED_ARRAY.h>

#include <EngineBackend/NONLINEAR_ELASTICITY.h>
#include <EngineBackend/SKINNING_NONLINEAR_ELASTICITY.h>
#include <EngineBackend/HYBRID_NONLINEAR_ELASTICITY.h>

#include <EngineBackend/CG_SYSTEM.h>
#include <EngineBackend/CG_VECTOR.h>

#include <EngineBackend/ELASTIC_SOLVER.h>
#include <EngineBackend/Write_Output.h>

#include "ENGINE_INTERFACE_LOCAL.h"


namespace PhysBAM{

    ENGINE_INTERFACE_LOCAL::ENGINE_INTERFACE_LOCAL() :
        engine( NULL ), solver(NULL), engineCreated( false ), engineInitialized( false ) 
    {
    }

    ENGINE_INTERFACE_LOCAL::~ENGINE_INTERFACE_LOCAL()
    {
    }


void ENGINE_INTERFACE_LOCAL::
CreateEngine(T_INDEX bounds, T dx, TV min_corner, T mu, T lambda, T alpha, T cutoff_value, T stabilization){
    
    if( ! engineCreated ){
        engine = new T_DISCRETIZATION(bounds, dx, min_corner);
        engine->Initialize_Parameters(mu, lambda, alpha, cutoff_value, stabilization);              
        engineCreated = true;
    }
}

void  ENGINE_INTERFACE_LOCAL::
Initialize_Muscles(){
    if( engineCreated ){
        engine->Initialize_Muscles();
        int num_of_muscles = 0;
        engine->muscle_fiber_max_stresses.Resize(num_of_muscles);
        engine->muscle_activations.Resize(num_of_muscles);
    }
}

void  ENGINE_INTERFACE_LOCAL::
Initialize_Mesh(int mesh_cell_count, int mesh_node_count){
    if( engineCreated ){
        engine->Initialize_Mesh(mesh_cell_count, mesh_node_count);
    }
}

void  ENGINE_INTERFACE_LOCAL::
InitializeEngine(){
    if( engineCreated ){
        solver = new ELASTIC_SOLVER<T_STATE_BINDER, T_STATE, T_DISCRETIZATION>(engine, false);
        solver->Initialize(16);
        engineInitialized = true;
    }
}

void  ENGINE_INTERFACE_LOCAL::
Exact_Solve(const int krylov_iterations,const int newton_iterations,const T krylov_tolerance,const T newton_tolerance,const bool no_cut_cells, float& result){
    if( engineInitialized )
        result = solver->Exact_Solve(krylov_iterations, newton_iterations, krylov_tolerance,
                                     newton_tolerance, no_cut_cells);
}

void  ENGINE_INTERFACE_LOCAL::
DestroyEngine(){
    if( engineCreated ){
        if( engineInitialized ){
            delete solver;
            solver = NULL;
            engineInitialized = false;
        }
        delete engine;
        engine = NULL;
        engineCreated = false;
    }
}

ENGINE_INTERFACE_LOCAL::TV  ENGINE_INTERFACE_LOCAL::
Displacement_Grid( T_INDEX cell, TV weights ) const {
    if( engineInitialized ){
        return engine->Displacement_Grid( cell, weights,
                                          T_DISCRETIZATION::BASE::View_Convert(solver->U().x) );
    }
    return TV();
}

ENGINE_INTERFACE_LOCAL::TV  ENGINE_INTERFACE_LOCAL::
Displacement_Mesh( int cell, TV weights ) const {
    if( engineInitialized ){
        return engine->Displacement_Mesh( cell, weights,
                                          T_DISCRETIZATION::BASE::View_Convert(solver->U().x),
                                          T_DISCRETIZATION::View_Convert(solver->U().x_mesh) );
    }
    return TV();
}

void ENGINE_INTERFACE_LOCAL::
Displacement( ARRAY<TRIPLE<int, T_INDEX, TV> >& queue, ARRAY<TV>& updates) const {
    updates.Resize( queue.m );

    if( engineInitialized ){
        for( int q =1; q <= queue.m; q++){
            if(queue(q).x == 0)
                updates(q) = engine->Displacement_Grid( queue(q).y, queue(q).z,
                                                        T_DISCRETIZATION::BASE::View_Convert(solver->U().x) );
            else
                updates(q) = engine->Displacement_Mesh( queue(q).x, queue(q).z,
                                                        T_DISCRETIZATION::BASE::View_Convert(solver->U().x),
                                                        T_DISCRETIZATION::View_Convert(solver->U().x_mesh));
        }
    } 
}

void ENGINE_INTERFACE_LOCAL::
Deformation( ARRAY<TRIPLE<int, T_INDEX, TV> >& queue, ARRAY<TV>& updates) const {
    updates.Resize( queue.m );

    if( engineInitialized ){
        for( int q =1; q <= queue.m; q++){
            if(queue(q).x == 0)
                updates(q) = engine->Deformation_Grid( queue(q).y, queue(q).z,
                                                        T_DISCRETIZATION::BASE::View_Convert(solver->U().x) );
            else
                updates(q) = engine->Deformation_Mesh( queue(q).x, queue(q).z,
                                                        T_DISCRETIZATION::BASE::View_Convert(solver->U().x),
                                                        T_DISCRETIZATION::View_Convert(solver->U().x_mesh));
        }
    }
}


void ENGINE_INTERFACE_LOCAL::
Stress( ARRAY<TRIPLE<int, T_INDEX, TV> >& queue, ARRAY<TV>& updates) const {
    updates.Resize( queue.m );

    if( engineInitialized ){
        for( int q =1; q <= queue.m; q++){
            if(queue(q).x == 0){
                updates(q) = engine->Stress_Grid(queue(q).y, queue(q).z);
            }
            else{
                updates(q) = engine->Stress_Mesh(queue(q).x, queue(q).z);
            }
        }
    }
}

void ENGINE_INTERFACE_LOCAL::
Strain( ARRAY<TRIPLE<int, T_INDEX, TV> >& queue, ARRAY<TV>& updates) const {
    updates.Resize( queue.m );

    if( engineInitialized ){
        for( int q =1; q <= queue.m; q++){
            if(queue(q).x == 0){
                updates(q) = engine->Strain_Grid(queue(q).y, queue(q).z);
            }
            else{
                updates(q) = engine->Strain_Mesh(queue(q).x, queue(q).z);
            }
        }
    }
}


ENGINE_INTERFACE_LOCAL::T   ENGINE_INTERFACE_LOCAL::
h() const {
    if( engineCreated )
        return engine->h;
    return 1.0;
}

ENGINE_INTERFACE_LOCAL::T_INDEX  ENGINE_INTERFACE_LOCAL::
Cell(TV& location) const {
    if( engineCreated )
        return engine->grid.Cell( location, 0 );
    return T_INDEX();
}

ENGINE_INTERFACE_LOCAL::TV  ENGINE_INTERFACE_LOCAL::
Node(T_INDEX& location) const {
    if( engineCreated )
        return engine->grid.Node( location );
    return TV();    
}

RANGE<ENGINE_INTERFACE_LOCAL::T_INDEX>  ENGINE_INTERFACE_LOCAL::
Padded_Node_Domain() const {
    if( engineCreated )
        return engine->padded_node_domain;
    return RANGE<T_INDEX>();        
}

RANGE<ENGINE_INTERFACE_LOCAL::T_INDEX>  ENGINE_INTERFACE_LOCAL::
Unpadded_Node_Domain() const {
    if( engineCreated )
        return engine->unpadded_node_domain;
    return RANGE<T_INDEX>();    
}

RANGE<ENGINE_INTERFACE_LOCAL::T_INDEX>  ENGINE_INTERFACE_LOCAL::
Padded_Cell_Domain() const {
    if( engineCreated )
        return engine->padded_cell_domain;
    return RANGE<T_INDEX>();
}

RANGE<ENGINE_INTERFACE_LOCAL::T_INDEX>  ENGINE_INTERFACE_LOCAL::
Unpadded_Cell_Domain() const {
    if( engineCreated )
        return engine->unpadded_cell_domain;
    return RANGE<T_INDEX>();
}

GRID<ENGINE_INTERFACE_LOCAL::TV>  ENGINE_INTERFACE_LOCAL::
Grid() const {
    if( engineCreated )
        return engine->grid;
    return GRID<TV>();    
}

void  ENGINE_INTERFACE_LOCAL::
GetCoarseGrid( VECTOR< ARRAY< T, T_INDEX >, d>& grid_u ) const {
    if( engineInitialized )
        for( int v=1; v<=d; v++)
            grid_u(v) = solver->U().x(v);
}

void  ENGINE_INTERFACE_LOCAL::
GetCellType(ARRAY<CELL_TYPE, T_INDEX>& cell_type) const {
    if( engineCreated )
        cell_type = engine->cell_type;
}

void  ENGINE_INTERFACE_LOCAL::
SetCellType(ARRAY<CELL_TYPE, T_INDEX>& cell_type) {
    if( engineCreated )
        engine->cell_type = cell_type;
}

void  ENGINE_INTERFACE_LOCAL::
GetCellTypeMesh(ARRAY<CELL_TYPE, int>& cell_type_mesh) const {
    if( engineCreated )
        cell_type_mesh = engine->cell_type_mesh;
}

void  ENGINE_INTERFACE_LOCAL::
SetCellTypeMesh(ARRAY<CELL_TYPE, int>& cell_type_mesh) {
    if( engineCreated )
        engine->cell_type_mesh = cell_type_mesh;
}

void  ENGINE_INTERFACE_LOCAL::
GetCellIndicesMesh(ARRAY<T_INDEX,  int>& cell_indices_mesh) const {
    if( engineCreated )
        cell_indices_mesh = engine->cell_indices_mesh;
}

void  ENGINE_INTERFACE_LOCAL::
SetCellIndicesMesh(ARRAY<T_INDEX, int>& cell_indices_mesh) {
    if( engineCreated )
        engine->cell_indices_mesh = cell_indices_mesh;
}

void  ENGINE_INTERFACE_LOCAL::
GetCellsMesh(ARRAY<VECTOR<int, 8>, int>& cells_mesh) const {
    if( engineCreated )
        cells_mesh = engine->cells_mesh;
}

void  ENGINE_INTERFACE_LOCAL::
SetCellsMesh(ARRAY<VECTOR<int, 8>, int>& cells_mesh) {
    if( engineCreated )
        engine->cells_mesh = cells_mesh;
}

int  ENGINE_INTERFACE_LOCAL::
Number_Of_Mesh_Nodes() const {
    if( engineCreated )
        return engine->Number_Of_Mesh_Nodes();
    return 0;
}

int  ENGINE_INTERFACE_LOCAL::
Number_Of_Mesh_Cells() const {
    if( engineCreated )
        return engine->Number_Of_Mesh_Cells();
    return 0;
}


bool ENGINE_INTERFACE_LOCAL::
Node_Is_Dirichlet(const T_INDEX& index) const {
    if( engineCreated )
        return engine->node_is_dirichlet(index);
    return false;
}

bool ENGINE_INTERFACE_LOCAL::
Node_Is_Active(const T_INDEX& index) const {
    if( engineCreated )
        return engine->node_is_active(index);
    return false;
}

bool ENGINE_INTERFACE_LOCAL::
Node_Is_Dirichlet_Mesh(const int& index) const {
    if( engineCreated )
        return engine->node_is_dirichlet_mesh(index);
    return false;
}

bool ENGINE_INTERFACE_LOCAL::
Node_Is_Active_Mesh(const int& index) const {
    if( engineCreated )
        return engine->node_is_active_mesh(index);
    return false;
}

VECTOR< ARRAY_VIEW< ENGINE_INTERFACE_LOCAL::T, ENGINE_INTERFACE_LOCAL::T_INDEX >, ENGINE_INTERFACE_LOCAL::d >& ENGINE_INTERFACE_LOCAL::
U(){
    if( engineCreated )
        return solver->U().x;
}

VECTOR< ARRAY_VIEW< ENGINE_INTERFACE_LOCAL::T, int >, ENGINE_INTERFACE_LOCAL::d >& ENGINE_INTERFACE_LOCAL::
U_Mesh(){
    if( engineCreated )
        return solver->U().x_mesh;
}

void ENGINE_INTERFACE_LOCAL::
Output_Structures(GEOMETRY_PARTICLES<TV>& particles,ARRAY<STRUCTURE<TV>*>& collection_structures) const{
    if( engineCreated ){
#ifdef ENABLE_PHYSBAM_IO
        Create_Output_Data_Lattice(*engine, solver->U(), particles, collection_structures);
#endif
    }
}

int  ENGINE_INTERFACE_LOCAL::
CreateNewConstraint(CONSTRAINT_TYPE ctype) {
    if( engineCreated ){
        int cid;
        switch( ctype ){
        case DYNAMIC:
            cid = engine->dynamic_point_constraints.Append(CONSTRAINT_SEGMENT<T,d>());
            break;
        case STATIC:
            cid = engine->static_point_constraints.Append(CONSTRAINT_SEGMENT<T,d>());
            break;
        case COLLISION:
            cid = engine->collision_constraints.Append(CONSTRAINT_SEGMENT<T,d>());
            break;
        }
        return cid;
    }
    return -1;
}

void  ENGINE_INTERFACE_LOCAL::
GetConstraint(CONSTRAINT_TYPE ctype, int cid, CONSTRAINT_SEGMENT<T,d>& cs) const {
    if(engineCreated){
        switch( ctype ){
        case DYNAMIC:
            cs = engine->dynamic_point_constraints( cid );
            break;
        case  STATIC:
            cs = engine->static_point_constraints( cid );
            break;
        case COLLISION:
            cs = engine->collision_constraints( cid );
            break;
        }
    }
}

void  ENGINE_INTERFACE_LOCAL::
SetConstraint(CONSTRAINT_TYPE ctype, int cid, CONSTRAINT_SEGMENT<T,d>& cs) {
    if(engineCreated){
        switch( ctype ){
        case DYNAMIC:
            engine->dynamic_point_constraints( cid ) = cs;
            break;
        case  STATIC:
            engine->static_point_constraints( cid ) = cs;
            break;
        case COLLISION:
            engine->collision_constraints( cid ) = cs;
            break;
        }
    }
}

void  ENGINE_INTERFACE_LOCAL::
RemoveConstraint(CONSTRAINT_TYPE ctype, int& cid) {
    if( engineCreated ){
        switch( ctype ){
        case DYNAMIC:
            {ARRAY<CONSTRAINT_SEGMENT<T,d> >& point_constraints = engine->dynamic_point_constraints;
            if(cid!=point_constraints.m){
                int other_cid=point_constraints.m;
                PHYSBAM_ASSERT(other_cid!=0);
                exchange(point_constraints(cid),point_constraints(other_cid));
                cid = other_cid;
            }
            // Eliminate last constraint (and hook id)
            point_constraints.Remove_End();            
            }
            break;
        case STATIC:
            {ARRAY<CONSTRAINT_SEGMENT<T,d> >& point_constraints = engine->static_point_constraints;
            if(cid!=point_constraints.m){
                int other_cid=point_constraints.m;
                PHYSBAM_ASSERT(other_cid!=0);
                exchange(point_constraints(cid),point_constraints(other_cid));
                cid = other_cid;                            
            }
            // Eliminate last constraint (and hook id)
            point_constraints.Remove_End();
            }
            break;
        case COLLISION:
            {ARRAY<CONSTRAINT_SEGMENT<T,d> >& point_constraints = engine->collision_constraints;
            if(cid!=point_constraints.m){
                int other_cid=point_constraints.m;
                PHYSBAM_ASSERT(other_cid!=0);
                exchange(point_constraints(cid),point_constraints(other_cid));
                cid = other_cid;
            }
            // Eliminate last constraint (and hook id)
            point_constraints.Remove_End();
            }
            break;
        }
    }
}

int  ENGINE_INTERFACE_LOCAL::
NumberOfConstraints(CONSTRAINT_TYPE ctype) const {
    if( engineCreated ){
        switch( ctype ){
        case DYNAMIC:
            return engine->dynamic_point_constraints.m;
            break;
        case  STATIC:
            return engine->static_point_constraints.m;
            break;
        case COLLISION:
            return engine->collision_constraints.m;
            break;
        } 
    }
    return 0;
}

void  ENGINE_INTERFACE_LOCAL::
GetCollisionConstants(ARRAY<T>& collision_constants) const {
    if( engineCreated )
        collision_constants = engine->collision_spring_constants;
}

void  ENGINE_INTERFACE_LOCAL::
SetCollisionConstants(ARRAY<T>& collision_constants) {
    if( engineCreated )
        engine->collision_spring_constants = collision_constants;
}

void  ENGINE_INTERFACE_LOCAL::
GetCollisionLocations(ARRAY<TV>& collision_locations) const {
    if( engineCreated )
        collision_locations = engine->collision_spring_locations;
}

void  ENGINE_INTERFACE_LOCAL::
SetCollisionLocations(ARRAY<TV>& collision_locations) {
    if( engineCreated )
        engine->collision_spring_locations = collision_locations;
}

void ENGINE_INTERFACE_LOCAL::
GetMuscleData(int& max_muscle,
              ARRAY< ARRAY<int>, T_INDEX >& muscle_ids,
              ARRAY< ARRAY<int>, int >& muscle_ids_mesh, 
              ARRAY< ARRAY<T>, T_INDEX >& muscle_density,
              ARRAY< ARRAY<T>, int >& muscle_density_mesh, 
              ARRAY< ARRAY<TV>, T_INDEX >& muscle_fiber,
              ARRAY< ARRAY<TV>, int >& muscle_fiber_mesh ){

    if( engineCreated ){
        max_muscle = engine->muscle_activations.m;
        muscle_ids = engine->cell_muscles;
        muscle_ids_mesh = engine->cell_muscles_mesh;
        muscle_density = engine->cell_densities;
        muscle_density_mesh = engine->cell_densities_mesh;
        muscle_fiber = engine->cell_fibers;
        muscle_fiber_mesh = engine->cell_fibers_mesh;
    }
}

void ENGINE_INTERFACE_LOCAL::
SetMuscleData(int& max_muscle,
              ARRAY< ARRAY<int>, T_INDEX>& muscle_ids,
              ARRAY< ARRAY<int>, int>& muscle_ids_mesh, 
              ARRAY< ARRAY<T>, T_INDEX>& muscle_density,
              ARRAY< ARRAY<T>, int>& muscle_density_mesh, 
              ARRAY< ARRAY<TV>, T_INDEX>& muscle_fiber,
              ARRAY< ARRAY<TV>, int>& muscle_fiber_mesh)
{
    
    if( engineCreated ){
        engine->muscle_activations.Resize( max_muscle );
        engine->muscle_fiber_max_stresses.Resize( max_muscle );
        engine->cell_muscles          = muscle_ids;
        engine->cell_muscles_mesh     = muscle_ids_mesh;
        engine->cell_densities        = muscle_density;      
        engine->cell_densities_mesh   = muscle_density_mesh;
        engine->cell_fibers           = muscle_fiber;        
        engine->cell_fibers_mesh      = muscle_fiber_mesh;
    }

}

void ENGINE_INTERFACE_LOCAL::
GetMuscleActivations( ARRAY<T, int>& muscle_activations ){
    if(engineCreated)
        muscle_activations = engine->muscle_activations;
}

void ENGINE_INTERFACE_LOCAL::
SetMuscleActivations( ARRAY<T, int>& muscle_activations ){
    if(engineCreated)
        engine->muscle_activations = muscle_activations;
}

void ENGINE_INTERFACE_LOCAL::
GetMuscleFiberMaxStress( ARRAY<T, int>& fiber_max_stress ){
    if(engineCreated)
        fiber_max_stress = engine->muscle_fiber_max_stresses;
}

void ENGINE_INTERFACE_LOCAL::
SetMuscleFiberMaxStress( ARRAY<T, int>& fiber_max_stress ){
    if(engineCreated)
        engine->muscle_fiber_max_stresses = fiber_max_stress;
}

bool ENGINE_INTERFACE_LOCAL::EngineReady() const{
    return engineCreated;   
}
    
bool ENGINE_INTERFACE_LOCAL::EngineInitialized(){
    return engineInitialized;
}

}

#ifndef __ENGINE_INTERFACE_H__
#define __ENGINE_INTERFACE_H__
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Data_Structures/TRIPLE.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/STRUCTURE.h>
#include <PhysBAM_Geometry/Geometry_Particles/GEOMETRY_PARTICLES.h>
#include "CONSTRAINTS.h"
#include "CELL_TYPE.h"

namespace PhysBAM{

    class ENGINE_INTERFACE
    {

    protected:
        typedef float T;
        static const int d = 3;
        
        typedef VECTOR<T,d> TV;
        typedef VECTOR<int,d> T_INDEX;
        typedef bool T_FLAG;

        // Datatypes
    public:
        typedef enum {DYNAMIC=1, STATIC=2, COLLISION=3} CONSTRAINT_TYPE;
        
        // Methods
    public:
        ENGINE_INTERFACE();
        virtual ~ENGINE_INTERFACE();
        
        virtual void CreateEngine(T_INDEX bounds, T dx, TV min_corner, T mu, T lambda, T alpha, T cutoff_value, T stabilization) = 0;
        virtual void Initialize_Muscles() = 0;
        virtual void Initialize_Mesh(int mesh_cell_count, int mesh_node_count) = 0;
        virtual void InitializeEngine() = 0;
        virtual void Exact_Solve( int krylov_iterations, int newton_iterations, T krylov_tolerance, T newton_tolerance, bool no_cut_cells, float& result) = 0;
        virtual void DestroyEngine() = 0;

        virtual TV Displacement_Grid( T_INDEX cell, TV weights ) const = 0;
        virtual TV Displacement_Mesh( int cell, TV weights ) const  = 0;
        virtual void Displacement( ARRAY<TRIPLE<int, T_INDEX, TV> >& queue, ARRAY<TV>& updates) const = 0;
        virtual void Deformation( ARRAY<TRIPLE<int, T_INDEX, TV> >& queue, ARRAY<TV>& updates) const = 0;
        virtual void Stress( ARRAY<TRIPLE<int, T_INDEX, TV> >& queue, ARRAY<TV>& updates) const = 0;
        virtual void Strain( ARRAY<TRIPLE<int, T_INDEX, TV> >& queue, ARRAY<TV>& updates) const = 0;

        //Accessors
    public:
        virtual T h() const  = 0;

        virtual T_INDEX Cell( TV& location) const  = 0;
        virtual TV Node( T_INDEX& location) const  = 0;
    
        virtual RANGE<T_INDEX> Padded_Node_Domain() const = 0;
        virtual RANGE<T_INDEX> Unpadded_Node_Domain() const = 0;
        virtual RANGE<T_INDEX> Padded_Cell_Domain() const = 0;
        virtual RANGE<T_INDEX> Unpadded_Cell_Domain() const = 0;
        virtual GRID<TV> Grid() const  = 0;

        virtual void GetCoarseGrid( VECTOR< ARRAY< T, T_INDEX >, d>& grid_u ) const = 0;
    
        virtual void GetCellType(ARRAY<CELL_TYPE, T_INDEX>& cell_type) const = 0;
        virtual void SetCellType( ARRAY<CELL_TYPE, T_INDEX>& cell_type) = 0;

        virtual void GetCellTypeMesh(ARRAY<CELL_TYPE, int>& cell_type_mesh) const = 0;
        virtual void SetCellTypeMesh( ARRAY<CELL_TYPE, int>& cell_type_mesh) = 0;

        virtual void GetCellIndicesMesh(ARRAY<T_INDEX,  int>& cell_indices_mesh) const = 0;
        virtual void SetCellIndicesMesh( ARRAY<T_INDEX, int>& cell_indices_mesh) = 0;

        virtual void GetCellsMesh(ARRAY<VECTOR<int, 8>, int>& cells_mesh) const = 0;
        virtual void SetCellsMesh( ARRAY<VECTOR<int, 8>, int>& cells_mesh) = 0;

        virtual int Number_Of_Mesh_Nodes() const = 0;
        virtual int Number_Of_Mesh_Cells() const = 0;

        virtual void Output_Structures(GEOMETRY_PARTICLES<TV>& particles,ARRAY<STRUCTURE<TV>*>& collection_structures) const = 0;

        virtual int CreateNewConstraint(CONSTRAINT_TYPE ctype) = 0;
        virtual void GetConstraint(CONSTRAINT_TYPE ctype, int cid, CONSTRAINT_SEGMENT<T,d>& cs) const = 0;
        virtual void SetConstraint(CONSTRAINT_TYPE ctype, int cid,  CONSTRAINT_SEGMENT<T,d>& cs) = 0;
        virtual void RemoveConstraint(CONSTRAINT_TYPE ctype, int cid) = 0;
        virtual int NumberOfConstraints(CONSTRAINT_TYPE ctype) const = 0;

        virtual void GetCollisionConstants(ARRAY<T>& collision_constants) const = 0;
        virtual void SetCollisionConstants( ARRAY<T>& collision_constants) = 0;

        virtual void GetCollisionLocations(ARRAY<TV>& collision_locations) const = 0;
        virtual void SetCollisionLocations( ARRAY<TV>& collision_locations) = 0;

        virtual void GetMuscleData(int& max_muscle,
                                   ARRAY< ARRAY<int>, T_INDEX>& muscle_ids,
                                   ARRAY< ARRAY<int>, int>& muscle_ids_mesh, 
                                   ARRAY< ARRAY<T>, T_INDEX>& muscle_density,
                                   ARRAY< ARRAY<T>, int>& muscle_density_mesh, 
                                   ARRAY< ARRAY<TV>, T_INDEX>& muscle_fiber,
                                   ARRAY< ARRAY<TV>, int>& muscle_fiber_mesh) = 0;

        virtual void SetMuscleData(int& max_muscle,
                                   ARRAY< ARRAY<int>, T_INDEX>& muscle_ids,
                                   ARRAY< ARRAY<int>, int>& muscle_ids_mesh, 
                                   ARRAY< ARRAY<T>, T_INDEX>& muscle_density,
                                   ARRAY< ARRAY<T>, int>& muscle_density_mesh, 
                                   ARRAY< ARRAY<TV>, T_INDEX>& muscle_fiber,
                                   ARRAY< ARRAY<TV>, int>& muscle_fiber_mesh) = 0;

        virtual void GetMuscleActivations( ARRAY<T, int>& muscle_activations ) = 0;
        virtual void SetMuscleActivations( ARRAY<T, int>& muscle_activations ) = 0;

        virtual void GetMuscleFiberMaxStress( ARRAY<T, int>& fiber_max_stress ) = 0;
        virtual void SetMuscleFiberMaxStress( ARRAY<T, int>& fiber_max_stress ) = 0;

        virtual bool EngineReady() const = 0;
        virtual bool EngineInitialized() = 0;

    };

}
#endif

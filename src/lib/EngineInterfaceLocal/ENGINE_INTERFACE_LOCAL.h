#ifndef __ENGINE_INTERFACE_LOCAL_H__
#define __ENGINE_INTERFACE_LOCAL_H__
#include <EngineInterface/ENGINE_INTERFACE.h>
#include <Common_Tools/Math_Tools/RANGE_ITERATOR.h>
#include <Common/STENCIL.h>
#include <Common/STENCIL_ITERATOR.h>

namespace PhysBAM{

    template <class T, int d, bool enable_constraints, bool enable_muscles> class SKINNING_NONLINEAR_ELASTICITY;
    template <class ELASTICITY> class HYBRID_NONLINEAR_ELASTICITY;
    template <class ELASTICITY> struct HYBRID_NONLINEAR_ELASTICITY_STATE;
    template <class ARRAY, int d, int stride> struct HYBRID_BINDER;
    template <class T, class ID, unsigned int align> class ALIGNED_ARRAY;
    template <class BINDER, class STATE, class DISCRETIZATION> class ELASTIC_SOLVER;

    class ENGINE_INTERFACE_LOCAL : public ENGINE_INTERFACE
    {
        using ENGINE_INTERFACE::T;
        using ENGINE_INTERFACE::d;
        using ENGINE_INTERFACE::TV;
        using ENGINE_INTERFACE::T_INDEX;
        using ENGINE_INTERFACE::T_FLAG;
        
        typedef STENCIL<T,d> T_STENCIL;
        typedef STENCIL_ITERATOR<T,d> T_STENCIL_ITERATOR;
        typedef STENCIL_ITERATOR<const T,d> T_CONST_STENCIL_ITERATOR;

        typedef HYBRID_NONLINEAR_ELASTICITY<SKINNING_NONLINEAR_ELASTICITY<T,d,true,true> > T_DISCRETIZATION;
        typedef HYBRID_NONLINEAR_ELASTICITY_STATE<SKINNING_NONLINEAR_ELASTICITY<T,d,true,true> > T_STATE;
        typedef HYBRID_BINDER<ALIGNED_ARRAY<T, int, 64>, d, 16> T_STATE_BINDER;
        
        T_DISCRETIZATION* engine;
        ELASTIC_SOLVER<T_STATE_BINDER, T_STATE, T_DISCRETIZATION>* solver;

        bool engineCreated;
        bool engineInitialized;

        // Datatypes
    public:
        using ENGINE_INTERFACE::CONSTRAINT_TYPE;
        
        // Methods
    public:
        ENGINE_INTERFACE_LOCAL();
        virtual ~ENGINE_INTERFACE_LOCAL();
        
        virtual void CreateEngine(T_INDEX bounds, T dx, TV min_corner, T mu, T lambda, T alpha, T cutoff_value, T stabilization);
        virtual void Initialize_Muscles();
        virtual void Initialize_Mesh(int mesh_cell_count, int mesh_node_count);
        virtual void InitializeEngine();
        virtual void Exact_Solve( int krylov_iterations, int newton_iterations, T krylov_tolerance, T newton_tolerance, bool no_cut_cells, float& result);
        virtual void DestroyEngine();

        virtual TV Displacement_Grid( T_INDEX cell, TV weights ) const;
        virtual TV Displacement_Mesh( int cell, TV weights ) const ;
        virtual void Displacement( ARRAY<TRIPLE<int, T_INDEX, TV> >& queue, ARRAY<TV>& updates) const;
        virtual void Deformation( ARRAY<TRIPLE<int, T_INDEX, TV> >& queue, ARRAY<TV>& updates) const;
        virtual void Stress( ARRAY<TRIPLE<int, T_INDEX, TV> >& queue, ARRAY<TV>& updates) const;
        virtual void Strain( ARRAY<TRIPLE<int, T_INDEX, TV> >& queue, ARRAY<TV>& updates) const;

        //Accessors
    public:

        virtual T h() const ;

        virtual T_INDEX Cell( TV& location) const ;
        virtual TV Node( T_INDEX& location) const ;
    
        virtual RANGE<T_INDEX> Padded_Node_Domain() const;
        virtual RANGE<T_INDEX> Unpadded_Node_Domain() const;
        virtual RANGE<T_INDEX> Padded_Cell_Domain() const;
        virtual RANGE<T_INDEX> Unpadded_Cell_Domain() const;
        virtual GRID<TV> Grid() const ;

        virtual void GetCoarseGrid( VECTOR< ARRAY< T, T_INDEX >, d>& grid_u ) const;
    
        virtual void GetCellType(ARRAY<CELL_TYPE, T_INDEX>& cell_type) const;
        virtual void SetCellType( ARRAY<CELL_TYPE, T_INDEX>& cell_type);

        virtual void GetCellTypeMesh(ARRAY<CELL_TYPE, int>& cell_type_mesh) const;
        virtual void SetCellTypeMesh( ARRAY<CELL_TYPE, int>& cell_type_mesh);

        virtual void GetCellIndicesMesh(ARRAY<T_INDEX,  int>& cell_indices_mesh) const;
        virtual void SetCellIndicesMesh( ARRAY<T_INDEX, int>& cell_indices_mesh);

        virtual void GetCellsMesh(ARRAY<VECTOR<int, 8>, int>& cells_mesh) const;
        virtual void SetCellsMesh( ARRAY<VECTOR<int, 8>, int>& cells_mesh);

        virtual int Number_Of_Mesh_Nodes() const;
        virtual int Number_Of_Mesh_Cells() const;

        bool Node_Is_Dirichlet(const T_INDEX& index) const;
        bool Node_Is_Active(const T_INDEX& index) const;
        bool Node_Is_Dirichlet_Mesh(const int& index) const;
        bool Node_Is_Active_Mesh(const int& index) const;
        VECTOR< ARRAY_VIEW< T, T_INDEX >, d >& U();
        VECTOR< ARRAY_VIEW< T, int >, d >& U_Mesh();
        virtual void Output_Structures(GEOMETRY_PARTICLES<TV>& particles,ARRAY<STRUCTURE<TV>*>& collection_structures) const;

        virtual int CreateNewConstraint(CONSTRAINT_TYPE ctype);
        virtual void GetConstraint(CONSTRAINT_TYPE ctype, int cid, CONSTRAINT_SEGMENT<T,d>& cs) const;
        virtual void SetConstraint(CONSTRAINT_TYPE ctype, int cid,  CONSTRAINT_SEGMENT<T,d>& cs);
        virtual void RemoveConstraint(CONSTRAINT_TYPE ctype, int& cid);
        virtual int NumberOfConstraints(CONSTRAINT_TYPE ctype) const;

        virtual void GetCollisionConstants(ARRAY<T>& collision_constants) const;
        virtual void SetCollisionConstants( ARRAY<T>& collision_constants);

        virtual void GetCollisionLocations(ARRAY<TV>& collision_locations) const;
        virtual void SetCollisionLocations( ARRAY<TV>& collision_locations);

        virtual void GetMuscleData(int& max_muscle,
                                   ARRAY< ARRAY<int>, T_INDEX>& muscle_ids,
                                   ARRAY< ARRAY<int>, int>& muscle_ids_mesh, 
                                   ARRAY< ARRAY<T>, T_INDEX>& muscle_density,
                                   ARRAY< ARRAY<T>, int>& muscle_density_mesh, 
                                   ARRAY< ARRAY<TV>, T_INDEX>& muscle_fiber,
                                   ARRAY< ARRAY<TV>, int>& muscle_fiber_mesh);

        virtual void SetMuscleData(int& max_muscle,
                                   ARRAY< ARRAY<int>, T_INDEX>& muscle_ids,
                                   ARRAY< ARRAY<int>, int>& muscle_ids_mesh, 
                                   ARRAY< ARRAY<T>, T_INDEX>& muscle_density,
                                   ARRAY< ARRAY<T>, int>& muscle_density_mesh, 
                                   ARRAY< ARRAY<TV>, T_INDEX>& muscle_fiber,
                                   ARRAY< ARRAY<TV>, int>& muscle_fiber_mesh);

        virtual void GetMuscleActivations( ARRAY<T, int>& muscle_activations );
        virtual void SetMuscleActivations( ARRAY<T, int>& muscle_activations );

        virtual void GetMuscleFiberMaxStress( ARRAY<T, int>& fiber_max_stress );
        virtual void SetMuscleFiberMaxStress( ARRAY<T, int>& fiber_max_stress );

        virtual bool EngineReady() const;
        virtual bool EngineInitialized();

    };

}
#endif

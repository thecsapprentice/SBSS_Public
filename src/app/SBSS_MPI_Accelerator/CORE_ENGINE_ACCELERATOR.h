#ifndef __CORE_ENGINE_ACCELERATOR_H__
#define __CORE_ENGINE_ACCELERATOR_H__

#include <PhysBAM_Tools/Krylov_Solvers/SYMMQMR.h>
#include <PhysBAM_Tools/Krylov_Solvers/CONJUGATE_GRADIENT.h>
#include <PhysBAM_Tools/Data_Structures/HASHTABLE.h>

#include <Common_Tools/Math_Tools/RANGE_ITERATOR.h>
#include <Common/STENCIL.h>
#include <Common/STENCIL_ITERATOR.h>
#include <Common/ALIGNED_ARRAY.h>

#include <EngineBackend/NONLINEAR_ELASTICITY.h>
#include <EngineBackend/SKINNING_NONLINEAR_ELASTICITY.h>
#include <EngineBackend/HYBRID_NONLINEAR_ELASTICITY.h>

#include <EngineBackend/CG_SYSTEM.h>
#include <EngineBackend/CG_VECTOR.h>

#include <EngineBackend/ELASTIC_SOLVER.h>

#include <EngineInterfaceLocal/ENGINE_INTERFACE_LOCAL.h>

namespace PhysBAM{

class CORE_ENGINE_ACCELERATOR
{
    typedef float T;
    static const int d = 3;

    typedef VECTOR<T,d> TV;
    typedef VECTOR<int,d> T_INDEX;
    typedef bool T_FLAG;

    ENGINE_INTERFACE_LOCAL engineInterface;
    // Methods
 public:
    CORE_ENGINE_ACCELERATOR(int partition_count_input);
    ~CORE_ENGINE_ACCELERATOR();

    void CreateEngine();
    void Initialize_Muscles();
    void Initialize_Mesh();
    void InitializeEngine();
    void Exact_Solve();
    void DestroyEngine();

    void Displacement_Grid() const;
    void Displacement_Mesh() const ;
    void Displacement() const;
    void Deformation() const;
    void Stress() const;
    void Strain() const;

    //Accessors
 public:
    void h() const ;

    void Cell() const ;
    void Node() const ;
    
    void Padded_Node_Domain() const;
    void Unpadded_Node_Domain() const;
    void Padded_Cell_Domain() const;
    void Unpadded_Cell_Domain() const;
    void Grid() const ;

    void GetCoarseGrid() const;

    void GetCellType() const;
    void SetCellType();

    void GetCellTypeMesh() const;
    void SetCellTypeMesh();

    void GetCellIndicesMesh() const;
    void SetCellIndicesMesh();

    void GetCellsMesh() const;
    void SetCellsMesh();

    void Number_Of_Mesh_Nodes() const;
    void Number_Of_Mesh_Cells() const;

    void CreateNewConstraint();
    void GetConstraint() const;
    void SetConstraint();
    void RemoveConstraint();
    void NumberOfConstraints() const;

    void GetCollisionConstants() const;
    void SetCollisionConstants();

    void GetCollisionLocations() const;
    void SetCollisionLocations();

    void GetMuscleData() const;
    void SetMuscleData();
    void GetMuscleActivations() const;
    void SetMuscleActivations();
    void GetMuscleMaxStress() const;
    void SetMuscleMaxStress();

    void EngineReady() const;
    void EngineInitialized();

};

}
#endif

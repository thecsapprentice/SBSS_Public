//#####################################################################
// Copyright 2010, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class HYBRID_NONLINEAR_ELASTICITY
//#####################################################################
#ifndef __HYBRID_NONLINEAR_ELASTICITY__
#define __HYBRID_NONLINEAR_ELASTICITY__

#include "OVERRIDES.h"

#include "SKINNING_NONLINEAR_ELASTICITY.h"
#include <PhysBAM_Tools/Data_Structures/TRIPLE.h>

namespace PhysBAM{


template<class T_ARRAY, int DIM, int stride=16> 
struct HYBRID_BINDER : BINDER<T_ARRAY, DIM, stride>
    {
        static const int STRIDE=stride;
        VECTOR<T_ARRAY,DIM> x_mesh;
        T_ARRAY p_mesh;
    };

template<class T_NONLINEAR_ELASTICITY>
struct HYBRID_NONLINEAR_ELASTICITY_STATE : public CG_POLICY<T_NONLINEAR_ELASTICITY>::STATE
{
    typedef typename CG_POLICY<T_NONLINEAR_ELASTICITY>::STATE T_BASE_STATE;
    typedef typename T_BASE_STATE::SCALAR SCALAR;
    static const int DIM = T_BASE_STATE::DIM;

    typedef ARRAY<SCALAR> T_SCALAR_VARIABLE_MESH_RAW;
    typedef VECTOR<T_SCALAR_VARIABLE_MESH_RAW,DIM> T_VECTOR_VARIABLE_MESH_RAW;

    typedef typename T_BASE_STATE::T_SCALAR_VARIABLE T_SCALAR_VARIABLE_VIEW;
    typedef typename T_BASE_STATE::T_VECTOR_VARIABLE T_VECTOR_VARIABLE_VIEW;
    typedef ARRAY_VIEW<SCALAR> T_SCALAR_VARIABLE_MESH_VIEW;
    typedef VECTOR<T_SCALAR_VARIABLE_MESH_VIEW,DIM> T_VECTOR_VARIABLE_MESH_VIEW;

    using T_BASE_STATE::x;
    using T_BASE_STATE::p;    
    T_VECTOR_VARIABLE_MESH_VIEW x_mesh;
    T_SCALAR_VARIABLE_MESH_VIEW p_mesh;

    HYBRID_NONLINEAR_ELASTICITY_STATE();

    HYBRID_NONLINEAR_ELASTICITY_STATE& operator+=( const HYBRID_NONLINEAR_ELASTICITY_STATE& other ){
        x += other.x;
        p += other.p;
        x_mesh += other.x_mesh;
        p_mesh += other.p_mesh;
        return *this;
    }

    HYBRID_NONLINEAR_ELASTICITY_STATE& operator*=( const HYBRID_NONLINEAR_ELASTICITY_STATE& other ){
        x *= other.x;
        p *= other.p;
        x_mesh *= other.x_mesh;
        p_mesh *= other.p_mesh;
        return *this;
    }

    template<class T_ARRAY, int stride>
    void Initialize(HYBRID_BINDER<T_ARRAY, DIM, stride>& binder, int number_of_mesh_cells, int number_of_mesh_nodes, bool resize_binder=true)
    {
        //LOG::SCOPE scope("HYBRID_NONLINEAR_ELASTICITY_STATE::Initialize()");
        //LOG::cout << "Initializing with data stride: " << stride << std::endl;
        int num_mcells = number_of_mesh_cells + stride*2;
        int num_mnodes = number_of_mesh_nodes + stride*2;

        if(resize_binder){
            for(int i=1;i<=DIM;i++){
                binder.x_mesh(i).Resize( num_mnodes,true,false);
                binder.x_mesh(i).Fill(SCALAR(0));
            }
            //binder.p_mesh.Resize( num_mcells,true,false, SCALAR());
            //binder.p_mesh.Fill(SCALAR(0));
            binder.p_mesh.Resize(0,true,false, SCALAR());
            binder.p_mesh.Fill(SCALAR(0));
        }

        for(int i=1;i<=DIM;i++){
            T_SCALAR_VARIABLE_MESH_VIEW t(number_of_mesh_nodes, binder.x_mesh(i).base_pointer);
            //x(i) = t;
            x_mesh(i).Exchange(t);
        }
        T_SCALAR_VARIABLE_MESH_VIEW t(number_of_mesh_cells, binder.p_mesh.base_pointer);
        //p = t;
        p_mesh.Exchange(t);
    }
};

template<class T_NONLINEAR_ELASTICITY>
class HYBRID_NONLINEAR_ELASTICITY:public T_NONLINEAR_ELASTICITY
{
    // Typedefs and static constants

public:
    typedef T_NONLINEAR_ELASTICITY BASE;
    typedef HYBRID_NONLINEAR_ELASTICITY_STATE<T_NONLINEAR_ELASTICITY> T_STATE;
    typedef typename BASE::SCALAR SCALAR;
    typedef typename BASE::SCALAR T;
    enum {d=BASE::DIM,DIM=BASE::DIM};
protected:
    // From base class : typedefs
    typedef typename BASE::TV TV;
    typedef typename BASE::T_INDEX T_INDEX;
    typedef typename BASE::T_MATERIAL T_MATERIAL;
    typedef typename BASE::T_STENCIL T_STENCIL;

    typedef typename BASE::T_VECTOR_VARIABLE_RAW T_VECTOR_VARIABLE_RAW;
    typedef typename BASE::T_SCALAR_VARIABLE_RAW T_SCALAR_VARIABLE_RAW;
    typedef typename BASE::T_VECTOR_VARIABLE_VIEW T_VECTOR_VARIABLE_VIEW;
    typedef typename BASE::T_SCALAR_VARIABLE_VIEW T_SCALAR_VARIABLE_VIEW;
    typedef typename BASE::T_VECTOR_VARIABLE_VIEW_CONST T_VECTOR_VARIABLE_VIEW_CONST;
    typedef typename BASE::T_SCALAR_VARIABLE_VIEW_CONST T_SCALAR_VARIABLE_VIEW_CONST;

    typedef ARRAY<T,int> T_SCALAR_VARIABLE_MESH_RAW;
    typedef VECTOR<T_SCALAR_VARIABLE_MESH_RAW,d> T_VECTOR_VARIABLE_MESH_RAW;
    typedef ARRAY_VIEW<T,int> T_SCALAR_VARIABLE_MESH_VIEW;
    typedef VECTOR<T_SCALAR_VARIABLE_MESH_VIEW,d> T_VECTOR_VARIABLE_MESH_VIEW;
    typedef ARRAY_VIEW<const T,int> T_SCALAR_VARIABLE_MESH_VIEW_CONST;
    typedef VECTOR<T_SCALAR_VARIABLE_MESH_VIEW_CONST,d> T_VECTOR_VARIABLE_MESH_VIEW_CONST;

    typedef ARRAY<MATRIX<T,d> > T_MATRIX_ARRAY_MESH;
    typedef ARRAY<DIAGONAL_MATRIX<T,d> > T_DIAGONAL_MATRIX_ARRAY_MESH;
    typedef ARRAY<ROTATED_STRESS_DERIVATIVE<T,d> > T_STRESS_DERIVATIVE_ARRAY_MESH;
    typedef STENCIL_ITERATOR<const T,d> T_CONST_STENCIL_ITERATOR; 

    using BASE::vertices_per_cell;
    using BASE::cells_per_block;

    using BASE::supports_constraints;
    using BASE::supports_muscles;

    // Inherited members

public:
    using BASE::h;
    using BASE::constant_partitions;
    using BASE::unpadded_cell_domain;
    using BASE::unpadded_node_domain;
    using BASE::padded_cell_domain;
    using BASE::padded_node_domain;
    using BASE::cell_type;
    using BASE::first_order;
    using BASE::allow_boundary_cells;
public:
    using BASE::G_One;

public:
    using BASE::node_is_active;
    using BASE::node_is_dirichlet;

    ARRAY<int> node_is_active_mesh;
    ARRAY<int> node_is_dirichlet_mesh;

private:
    int number_of_mesh_cells;
    int number_of_mesh_nodes;

public:
    ARRAY<CELL_TYPE> cell_type_mesh;
    ARRAY<T_INDEX> cell_indices_mesh;
    ARRAY<VECTOR<int,vertices_per_cell> > cells_mesh;


 protected:
    T_MATRIX_ARRAY_MESH U_mesh;
    T_MATRIX_ARRAY_MESH V_mesh;
    T_DIAGONAL_MATRIX_ARRAY_MESH Sigma_mesh;
    T_DIAGONAL_MATRIX_ARRAY_MESH Q_hat_mesh;
    T_STRESS_DERIVATIVE_ARRAY_MESH dP_dF_mesh;
 private:
    bool mesh_initialized;


protected:

    typedef T_INDEX HYBRID_BLOCK_BASE;
    typedef ARRAY<VECTOR<int,8> > HYBRID_BLOCK_CELLS;
    typedef ARRAY<VECTOR<int,27> > HYBRID_BLOCK_NODES;
    typedef TRIPLE<HYBRID_BLOCK_BASE, HYBRID_BLOCK_CELLS, HYBRID_BLOCK_NODES > HYBRID_BLOCK_STACK; 
    int number_of_blocks;
    ARRAY<T_INDEX> block_bases;
    HYBRID_BLOCK_CELLS hybrid_block_cells;
    HYBRID_BLOCK_NODES hybrid_block_nodes;
    ARRAY< VECTOR<int, 2>, T_INDEX > hybrid_block_grid_reverse_map;
    ARRAY< VECTOR<int, 2>, int > hybrid_block_mesh_reverse_map;


    // Flattened data structures
    
public:
    using BASE::cell_block_partition_offsets;
    ARRAY<int> cgblock_base_offsets;
    ARRAY<int> cg_offsets;
    ARRAY<int> cg_offsets_node_mesh;
    ARRAY<int> cg_offsets_cell_mesh;

#ifdef USE_SPECIALIZED_KERNELS
    // Hybrid Accleration structures
    ARRAY< PAIR< int, VECTOR<int, 2 > > > reference_constraint_block_ptr;
#ifdef ENABLE_MIC
    int* mic_block_offsets_low;
    int* mic_block_offsets_high;
    int* mic_block_offsets_grid_mask;
    int* mic_block_offsets_mesh_mask;
#endif

#endif

    using BASE::muscle_fiber_max_stresses;
    using BASE::muscle_activations;

    using BASE::cell_muscles;
    ARRAY< ARRAY<int>, int> cell_muscles_mesh;
    using BASE::cell_fibers;
    ARRAY< ARRAY<TV>, int> cell_fibers_mesh;
    using BASE::cell_F_fibers;
    ARRAY< ARRAY<TV>, int> cell_F_fibers_mesh;
    using BASE::cell_densities;
    ARRAY< ARRAY<T>, int> cell_densities_mesh;
    using BASE::cell_c1;
    ARRAY< ARRAY<T>, int> cell_c1_mesh;
    using BASE::cell_c2;
    ARRAY< ARRAY<T>, int> cell_c2_mesh;

#if 0
    // From base class : members

public:
    using BASE::h;
    using BASE::grid;
    using BASE::unpadded_cell_domain;
    using BASE::unpadded_node_domain;
    using BASE::padded_cell_domain;
    using BASE::padded_node_domain;
    using BASE::cell_centered_derivative_operator;
    using BASE::u;

    // New types
protected:
    typedef STENCIL_ITERATOR<const T,d> T_CONST_STENCIL_ITERATOR; 


public:
    // Constraints
    ARRAY<CONSTRAINT_SEGMENT<T,d> > dynamic_point_constraints;
    ARRAY<CONSTRAINT_SEGMENT<T,d> > static_point_constraints;

    // Collisions 
    ARRAY<CONSTRAINT_SEGMENT<T,d> > collision_constraints;
    ARRAY<T> collision_spring_constants;
    ARRAY<TV> collision_spring_locations;
#endif
//#####################################################################
public:
    // Overloaded base methods
    HYBRID_NONLINEAR_ELASTICITY(const T_INDEX& n_input,const T h_input,const TV origin);

    static T_VECTOR_VARIABLE_MESH_VIEW_CONST View_Convert(const T_VECTOR_VARIABLE_MESH_VIEW& in);

    void Initialize_Domain_Structures();    
    void Initialize_Mesh(const int number_of_mesh_cells_input,const int number_of_mesh_nodes_input);
    void Initialize_Muscles();
    void Initialize_Blocks(const int number_of_partitions);
    void Initialize_Blocks_Elasticity(const int number_of_partitions);
    void Initialize_Blocks_Muscles(const int number_of_partitions);
    void Initialize_Blocks_Constraints(const int number_of_partitions);
    bool Cell_Has_Constraint(const T_INDEX& grid_cell, int mesh_cell );

    template<class T_ARRAY, int stride> 
    void Initialize_State(HYBRID_BINDER<T_ARRAY, d, stride>& binder, T_STATE& state, bool resize_binder=true) const
        {
            PHYSBAM_ASSERT(mesh_initialized==true);
            BASE::Initialize_State(binder, state, resize_binder);
            state.Initialize(binder, number_of_mesh_cells, number_of_mesh_nodes, resize_binder);
        }

    void CompactData_Specialized(T_VECTOR_VARIABLE_VIEW_CONST x,T_SCALAR_VARIABLE_VIEW_CONST p,T_VECTOR_VARIABLE_MESH_VIEW_CONST x_mesh,T_SCALAR_VARIABLE_MESH_VIEW_CONST p_mesh) const;
    void UnCompactData_Specialized(T_VECTOR_VARIABLE_VIEW x,T_SCALAR_VARIABLE_VIEW p, T_VECTOR_VARIABLE_MESH_VIEW x_mesh,T_SCALAR_VARIABLE_MESH_VIEW p_mesh) const;
    void CompactData_Specialized(T_VECTOR_VARIABLE_VIEW_CONST x,T_VECTOR_VARIABLE_MESH_VIEW_CONST x_mesh) const;
    void UnCompactData_Specialized(T_VECTOR_VARIABLE_VIEW x, T_VECTOR_VARIABLE_MESH_VIEW x_mesh) const;


    static bool Supports_Constraints() {return supports_constraints;}
    static bool Supports_Muscles() {return supports_muscles;}

    void Update_Position_Based_State(const T_STATE& state_u, T_STATE& state_d);
    void Update_Position_Based_State_Muscles_And_Constraints(const T_STATE& state_u);
    void Update_Position_Based_State_Elasticity_Mesh(T_VECTOR_VARIABLE_VIEW_CONST u,T_SCALAR_VARIABLE_VIEW_CONST p,T_VECTOR_VARIABLE_VIEW diag, T_VECTOR_VARIABLE_MESH_VIEW_CONST u_mesh,T_SCALAR_VARIABLE_MESH_VIEW_CONST p_mesh, T_VECTOR_VARIABLE_MESH_VIEW diag_mesh);
    void Update_Position_Based_State_Constraints_Mesh(T_VECTOR_VARIABLE_VIEW_CONST u,T_SCALAR_VARIABLE_VIEW_CONST p,T_VECTOR_VARIABLE_MESH_VIEW_CONST u_mesh,T_SCALAR_VARIABLE_MESH_VIEW_CONST p_mesh);
    void Update_Position_Based_State_Muscles_Mesh(T_VECTOR_VARIABLE_VIEW_CONST u,T_SCALAR_VARIABLE_VIEW_CONST p,T_VECTOR_VARIABLE_MESH_VIEW_CONST u_mesh,T_SCALAR_VARIABLE_MESH_VIEW_CONST p_mesh);
    void Update_Position_Based_State_Specialized(T_VECTOR_VARIABLE_VIEW_CONST u,T_SCALAR_VARIABLE_VIEW_CONST p,T_VECTOR_VARIABLE_VIEW diag,T_VECTOR_VARIABLE_MESH_VIEW_CONST u_mesh,T_SCALAR_VARIABLE_MESH_VIEW_CONST p_mesh,T_VECTOR_VARIABLE_MESH_VIEW diag_mesh);
    //void Update_Collision_Constraints(T_VECTOR_VARIABLE_VIEW_CONST u, T_VECTOR_VARIABLE_MESH_VIEW_CONST u_mesh);

    void Add_Force(const T_STATE& state_u, T_STATE& state_f);
    void Add_Force_Muscles_And_Constraints(const T_STATE& state_u, T_STATE& state_f) const;
    void Add_Force_First_Order_Elasticity_Mesh(T_VECTOR_VARIABLE_VIEW_CONST u,T_SCALAR_VARIABLE_VIEW_CONST p,T_VECTOR_VARIABLE_MESH_VIEW_CONST u_mesh,T_SCALAR_VARIABLE_MESH_VIEW_CONST p_mesh,
                                               T_VECTOR_VARIABLE_VIEW f,T_SCALAR_VARIABLE_VIEW q,T_VECTOR_VARIABLE_MESH_VIEW f_mesh,T_SCALAR_VARIABLE_MESH_VIEW q_mesh) const;
    void Add_Force_Constraints_Mesh(const ARRAY<CONSTRAINT_SEGMENT<T,d> >& constraints,
                                    T_VECTOR_VARIABLE_VIEW_CONST u,T_SCALAR_VARIABLE_VIEW_CONST p,T_VECTOR_VARIABLE_MESH_VIEW_CONST u_mesh,T_SCALAR_VARIABLE_MESH_VIEW_CONST p_mesh,
                                    T_VECTOR_VARIABLE_VIEW f,T_SCALAR_VARIABLE_VIEW q,T_VECTOR_VARIABLE_MESH_VIEW f_mesh,T_SCALAR_VARIABLE_MESH_VIEW q_mesh) const;
    void Add_Force_Muscles_Mesh(T_VECTOR_VARIABLE_VIEW_CONST u,T_SCALAR_VARIABLE_VIEW_CONST p,T_VECTOR_VARIABLE_MESH_VIEW_CONST u_mesh,T_SCALAR_VARIABLE_MESH_VIEW_CONST p_mesh,
                                T_VECTOR_VARIABLE_VIEW f,T_SCALAR_VARIABLE_VIEW q,T_VECTOR_VARIABLE_MESH_VIEW f_mesh,T_SCALAR_VARIABLE_MESH_VIEW q_mesh) const;
    void Add_Force_First_Order_Elasticity_Specialized(T_VECTOR_VARIABLE_VIEW_CONST u,T_SCALAR_VARIABLE_VIEW_CONST p,T_VECTOR_VARIABLE_MESH_VIEW_CONST u_mesh,T_SCALAR_VARIABLE_MESH_VIEW_CONST p_mesh,
                                                      T_VECTOR_VARIABLE_VIEW f,T_SCALAR_VARIABLE_VIEW q,T_VECTOR_VARIABLE_MESH_VIEW f_mesh,T_SCALAR_VARIABLE_MESH_VIEW q_mesh) ;
    void Add_Force_Stab_Mesh(T_VECTOR_VARIABLE_VIEW_CONST dx,T_VECTOR_VARIABLE_VIEW df,T_VECTOR_VARIABLE_MESH_VIEW_CONST dx_mesh,T_VECTOR_VARIABLE_MESH_VIEW df_mesh) const;

    void Add_Force_Differential(const T_STATE& state_du, T_STATE& state_df) const;
    void Add_Force_Differential(const T_STATE& state_du, T_STATE& state_df, const int subdomain) const;
    void Add_Force_Differential_Muscles_And_Constraints(const T_STATE& state_du, T_STATE& state_df) const;
    void Add_Force_Differential_Elasticity_Mesh(T_VECTOR_VARIABLE_VIEW_CONST du, T_SCALAR_VARIABLE_VIEW_CONST dp,T_VECTOR_VARIABLE_MESH_VIEW_CONST du_mesh, T_SCALAR_VARIABLE_MESH_VIEW_CONST dp_mesh,
                                                T_VECTOR_VARIABLE_VIEW df,T_SCALAR_VARIABLE_VIEW dq,T_VECTOR_VARIABLE_MESH_VIEW df_mesh,T_SCALAR_VARIABLE_MESH_VIEW dq_mesh) const;
    void Add_Force_Differential_Muscles_Mesh(T_VECTOR_VARIABLE_VIEW_CONST du,T_SCALAR_VARIABLE_VIEW_CONST dp,T_VECTOR_VARIABLE_MESH_VIEW_CONST du_mesh, T_SCALAR_VARIABLE_MESH_VIEW_CONST dp_mesh,T_VECTOR_VARIABLE_VIEW df,T_SCALAR_VARIABLE_VIEW dq,T_VECTOR_VARIABLE_MESH_VIEW df_mesh,T_SCALAR_VARIABLE_MESH_VIEW dq_mesh) const;
    void Add_Force_Differential_Constraints_Mesh(const ARRAY<CONSTRAINT_SEGMENT<T,d> >& constraints,T_VECTOR_VARIABLE_VIEW_CONST du,T_SCALAR_VARIABLE_VIEW_CONST dp,T_VECTOR_VARIABLE_MESH_VIEW_CONST du_mesh, T_SCALAR_VARIABLE_MESH_VIEW_CONST dp_mesh,T_VECTOR_VARIABLE_VIEW df,T_SCALAR_VARIABLE_VIEW dq,T_VECTOR_VARIABLE_MESH_VIEW df_mesh,T_SCALAR_VARIABLE_MESH_VIEW dq_mesh) const;
    void Add_Force_Differential_Elasticity_Specialized(T_VECTOR_VARIABLE_VIEW_CONST du,T_SCALAR_VARIABLE_VIEW_CONST dp,T_VECTOR_VARIABLE_MESH_VIEW_CONST du_mesh,T_SCALAR_VARIABLE_MESH_VIEW_CONST dp_mesh,T_VECTOR_VARIABLE_VIEW df,T_SCALAR_VARIABLE_VIEW dq, T_VECTOR_VARIABLE_MESH_VIEW df_mesh,T_SCALAR_VARIABLE_MESH_VIEW dq_mesh) const;
    void Add_Force_Differential_Elasticity_Specialized(T_VECTOR_VARIABLE_VIEW_CONST du,T_SCALAR_VARIABLE_VIEW_CONST dp,T_VECTOR_VARIABLE_MESH_VIEW_CONST du_mesh,T_SCALAR_VARIABLE_MESH_VIEW_CONST dp_mesh,T_VECTOR_VARIABLE_VIEW df,T_SCALAR_VARIABLE_VIEW dq, T_VECTOR_VARIABLE_MESH_VIEW df_mesh,T_SCALAR_VARIABLE_MESH_VIEW dq_mesh, const int subdomain) const;



    int Number_Of_Mesh_Cells() const { return number_of_mesh_cells; };
    int Number_Of_Mesh_Nodes() const { return number_of_mesh_nodes; };

    TV Deformation(const T_INDEX& cell_index,const TV& multilinear_coordinates, const T_STATE& state) const;

    static T_STENCIL Multilinear_Interpolation_Stencil_Grid(const T_INDEX& cell_index,const TV& multilinear_coordinates);
    TV Deformation_Grid(const T_INDEX& cell_index,const TV& multilinear_coordinates,T_VECTOR_VARIABLE_VIEW_CONST u) const;
    static TV Displacement_Grid(const T_INDEX& cell_index,const TV& multilinear_coordinates,T_VECTOR_VARIABLE_VIEW_CONST du);

    static T_STENCIL Multilinear_Interpolation_Stencil_Mesh(const TV& multilinear_coordinates);
    TV Deformation_Mesh(const int& mesh_element,const TV& multilinear_coordinates,T_VECTOR_VARIABLE_VIEW_CONST u,T_VECTOR_VARIABLE_MESH_VIEW_CONST u_mesh) const;
    TV Displacement_Mesh(const int& mesh_element,const TV& multilinear_coordinates,T_VECTOR_VARIABLE_VIEW_CONST du,T_VECTOR_VARIABLE_MESH_VIEW_CONST du_mesh) const;
    void Initialize_Undeformed_Configuration(T_STATE& state);

    TV Stress_Mesh( const int& mesh_element, const TV& multilinear_coordinates );
    TV Stress_Grid( const T_INDEX& cell_index, const TV& multilinear_coordinates );

    TV Strain_Mesh( const int& mesh_element, const TV& multilinear_coordinates );
    TV Strain_Grid( const T_INDEX& cell_index, const TV& multilinear_coordinates );


#if 0

    void Initialize_Blocks(const int number_of_partitions);
    void Initialize_Blocks_Muscles(const int number_of_partitions);
    void Initialize_Blocks_Constraints(const int number_of_partitions);

    virtual void Update_Position_Based_State(const bool allow_boundary_cells);
    void Update_Position_Based_State_Muscles(const bool allow_boundary_cells);
    void Update_Position_Based_State_Constraints(const bool allow_boundary_cells);
    void Update_Collision_Constraints();

    virtual void Add_Force(T_VECTOR_VARIABLE_VIEW f,T_SCALAR_VARIABLE_VIEW q,const bool first_order) const;
    void Add_Force_Second_Order_Muscles(T_VECTOR_VARIABLE_VIEW f,T_SCALAR_VARIABLE_VIEW q) const;
    void Add_Force_First_Order_Muscles(T_VECTOR_VARIABLE_VIEW f,T_SCALAR_VARIABLE_VIEW q) const;
    void Add_Force_Constraints(const ARRAY<CONSTRAINT_SEGMENT<T,d> >& constraints, T_VECTOR_VARIABLE_VIEW f,T_SCALAR_VARIABLE_VIEW q) const;

    virtual void Add_Force_Differential(T_VECTOR_VARIABLE_VIEW_CONST du,T_SCALAR_VARIABLE_VIEW_CONST dp,T_VECTOR_VARIABLE_VIEW df,T_SCALAR_VARIABLE_VIEW dq) const;
    void Add_Force_Differential_Muscles(T_VECTOR_VARIABLE_VIEW_CONST du,T_SCALAR_VARIABLE_VIEW_CONST dp,T_VECTOR_VARIABLE_VIEW df,T_SCALAR_VARIABLE_VIEW dq) const;
    void Add_Force_Differential_Constraints(const ARRAY<CONSTRAINT_SEGMENT<T,d> >& constraints,T_VECTOR_VARIABLE_VIEW_CONST du,T_SCALAR_VARIABLE_VIEW_CONST dp,T_VECTOR_VARIABLE_VIEW df,T_SCALAR_VARIABLE_VIEW dq) const;

    // New Methods
    void Initialize_Undeformed_Configuration();
    void Initialize_Muscles();
    static T Tension(const T stretch, const T activation, const T density, const T fiber_max_stress);
    static T Tension_Derivative(const T stretch, const T activation, const T density, const T fiber_max_stress);
    //static T Tension(const T_INDEX& cell_index,const int cell_muscle_index,const T stretch);
    //static T Tension_Derivative(const T_INDEX& cell_index,const int cell_muscle_index,const T stretch);
    //void Build_Constraint_Matrix();
    int Add_Embedded_Point_To_Fixed_Point_Spring_Constraint(const T spring_coefficient,const TV& embedded_point_material_space_location,const TV& fixed_point_world_space_location, bool is_static);
    int Add_Two_Embedded_Point_Spring_Constraint(const T spring_coefficient,const TV& embedded_point_material_space_location1,const TV& embedded_point_material_space_location2, bool is_static);
    int Constraint_Count() {return static_point_constraints.m + dynamic_point_constraints.m; };

    static T_STENCIL Multilinear_Interpolation_Stencil(const T_INDEX& cell_index,const TV& multilinear_coordinates);
    TV Deformation(const T_INDEX& cell_index,const TV& multilinear_coordinates) const;
    static TV Displacement(const T_INDEX& cell_index,const TV& multilinear_coordinates,T_VECTOR_VARIABLE_VIEW_CONST du);

#endif
//#####################################################################
};
}
#endif

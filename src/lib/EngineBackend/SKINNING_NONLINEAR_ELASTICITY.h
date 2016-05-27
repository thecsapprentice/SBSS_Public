//#####################################################################
// Copyright 2010, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SKINNING_NONLINEAR_ELASTICITY
//#####################################################################
#ifndef __SKINNING_NONLINEAR_ELASTICITY__
#define __SKINNING_NONLINEAR_ELASTICITY__

#include "OVERRIDES.h"

#include "NONLINEAR_ELASTICITY.h"
#include <EngineInterface/CONSTRAINTS.h>
//#include "COLLISION_INTERFACE.h"

namespace PhysBAM{

template<class T,int d, bool enable_constraints, bool enable_muscles>
class SKINNING_NONLINEAR_ELASTICITY:public NONLINEAR_ELASTICITY<T,d>
{
    typedef NONLINEAR_ELASTICITY<T,d> BASE;

public:
    typedef typename NONLINEAR_ELASTICITY<T,d>::T_STATE T_STATE;
    enum {DIM=d};
    typedef T SCALAR;
    // From base class : typedefs

    static const bool supports_constraints = enable_constraints;
    static const bool supports_muscles = ((d==3) && enable_muscles);

protected:
    typedef typename BASE::TV TV;
    typedef typename BASE::T_INDEX T_INDEX;
    typedef typename BASE::T_STENCIL T_STENCIL;
    typedef typename BASE::T_VECTOR_VARIABLE_VIEW T_VECTOR_VARIABLE_VIEW;
    typedef typename BASE::T_SCALAR_VARIABLE_VIEW T_SCALAR_VARIABLE_VIEW;
    typedef typename BASE::T_VECTOR_VARIABLE_VIEW_CONST T_VECTOR_VARIABLE_VIEW_CONST;
    typedef typename BASE::T_SCALAR_VARIABLE_VIEW_CONST T_SCALAR_VARIABLE_VIEW_CONST;
    // From base class : members

public:
    using BASE::h;
    using BASE::grid;
    using BASE::unpadded_cell_domain;
    using BASE::unpadded_node_domain;
    using BASE::padded_cell_domain;
    using BASE::padded_node_domain;
    using BASE::cell_centered_derivative_operator;
    using BASE::first_order;
    using BASE::allow_boundary_cells;

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

 protected:
    //COLLISION_INTERFACE<T,d>* collision_shape;

//#####################################################################
public:
    // Overloaded base methods
    SKINNING_NONLINEAR_ELASTICITY(const T_INDEX& n_input,const T h_input,const TV origin);
    ~SKINNING_NONLINEAR_ELASTICITY();

    void Initialize_Blocks(const int number_of_partitions);
    void Initialize_Blocks_Muscles(const int number_of_partitions);
    void Initialize_Blocks_Constraints(const int number_of_partitions);
    template<class T_ARRAY, int stride>
        void Initialize_State(BINDER<T_ARRAY, d, stride>& binder, T_STATE& state, bool resize_binder=true) const
        {BASE::Initialize_State(binder, state, resize_binder);}

    static bool Supports_Constraints() {return supports_constraints;}
    static bool Supports_Muscles() {return supports_muscles;}

    void Update_Position_Based_State(const T_STATE& state_u, T_STATE& state_d);
    void Update_Position_Based_State_Muscles(T_VECTOR_VARIABLE_VIEW_CONST u,T_SCALAR_VARIABLE_VIEW_CONST p);
    void Update_Position_Based_State_Constraints(T_VECTOR_VARIABLE_VIEW_CONST u,T_SCALAR_VARIABLE_VIEW_CONST p);
    //void Update_Collision_Constraints(T_VECTOR_VARIABLE_VIEW_CONST u);

    void Add_Force(const T_STATE& state_u, T_STATE& state_f) ;
    void Add_Force_Second_Order_Muscles(T_VECTOR_VARIABLE_VIEW_CONST u,T_SCALAR_VARIABLE_VIEW_CONST p, T_VECTOR_VARIABLE_VIEW f,T_SCALAR_VARIABLE_VIEW q) const;
    void Add_Force_First_Order_Muscles(T_VECTOR_VARIABLE_VIEW_CONST u,T_SCALAR_VARIABLE_VIEW_CONST p, T_VECTOR_VARIABLE_VIEW f,T_SCALAR_VARIABLE_VIEW q) const;
    void Add_Force_Constraints(const ARRAY<CONSTRAINT_SEGMENT<T,d> >& constraints, T_VECTOR_VARIABLE_VIEW_CONST u,T_SCALAR_VARIABLE_VIEW_CONST p, T_VECTOR_VARIABLE_VIEW f,T_SCALAR_VARIABLE_VIEW q) const;

    void Add_Force_Differential(const T_STATE& state_u, T_STATE& state_f) const;
    void Add_Force_Differential_Muscles(T_VECTOR_VARIABLE_VIEW_CONST du,T_SCALAR_VARIABLE_VIEW_CONST dp,T_VECTOR_VARIABLE_VIEW df,T_SCALAR_VARIABLE_VIEW dq) const;
    void Add_Force_Differential_Constraints(const ARRAY<CONSTRAINT_SEGMENT<T,d> >& constraints,T_VECTOR_VARIABLE_VIEW_CONST du,T_SCALAR_VARIABLE_VIEW_CONST dp,T_VECTOR_VARIABLE_VIEW df,T_SCALAR_VARIABLE_VIEW dq) const;

    // New Methods
    void Initialize_Undeformed_Configuration(T_STATE& state);
    void Initialize_Muscles();
    //void Set_Collision_Object(COLLISION_INTERFACE<T,d>* new_collision_object);
    static T Tension(const T stretch, const T activation, const T density, const T fiber_max_stress);
    static T Tension_Derivative(const T stretch, const T activation, const T density, const T fiber_max_stress);
    //static T Tension(const T_INDEX& cell_index,const int cell_muscle_index,const T stretch);
    //static T Tension_Derivative(const T_INDEX& cell_index,const int cell_muscle_index,const T stretch);
    //void Build_Constraint_Matrix();
    int Add_Embedded_Point_To_Fixed_Point_Spring_Constraint(const T spring_coefficient,const TV& embedded_point_material_space_location,const TV& fixed_point_world_space_location, bool is_static);
    int Add_Two_Embedded_Point_Spring_Constraint(const T spring_coefficient,const TV& embedded_point_material_space_location1,const TV& embedded_point_material_space_location2, bool is_static);
    int Constraint_Count() {return static_point_constraints.m + dynamic_point_constraints.m; };
    static T_STENCIL Multilinear_Interpolation_Stencil(const T_INDEX& cell_index,const TV& multilinear_coordinates);

    TV Deformation(const T_INDEX& cell_index, const TV& multilinear_coordinates, const T_STATE& state) const;
    TV Deformation_Grid(const T_INDEX& cell_index,const TV& multilinear_coordinates,T_VECTOR_VARIABLE_VIEW_CONST u) const;
    static TV Displacement_Grid(const T_INDEX& cell_index,const TV& multilinear_coordinates,T_VECTOR_VARIABLE_VIEW_CONST du);

//#####################################################################
};
}
#endif

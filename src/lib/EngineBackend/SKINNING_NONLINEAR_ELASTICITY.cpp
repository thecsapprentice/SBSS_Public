//#####################################################################
// Copyright 2010-2012, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SKINNING_NONLINEAR_ELASTICITY
//#####################################################################
#include "SKINNING_NONLINEAR_ELASTICITY.h"
#include <Common/RANGE_ITERATOR.h>

using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T,int d,bool enable_constraints,bool enable_muscles> SKINNING_NONLINEAR_ELASTICITY<T,d,enable_constraints,enable_muscles>::
SKINNING_NONLINEAR_ELASTICITY(const T_INDEX& n_input,const T h_input,const TV origin)
    :BASE(n_input,h_input,origin)//, collision_shape(NULL)
{}
//#####################################################################
// Destructor
//#####################################################################
template<class T,int d,bool enable_constraints,bool enable_muscles> SKINNING_NONLINEAR_ELASTICITY<T,d,enable_constraints,enable_muscles>::
~SKINNING_NONLINEAR_ELASTICITY(){
//    if( collision_shape != NULL )
//        delete collision_shape; // We can do this as Set_Collision_Shape is what has access to this.
}
//#####################################################################
// Function Initialize_Blocks
//#####################################################################
template<class T,int d,bool enable_constraints,bool enable_muscles> void SKINNING_NONLINEAR_ELASTICITY<T,d,enable_constraints,enable_muscles>::
Initialize_Blocks(const int number_of_partitions)
{
    LOG::SCOPE scope("SKINNING_NONLINEAR_ELASTICITY::Initialize_Blocks()");
    BASE::Initialize_Blocks(number_of_partitions);
    Initialize_Blocks_Constraints(number_of_partitions);
    Initialize_Blocks_Muscles(number_of_partitions);
}

//#####################################################################
// Function Update_Position_Based_State
//#####################################################################
template<class T,int d,bool enable_constraints,bool enable_muscles> void SKINNING_NONLINEAR_ELASTICITY<T,d,enable_constraints,enable_muscles>::
Update_Position_Based_State(const T_STATE& state_u, T_STATE& state_d)
{
    LOG::SCOPE scope("SKINNING_NONLINEAR_ELASTICITY::Update_Position_Based_State()");
#ifdef USE_SPECIALIZED_KERNELS
    NONLINEAR_ELASTICITY<T,d>::template Update_Position_Based_State_Specialized<false,enable_muscles>(BASE::View_Convert(state_u.x), state_u.p, state_d.x);
#else
    NONLINEAR_ELASTICITY<T,d>::Update_Position_Based_State(state_u, state_d);
    Update_Position_Based_State_Muscles(BASE::View_Convert(state_u.x), state_u.p);
#endif
    Update_Position_Based_State_Constraints(BASE::View_Convert(state_u.x), state_u.p);
}
/*
//#####################################################################
// Function Add_Force_First_Order_Nostab
//#####################################################################
template<class T,int d,bool enable_constraints,bool enable_muscles> void SKINNING_NONLINEAR_ELASTICITY<T,d,enable_constraints,enable_muscles>::
Update_Collision_Constraints(T_VECTOR_VARIABLE_VIEW_CONST u)
{
    if( !collision_shape )
        return;

    PHYSBAM_ASSERT( collision_spring_locations.m == collision_constraints.m );
    PHYSBAM_ASSERT( collision_spring_constants.m == collision_constraints.m );

    TV collision_center;
    collision_center(2) = -0.055f;
    //T collision_radius = 0.05f;
    TV collision_normal;
    collision_normal(2) = 1.0f; 

    for( int m=1; m<=collision_constraints.m; m++){
        CONSTRAINT_SEGMENT<T,d>& ic = collision_constraints(m);
        TV p;

        if(ic.endpoints[1].type == CONSTRAINT_NODE<T,d>::GRID_FIXED)
            p = Deformation_Grid(ic.endpoints[1].grid_index(),ic.endpoints[1].multilinear_coordinates(), u);
        //else if(ic.endpoints[1].type == CONSTRAINT_NODE<T,d>::MESH_FIXED)
        //    p = Deformation(ic.endpoints[1].mesh_index(),ic.endpoints[1].multilinear_coordinates);
        else
            PHYSBAM_FATAL_ERROR();

        T depth = collision_shape->Phi( p );

        if( depth < 0 ) { // We are inside the collision shape!
            //LOG::cout << "COLLISION DETECTED!!!" << std::endl;
            //LOG::cout << "Point of Collision: " << p << std::endl;
            //LOG::cout << "New Constraint Point: " << collision_shape->ClosestPoint(p) << std::endl;
#if 1
            collision_spring_locations(m) = collision_shape->ClosestPoint(p);
            collision_spring_constants(m) = 1e6;
#else
            collision_spring_locations(m) = TV();
            collision_spring_constants(m) = 0.0;
#endif
        }
        else{
            collision_spring_locations(m) = TV();
            collision_spring_constants(m) = 0.0;
        }
    }
}

//#####################################################################
// Function Set_Collision_Object
//#####################################################################
template<class T,int d,bool enable_constraints,bool enable_muscles> void SKINNING_NONLINEAR_ELASTICITY<T,d,enable_constraints,enable_muscles>::
Set_Collision_Object(COLLISION_INTERFACE<T,d>* new_collision_object)
{
    if( collision_shape )
        delete collision_shape;
    
    collision_shape = new_collision_object;
}
*/

//#####################################################################
// Function Add_Force
//#####################################################################
template<class T,int d,bool enable_constraints,bool enable_muscles> void SKINNING_NONLINEAR_ELASTICITY<T,d,enable_constraints,enable_muscles>::
Add_Force(const T_STATE& state_u, T_STATE& state_f) 
{
    LOG::SCOPE scope("SKINNING_NONLINEAR_ELASTICITY::Add_Force");

    if(first_order){
#ifdef USE_SPECIALIZED_KERNELS
        NONLINEAR_ELASTICITY<T,d>::template Add_Force_First_Order_Elasticity_Specialized<false,enable_muscles>(BASE::View_Convert(state_u.x), state_u.p, state_f.x,state_f.p);
#else
        NONLINEAR_ELASTICITY<T,d>::Add_Force(state_u, state_f);
        Add_Force_First_Order_Muscles(BASE::View_Convert(state_u.x), state_u.p, state_f.x,state_f.p);
#endif
        Add_Force_Constraints(dynamic_point_constraints, BASE::View_Convert(state_u.x), state_u.p, state_f.x,state_f.p);
        Add_Force_Constraints(static_point_constraints,  BASE::View_Convert(state_u.x), state_u.p, state_f.x,state_f.p);
        Add_Force_Constraints(collision_constraints,     BASE::View_Convert(state_u.x), state_u.p, state_f.x,state_f.p);

    }
    else{
        NONLINEAR_ELASTICITY<T,d>::Add_Force(state_u, state_f);
        Add_Force_Constraints(dynamic_point_constraints, BASE::View_Convert(state_u.x), state_u.p, state_f.x,state_f.p);
        Add_Force_Constraints(static_point_constraints,  BASE::View_Convert(state_u.x), state_u.p, state_f.x,state_f.p);
        Add_Force_Constraints(collision_constraints,     BASE::View_Convert(state_u.x), state_u.p, state_f.x,state_f.p);
        Add_Force_Second_Order_Muscles(BASE::View_Convert(state_u.x), state_u.p, state_f.x,state_f.p);
    }
}
//#####################################################################
// Function Add_Force_Differential
//#####################################################################
template<class T,int d,bool enable_constraints,bool enable_muscles> void SKINNING_NONLINEAR_ELASTICITY<T,d,enable_constraints,enable_muscles>::
Add_Force_Differential(const T_STATE& state_u,T_STATE& state_f) const
{
#ifdef LOG_DETAILED_PERFORMANCE_NE
    LOG::SCOPE scope("SKINNING_NONLINEAR_ELASTICITY::Add_Force_Differential");
#endif
    
#ifdef USE_SPECIALIZED_KERNELS
    NONLINEAR_ELASTICITY<T,d>::template Add_Force_Differential_Elasticity_Specialized<false,enable_muscles>(BASE::View_Convert(state_u.x),state_u.p,state_f.x,state_f.p);
#else
    NONLINEAR_ELASTICITY<T,d>::Add_Force_Differential(state_u,state_f);
    Add_Force_Differential_Muscles(BASE::View_Convert(state_u.x),state_u.p,state_f.x,state_f.p);
#endif
    Add_Force_Differential_Constraints(dynamic_point_constraints,BASE::View_Convert(state_u.x),state_u.p,state_f.x,state_f.p);
    Add_Force_Differential_Constraints(static_point_constraints,BASE::View_Convert(state_u.x),state_u.p,state_f.x,state_f.p);
    Add_Force_Differential_Constraints(collision_constraints,BASE::View_Convert(state_u.x),state_u.p,state_f.x,state_f.p);

}
//#####################################################################
// Function Initialize_Undeformed_Configuration
//#####################################################################
template<class T,int d,bool enable_constraints,bool enable_muscles> void SKINNING_NONLINEAR_ELASTICITY<T,d,enable_constraints,enable_muscles>::
Initialize_Undeformed_Configuration(T_STATE& state)
{
    for(RANGE_ITERATOR<d> iterator(unpadded_node_domain);iterator.Valid();iterator.Next()){
        const T_INDEX& index=iterator.Index();
        if(this->node_is_active(index) || this->node_is_dirichlet(index))
            for(int v=1;v<=d;v++)
                state.x(v)(index)=T(0);
    }
}

//#####################################################################
// Function Initialize_Muscles
//#####################################################################
template<class T,int d,bool enable_constraints,bool enable_muscles> void SKINNING_NONLINEAR_ELASTICITY<T,d,enable_constraints,enable_muscles>::
Initialize_Muscles()
{
    LOG::SCOPE scope("SKINNING_NONLINEAR_ELASTICITY::Initialize_Muscles()");

    BASE::cell_muscles.Resize(padded_cell_domain,true,false);
    BASE::cell_fibers.Resize(padded_cell_domain,true,false);
    BASE::cell_F_fibers.Resize(padded_cell_domain,true,false);
    BASE::cell_densities.Resize(padded_cell_domain,true,false);
    BASE::cell_c1.Resize(padded_cell_domain,true,false);
    BASE::cell_c2.Resize(padded_cell_domain,true,false);
}
//#####################################################################
// Function Tension/Tension_Derivative
//#####################################################################
template<class T,int d,bool enable_constraints,bool enable_muscles> T SKINNING_NONLINEAR_ELASTICITY<T,d,enable_constraints,enable_muscles>::
Tension(const T stretch, const T activation, const T density, const T fiber_max_stress)
{
    T fiber_p1=(T).05;
    T fiber_p2=(T)6.6;
    T fiber_cutoff=(T)1.4;
    T cutoff_scaled=fiber_p2*(fiber_cutoff-1);
    T fiber_p3=fiber_p1*fiber_p2*(exp(cutoff_scaled)-1);
    T fiber_p4=fiber_p1*(exp(cutoff_scaled)*(1-fiber_p2*fiber_cutoff)+fiber_p2-1);

    T strain=stretch-1,strain_abs=abs(strain),active_tension=0,passive_tension=0,scale=(T)25/(T)6;
    if(stretch>fiber_cutoff)passive_tension=fiber_p3*stretch+fiber_p4;else if(stretch>1)passive_tension=fiber_p1*(exp(fiber_p2*strain)-fiber_p2*strain-1);
    if(strain_abs<.4)active_tension=activation*density*(1-scale*sqr(strain));else if(strain_abs<.6)active_tension=2*scale*activation*density*sqr(strain_abs-(T).6);
    return fiber_max_stress*(active_tension+passive_tension);
}      
template<class T,int d,bool enable_constraints,bool enable_muscles> T SKINNING_NONLINEAR_ELASTICITY<T,d,enable_constraints,enable_muscles>::
Tension_Derivative(const T stretch, const T activation, const T density, const T fiber_max_stress)
{
    T fiber_p1=(T).05;
    T fiber_p2=(T)6.6;
    T fiber_cutoff=(T)1.4;
    T cutoff_scaled=fiber_p2*(fiber_cutoff-1);
    T fiber_p3=fiber_p1*fiber_p2*(exp(cutoff_scaled)-1);

    T strain=stretch-1,strain_abs=abs(strain),active_tension_derivative=0,passive_tension_derivative=0,scale=(T)25/(T)6;
    if(stretch>fiber_cutoff)passive_tension_derivative=fiber_p3;else if(stretch>1)passive_tension_derivative=fiber_p1*fiber_p2*(exp(fiber_p2*strain)-1);
    if(strain_abs<.4)active_tension_derivative=-2*scale*activation*density*strain;else if(strain_abs<.6)active_tension_derivative=4*scale*activation*density*(strain-sign(strain)*(T).6);
    return fiber_max_stress*(active_tension_derivative+passive_tension_derivative);
}
//#####################################################################
// Function Add_Embedded_Point_To_Fixed_Point_Spring_Constraint
//#####################################################################
template<class T,int d,bool enable_constraints,bool enable_muscles> int SKINNING_NONLINEAR_ELASTICITY<T,d,enable_constraints,enable_muscles>::
Add_Embedded_Point_To_Fixed_Point_Spring_Constraint(const T spring_coefficient,const TV& embedded_point_material_space_location,const TV& fixed_point_world_space_location, bool is_static)
{
    ARRAY<CONSTRAINT_SEGMENT<T,d> >* point_constraints=NULL;
    if( is_static ){
        point_constraints = &static_point_constraints;
    }
    else{
        point_constraints = &dynamic_point_constraints;
    }

    T_INDEX cell_index=grid.Cell(embedded_point_material_space_location,0);
    PHYSBAM_ASSERT(unpadded_cell_domain.Lazy_Inside(cell_index));

    if(this->cell_type(cell_index)!=INTERIOR_CELL_TYPE && this->cell_type(cell_index)!=DIRICHLET_CELL_TYPE && this->cell_type(cell_index)!=BOUNDARY_CELL_TYPE){

        int new_constraint_index=(*point_constraints).Append(CONSTRAINT_SEGMENT<T,d>());
        CONSTRAINT_SEGMENT<T,d>& new_constraint=(*point_constraints)(new_constraint_index);
        new_constraint.is_reference = false;
        new_constraint.endpoints[0].type=CONSTRAINT_NODE<T,d>::KINEMATIC;
        new_constraint.endpoints[1].type=CONSTRAINT_NODE<T,d>::KINEMATIC;
        new_constraint.endpoints[0].spatial_location()=fixed_point_world_space_location;
        new_constraint.endpoints[1].spatial_location()=fixed_point_world_space_location;
        return new_constraint_index;
    }


    TV multilinear_coordinates=(embedded_point_material_space_location-grid.Node(cell_index))/h;
    PHYSBAM_ASSERT(multilinear_coordinates.Min()>-1e-4 && multilinear_coordinates.Max()<1+1e-4);


    int new_constraint_index=(*point_constraints).Append(CONSTRAINT_SEGMENT<T,d>());
    CONSTRAINT_SEGMENT<T,d>& new_constraint=(*point_constraints)(new_constraint_index);

    new_constraint.is_reference = false;
    new_constraint.endpoints[0].type=CONSTRAINT_NODE<T,d>::KINEMATIC;
    new_constraint.endpoints[0].spatial_location()=fixed_point_world_space_location;   
    new_constraint.endpoints[1].type=CONSTRAINT_NODE<T,d>::GRID_FIXED;
    new_constraint.endpoints[1].grid_index()=cell_index;
    new_constraint.endpoints[1].multilinear_coordinates() = multilinear_coordinates;
    new_constraint.spring_coefficient = spring_coefficient;   

    return new_constraint_index;
}
//#####################################################################
// Function Add_Two_Embedded_Point_Spring_Constraint
//#####################################################################
template<class T,int d,bool enable_constraints,bool enable_muscles> int SKINNING_NONLINEAR_ELASTICITY<T,d,enable_constraints,enable_muscles>::
Add_Two_Embedded_Point_Spring_Constraint(const T spring_coefficient,const TV& embedded_point_material_space_location1,const TV& embedded_point_material_space_location2, bool is_static)
{
    ARRAY<CONSTRAINT_SEGMENT<T,d> >* point_constraints=NULL;

    if( is_static ){
        point_constraints = &static_point_constraints;
    }
    else{
        point_constraints = &dynamic_point_constraints;
    }

    T_INDEX cell_index1=grid.Cell(embedded_point_material_space_location1,0);
    PHYSBAM_ASSERT(unpadded_cell_domain.Lazy_Inside(cell_index1));
    PHYSBAM_ASSERT(this->cell_type(cell_index1)==INTERIOR_CELL_TYPE || this->cell_type(cell_index1)==DIRICHLET_CELL_TYPE || this->cell_type(cell_index1)==BOUNDARY_CELL_TYPE);
    TV multilinear_coordinates1=(embedded_point_material_space_location1-grid.Node(cell_index1))/h;
    PHYSBAM_ASSERT(multilinear_coordinates1.Min()>-1e-4 && multilinear_coordinates1.Max()<1+1e-4);

    T_INDEX cell_index2=grid.Cell(embedded_point_material_space_location2,0);
    PHYSBAM_ASSERT(unpadded_cell_domain.Lazy_Inside(cell_index2));
    PHYSBAM_ASSERT(this->cell_type(cell_index2)==INTERIOR_CELL_TYPE || this->cell_type(cell_index2)==DIRICHLET_CELL_TYPE || this->cell_type(cell_index2)==BOUNDARY_CELL_TYPE);
    TV multilinear_coordinates2=(embedded_point_material_space_location2-grid.Node(cell_index2))/h;
    PHYSBAM_ASSERT(multilinear_coordinates2.Min()>-1e-4 && multilinear_coordinates2.Max()<1+1e-4);

    int new_constraint_index=(*point_constraints).Append(CONSTRAINT_SEGMENT<T,d>());
    CONSTRAINT_SEGMENT<T,d>& new_constraint=(*point_constraints)(new_constraint_index);

    new_constraint.is_reference = false;     
    new_constraint.endpoints[0].type=CONSTRAINT_NODE<T,d>::GRID_FIXED;
    new_constraint.endpoints[0].grid_index()=cell_index1;
    new_constraint.endpoints[0].multilinear_coordinates() = multilinear_coordinates1;
    new_constraint.endpoints[1].type=CONSTRAINT_NODE<T,d>::GRID_FIXED;
    new_constraint.endpoints[1].grid_index()=cell_index2;
    new_constraint.endpoints[1].multilinear_coordinates() = multilinear_coordinates2;
    new_constraint.spring_coefficient = spring_coefficient;
  
    return new_constraint_index;
}

//#####################################################################
// Function Multilinear_Interpolation_Stencil
//#####################################################################
template<class T,int d,bool enable_constraints,bool enable_muscles> typename SKINNING_NONLINEAR_ELASTICITY<T,d,enable_constraints,enable_muscles>::T_STENCIL
SKINNING_NONLINEAR_ELASTICITY<T,d,enable_constraints,enable_muscles>::
Multilinear_Interpolation_Stencil(const T_INDEX& cell_index,const TV& multilinear_coordinates)
{
    T_STENCIL interpolation_stencil;
    for(RANGE_ITERATOR<d> iterator(RANGE<T_INDEX>(cell_index,cell_index+1));iterator.Valid();iterator.Next()){
        T_INDEX node_index=iterator.Index();
        T weight=(T)1.;
        for(int v=1;v<=d;v++) weight*=node_index(v)==cell_index(v)?(T)1.-multilinear_coordinates(v):multilinear_coordinates(v);
        interpolation_stencil.Insert(node_index,weight);}
    return interpolation_stencil;
}
//#####################################################################
// Function Deformation
//#####################################################################
template<class T,int d,bool enable_constraints,bool enable_muscles> typename SKINNING_NONLINEAR_ELASTICITY<T,d,enable_constraints,enable_muscles>::TV
SKINNING_NONLINEAR_ELASTICITY<T,d,enable_constraints,enable_muscles>::
Deformation(const T_INDEX& cell_index,const TV& multilinear_coordinates, const T_STATE& state) const
{
   return Deformation_Grid(cell_index,multilinear_coordinates, BASE::View_Convert(state.x));
}
//#####################################################################
// Function Deformation
//#####################################################################
template<class T,int d,bool enable_constraints,bool enable_muscles> typename SKINNING_NONLINEAR_ELASTICITY<T,d,enable_constraints,enable_muscles>::TV
SKINNING_NONLINEAR_ELASTICITY<T,d,enable_constraints,enable_muscles>::
Deformation_Grid(const T_INDEX& cell_index,const TV& multilinear_coordinates, T_VECTOR_VARIABLE_VIEW_CONST u) const
{
#if 0
    T_STENCIL interpolation_stencil=Multilinear_Interpolation_Stencil(cell_index,multilinear_coordinates);
    TV result;
    for(int v=1;v<=d;v++) result(v)=interpolation_stencil*u(v);
    return result+grid.Node(cell_index)+grid.dX*multilinear_coordinates;
#else
    return Displacement_Grid(cell_index, multilinear_coordinates, u) + grid.dX*multilinear_coordinates;
#endif
}
//#####################################################################
// Function Deformation
//#####################################################################
template<class T,int d,bool enable_constraints,bool enable_muscles> typename SKINNING_NONLINEAR_ELASTICITY<T,d,enable_constraints,enable_muscles>::TV
SKINNING_NONLINEAR_ELASTICITY<T,d,enable_constraints,enable_muscles>::
Displacement_Grid(const T_INDEX& cell_index,const TV& multilinear_coordinates, T_VECTOR_VARIABLE_VIEW_CONST du)
{
    T_STENCIL interpolation_stencil=Multilinear_Interpolation_Stencil(cell_index,multilinear_coordinates);
    TV result;
    for(int v=1;v<=d;v++) result(v)=interpolation_stencil*du(v);
    return result;
}
//#####################################################################
template class SKINNING_NONLINEAR_ELASTICITY<float,2,true,true>;
template class SKINNING_NONLINEAR_ELASTICITY<float,2,true,false>;
template class SKINNING_NONLINEAR_ELASTICITY<float,3,true,true>;
template class SKINNING_NONLINEAR_ELASTICITY<float,3,true,false>;
#ifndef USE_SPECIALIZED_KERNELS
template class SKINNING_NONLINEAR_ELASTICITY<double,2,true,true>;
template class SKINNING_NONLINEAR_ELASTICITY<double,2,true,false>;
template class SKINNING_NONLINEAR_ELASTICITY<double,3,true,true>;
template class SKINNING_NONLINEAR_ELASTICITY<double,3,true,false>;
#endif

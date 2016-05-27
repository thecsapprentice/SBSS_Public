#include "SKINNING_NONLINEAR_ELASTICITY.h"
#include <PhysBAM_Tools/Math_Tools/cube.h>
#include <PhysBAM_Tools/Matrices/SYMMETRIC_MATRIX_2X2.h>
#include <PhysBAM_Tools/Matrices/SYMMETRIC_MATRIX_3X3.h>
//#include <PhysBAM_Tools/Read_Write/Grids_Uniform_Arrays/READ_WRITE_ARRAYS.h>
//#include <PhysBAM_Tools/Read_Write/Matrices_And_Vectors/READ_WRITE_MATRIX.h>
//#include <PhysBAM_Tools/Read_Write/Arrays/READ_WRITE_ARRAY.h>
//#include <PhysBAM_Tools/Read_Write/Vectors/READ_WRITE_VECTOR_3D.h>
#include <Common/RANGE_ITERATOR.h>
#include <Common/STENCIL_ITERATOR.h>

#include <iostream>
#include <fstream>

using namespace PhysBAM;


//#####################################################################
// Function Update_Position_Based_State_Constraints
//#####################################################################

#ifndef USE_SPECIALIZED_KERNELS
#include <SIMD_Optimized_Kernels/Kernel_Wrappers/Muscle_Tension/Muscle_Tension.h>
#include <SIMD_Optimized_Kernels/Kernel_Wrappers/Central_Gradient/Central_Gradient.h>
#endif

template<class T,int d,bool enable_constraints,bool enable_muscles> void SKINNING_NONLINEAR_ELASTICITY<T,d,enable_constraints,enable_muscles>::
Update_Position_Based_State_Constraints(T_VECTOR_VARIABLE_VIEW_CONST u,T_SCALAR_VARIABLE_VIEW_CONST p)
{
    if(!enable_constraints) return;
    LOG::SCOPE scope("SKINNING_NONLINEAR_ELASTICITY::Update_Position_Based_State_Constraints()");

    // Constraints
    //Build_Constraint_Matrix();

    //Update_Collision_Constraints(u);

}
//#####################################################################
// Function Add_Force_Constraints
//#####################################################################
template<class T,int d,bool enable_constraints,bool enable_muscles> void SKINNING_NONLINEAR_ELASTICITY<T,d,enable_constraints,enable_muscles>::
Add_Force_Constraints(const ARRAY<CONSTRAINT_SEGMENT<T,d> >& constraints,T_VECTOR_VARIABLE_VIEW_CONST u,T_SCALAR_VARIABLE_VIEW_CONST p, T_VECTOR_VARIABLE_VIEW f,T_SCALAR_VARIABLE_VIEW q) const
{
    if(!enable_constraints) return;
    LOG::SCOPE scope("SKINNING_NONLINEAR_ELASTICITY::Add_Force_Constraints()");

    VECTOR<TV,2> X, G;  // G is the vector of endpoint forces
//    VECTOR<TV,2> computed_velocity;
        
        for(int c=1;c<=constraints.m;c++){

            //LOG::cout << "Constraint " << c << std::endl;
            
            // PHASE 1
            for(int e=1; e<=2; e++){
                    const CONSTRAINT_NODE<T,d>& endpoint = constraints(c).endpoints[e-1];
                    //LOG::cout << "Endpoint " << e << " is ";
                    switch(endpoint.type){
                    case CONSTRAINT_NODE<T,d>::GRID_FIXED:
                        X(e) = Deformation_Grid(endpoint.grid_index(),endpoint.multilinear_coordinates(), u);
                        //LOG::cout << "grid " << p << std::endl;
                        break;
                    case CONSTRAINT_NODE<T,d>::MESH_FIXED:
                        PHYSBAM_FATAL_ERROR("No mesh allowed in traditional skinning.");
                        //p = Deformation_Mesh(endpoint.mesh_index(), endpoint.multilinear_coordinates);
                        //LOG::cout << "mesh " << p << std::endl;
                        //break;
                    case CONSTRAINT_NODE<T,d>::KINEMATIC:
                        X(e) = endpoint.spatial_location();
                        //LOG::cout << "kinematic " << p << std::endl;
                        break;
                    case CONSTRAINT_NODE<T,d>::KINEMATIC_PTR:
                        X(e) = endpoint.spatial_location_ptr(collision_spring_locations);
                        break;
                    }
            }
            // PHASE 2
            TV force = X(1) - X(2);
            if( constraints(c).is_reference )
                force *= collision_spring_constants(constraints(c).spring_coefficient_ptr);
            else
                force *= constraints(c).spring_coefficient;
            G(1) = -force;  //Each endpoint undergoes equal and opposite force
            G(2) = force;
            
            //LOG::cout << "Force: " << force << std::endl;

            // PHASE 3
            for(int e=1; e<=2; e++){
                    const CONSTRAINT_NODE<T,d>& endpoint = constraints(c).endpoints[e-1];
                    switch(endpoint.type){                       
                    case CONSTRAINT_NODE<T,d>::GRID_FIXED:
                        {
                            const T_STENCIL interpolation_stencil=Multilinear_Interpolation_Stencil(endpoint.grid_index(),endpoint.multilinear_coordinates());
                            for(T_CONST_STENCIL_ITERATOR stencil_iterator(interpolation_stencil);stencil_iterator.Valid();stencil_iterator.Next())
                                for(int v=1;v<=d;v++)
                                    f(v)(stencil_iterator.Key())+=G(e)(v)*stencil_iterator.Data();}
                        break;
                    
                    case CONSTRAINT_NODE<T,d>::MESH_FIXED:
                        {
                            PHYSBAM_FATAL_ERROR("No mesh allowed in traditional skinning.");
                            /*
                              const T_INDEX& grid_index = BASE::mesh_to_cell_map( endpoint.mesh_index() );
                              const T_STENCIL interpolation_stencil=Multilinear_Interpolation_Stencil(endpoint.mesh_index(),endpoint.multilinear_coordinates);
                              int flat_index=1;
                              for(T_CONST_STENCIL_ITERATOR stencil_iterator(interpolation_stencil);stencil_iterator.Valid();stencil_iterator.Next(),flat_index++){
                              PHYSBAM_ASSERT( BASE::mesh_to_node_rmap(stencil_iterator.Key() + grid_index).m > 0 ); // This should always be true
                              int mesh_node = BASE::elements(endpoint.mesh_index())(flat_index);
                              if( BASE::is_mesh_node_shared(mesh_node) >= 0){
                              //LOG::cout << "Applying Mesh force to grid node " << stencil_iterator.Key() << std::endl;
                              for(int v=1;v<=d;v++)  f(v)(stencil_iterator.Key()+grid_index)+=endpoint.computed_force(v)*stencil_iterator.Data();
                              }
                              else{
                              //LOG::cout << "Applying Mesh force to mesh node " << mesh_node << std::endl;
                              for(int v=1;v<=d;v++)  f_mesh(v)(mesh_node)+=endpoint.computed_force(v)*stencil_iterator.Data();                            
                              }
                              }
                            */
                        }
                        break;
                    
                    case CONSTRAINT_NODE<T,d>::KINEMATIC:
                        // Do nothing. Free space locations do not receive forces.
                        break;
                    case CONSTRAINT_NODE<T,d>::KINEMATIC_PTR:
                        // Do nothing. Free space locations do not receive forces.
                        break;
                    }
            }
        }
    

}
//#####################################################################
// Function Add_Force_Differential_Constraints
//#####################################################################
template<class T,int d,bool enable_constraints,bool enable_muscles> void SKINNING_NONLINEAR_ELASTICITY<T,d,enable_constraints,enable_muscles>::
Add_Force_Differential_Constraints(const ARRAY<CONSTRAINT_SEGMENT<T,d> >& constraints,T_VECTOR_VARIABLE_VIEW_CONST du,T_SCALAR_VARIABLE_VIEW_CONST dp,T_VECTOR_VARIABLE_VIEW df,T_SCALAR_VARIABLE_VIEW dq) const
{
    if(!enable_constraints) return;
    // LOG::SCOPE scope("SKINNING_NONLINEAR_ELASTICITY::Add_Force_Differential_Constraints()");

    VECTOR<TV,2> dX,dG; // dG is the vector of endpoint forces
        
    for(int c=1;c<=constraints.m;c++){
            
        // PHASE 1 - Interpolate or collect displacements
        for(int e=1; e<=2; e++){
            const CONSTRAINT_NODE<T,d>& endpoint = constraints(c).endpoints[e-1];
            switch(endpoint.type){
                case CONSTRAINT_NODE<T,d>::GRID_FIXED:
                    dX(e) = Displacement_Grid(endpoint.grid_index(),endpoint.multilinear_coordinates(),du);
                    break;
                case CONSTRAINT_NODE<T,d>::KINEMATIC:
                case CONSTRAINT_NODE<T,d>::KINEMATIC_PTR:
                    dX(e) = TV(); // This is differential, fixed_points evaluate to zero
                    break;
                case CONSTRAINT_NODE<T,d>::MESH_FIXED:
                    PHYSBAM_FATAL_ERROR("Mesh embeddings not supported");
            }}

        // PHASE 2
        TV dforce=dX(1)-dX(2);
        if(constraints(c).is_reference)
            dforce *= collision_spring_constants(constraints(c).spring_coefficient_ptr);
        else
            dforce *= constraints(c).spring_coefficient;
        dG(1) = -dforce;  //Each endpoint undergoes equal and opposite force
        dG(2) = dforce;
            
        // PHASE 3
        for(int e=1; e<=2; e++){
            const CONSTRAINT_NODE<T,d>& endpoint = constraints(c).endpoints[e-1];
            switch(endpoint.type){
                case CONSTRAINT_NODE<T,d>::GRID_FIXED:
                    {
                        const T_STENCIL interpolation_stencil=Multilinear_Interpolation_Stencil(endpoint.grid_index(),endpoint.multilinear_coordinates());
                        for(T_CONST_STENCIL_ITERATOR stencil_iterator(interpolation_stencil);stencil_iterator.Valid();stencil_iterator.Next())
                            for(int v=1;v<=d;v++)
                                df(v)(stencil_iterator.Key())+=dG(e)(v)*stencil_iterator.Data();
                    }
                    break;
                case CONSTRAINT_NODE<T,d>::KINEMATIC:
                case CONSTRAINT_NODE<T,d>::KINEMATIC_PTR:
                    // Do nothing. Free space locations do not receive forces.
                    break;                    
                case CONSTRAINT_NODE<T,d>::MESH_FIXED:
                    PHYSBAM_FATAL_ERROR("Mesh embeddings not supported");
                // {
                // const T_INDEX& grid_index = BASE::mesh_to_cell_map( endpoint.mesh_index() );
                // const T_STENCIL interpolation_stencil=Multilinear_Interpolation_Stencil(endpoint.mesh_index(),endpoint.multilinear_coordinates);
                // int flat_index=1;
                // for(T_CONST_STENCIL_ITERATOR stencil_iterator(interpolation_stencil);stencil_iterator.Valid();stencil_iterator.Next(),flat_index++){
                // PHYSBAM_ASSERT( BASE::mesh_to_node_rmap(stencil_iterator.Key() + grid_index).m > 0 ); // This should always be true
                // int mesh_node = BASE::elements(endpoint.mesh_index())(flat_index);
                // if( BASE::is_mesh_node_shared(mesh_node) >= 0)
                // for(int v=1;v<=d;v++)  df(v)(stencil_iterator.Key()+grid_index)+=endpoint.dG(v)*stencil_iterator.Data();
                // else
                // for(int v=1;v<=d;v++)  df_mesh(v)(mesh_node)+=endpoint.dG(v)*stencil_iterator.Data();                            
                // }}
                    break;
            }}
    }
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

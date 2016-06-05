#include "HYBRID_NONLINEAR_ELASTICITY.h"
#include <Common/STENCIL_ITERATOR.h>
#include <Common/RANGE_ITERATOR.h>
#include <EngineInterface/CONSTRAINTS.h>
#include <Common/GENERIC_CELL.h>
using namespace PhysBAM;

#ifdef LOG_DETAILED_PERFORMANCE
#define LOG_DETAILED_PERFORMANCE_NE
#endif

//#####################################################################
// Function Initialize_Blocks_Constraints
//#####################################################################
template<class T_NONLINEAR_ELASTICITY> void HYBRID_NONLINEAR_ELASTICITY<T_NONLINEAR_ELASTICITY>::
Initialize_Blocks_Constraints(int number_of_partitions)
{
#ifdef USE_SPECIALIZED_KERNELS
    typedef TRIPLE< GENERIC_CELL<T,d>, CONSTRAINT_SEGMENT<T,d>, int > CELL_CONSTRAINT_TEMP;

    ARRAY< int, T_INDEX > c_counter;
    ARRAY< ARRAY< TRIPLE< GENERIC_CELL<T,d>, CONSTRAINT_SEGMENT<T,d>, int > >, T_INDEX > cell_constraints;
    cell_constraints.Resize(padded_cell_domain);
    c_counter.Resize(padded_cell_domain); c_counter.Fill(0);
    int block_constraint_offset = 1;

    reference_constraint_block_ptr.Resize( BASE::collision_constraints.m ); // TODO: This might need to change later.

    ARRAY<int> constraint_base_offsets;
    ARRAY<int> constraint_base_lengths;
    ARRAY<BLOCKED_TYPE<TV,8> > constraint_multilinear_coordinates;
    ARRAY<BLOCKED_TYPE<T,8> > constraint_spring_constants;

    ARRAY<BLOCKED_TYPE<TV,8> > constraint_node_positions;
    ARRAY<BLOCKED_TYPE<int,8> > spring_id;
    ARRAY<BLOCKED_TYPE<int,8> > spring_id_X;
    ARRAY<BLOCKED_TYPE<int,8> > spring_id_Y;
    ARRAY<BLOCKED_TYPE<int,8> > spring_id_Z;

    constraint_base_offsets.Resize(number_of_blocks);
    constraint_base_lengths.Resize(number_of_blocks);

    int collected_constraints = 0;
    for(int c=1;c<=BASE::collision_constraints.m;c++){
        const CONSTRAINT_SEGMENT<T,d>& cs = BASE::collision_constraints(c);
        
        if( cs.isSinglePoint() ){
            T_INDEX grid_index;
            int mesh_index;
            bool isMesh;
            CELL_CONSTRAINT_TEMP cct;

            const CONSTRAINT_NODE<T,d>&  embedded_point = cs.getEmbeddedPoint();
            if( embedded_point.type == CONSTRAINT_NODE<T,d>::GRID_FIXED ){
                grid_index = embedded_point.grid_index();
                isMesh = false;
                cct.x.is_grid=true;
                cct.x.grid_index() = grid_index;
            }
            else{
                PHYSBAM_ASSERT( (embedded_point.type == CONSTRAINT_NODE<T,d>::MESH_FIXED) );
                mesh_index = embedded_point.mesh_index();
                grid_index = cell_indices_mesh(mesh_index);
                isMesh=true;
                cct.x.is_grid=false;
                cct.x.mesh_index() = mesh_index;
            }
            cct.y = cs;
            cct.z = c;
            collected_constraints++;
            cell_constraints(grid_index).Append( cct ); 
            c_counter(grid_index)+=1;
            PHYSBAM_ASSERT( cell_constraints(grid_index).Last().y.endpoints[1].grid_index() == cs.endpoints[1].grid_index() );
        }
        else if( cs.isDualPoint() ) {
            // Maybe do something here!
            PHYSBAM_FATAL_ERROR("Cannot add Two point Embedded constraints to the static constraint list.");
        }
        else{
            PHYSBAM_FATAL_ERROR("Other Constraints are not supported." );
        }
    }
           
    LOG::cout << "Total collision constraints: " << collected_constraints << std::endl;
    PHYSBAM_ASSERT( collected_constraints == BASE::collision_constraints.m );
    int total_blocked_constraints = 0;

    for( int block = 1; block <= number_of_blocks; block++ ){
        const T_INDEX& base_index = block_bases(block);
       
        //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        //
        //              Setup Blocked Constraint Data
        //
        //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        
        int block_constraint_count = 0;
        
        int flat_index=1;
        VECTOR<int,8> block_cell_bins;
        VECTOR<ARRAY<TRIPLE< GENERIC_CELL<T,d>, CONSTRAINT_SEGMENT<T,d>, int > >, 8> constraint_block_bins;
        
        for(RANGE_ITERATOR<d> cell_iterator(RANGE<T_INDEX>(base_index,base_index+1));
            cell_iterator.Valid();
            cell_iterator.Next(),flat_index++){
            const T_INDEX& index = cell_iterator.Index();
            int is_mesh = hybrid_block_cells(block)(flat_index);

            for(int m = 1; m <= cell_constraints(index).m; m++){
                switch(is_mesh){
                case -1:
                    break;
                case 0:
                    if(cell_constraints(index)(m).x.is_grid){
                        block_cell_bins(flat_index) += 1;
                        constraint_block_bins(flat_index).Append(cell_constraints(index)(m));
                        c_counter(index)--;
                        //LOG::cout << cell_constraints(index)(m).y << std::endl;
                        //LOG::cout << "Grid Constraint Found in cell " << index << " within block." << std::endl;
                        total_blocked_constraints++;
                    }
                    break;
                default:
                    if(!(cell_constraints(index)(m).x.is_grid) && cell_constraints(index)(m).x.mesh_index() == is_mesh){
                        block_cell_bins(flat_index) += 1;                    
                        constraint_block_bins(flat_index).Append(cell_constraints(index)(m));
                        c_counter(index)--;
                        //LOG::cout << cell_constraints(index)(m).y << std::endl;    
                        //LOG::cout << "Mesh Constraint Found in cell " << is_mesh << " within block." << std::endl;
                        total_blocked_constraints++;
                    }
                    break;
                }
            }
                       
        }
                    
        block_constraint_count = block_cell_bins.Max();

        constraint_base_offsets(block) = block_constraint_offset;
        constraint_base_lengths(block) = block_constraint_count;
        block_constraint_offset += block_constraint_count;

        for( int block_constraint = 0; block_constraint < block_constraint_count; block_constraint++){               
            constraint_multilinear_coordinates.Append(BLOCKED_TYPE<TV,8>());
            constraint_spring_constants.Append(BLOCKED_TYPE<T,8>());
            constraint_node_positions.Append(BLOCKED_TYPE<TV,8>());
            spring_id.Append(BLOCKED_TYPE<int,8>());
            spring_id_X.Append(BLOCKED_TYPE<int,8>());
            spring_id_Y.Append(BLOCKED_TYPE<int,8>());
            spring_id_Z.Append(BLOCKED_TYPE<int,8>());

            int cell=0;
            int flat_index=1;
            for(RANGE_ITERATOR<d> cell_iterator(RANGE<T_INDEX>(base_index,base_index+1));
                cell_iterator.Valid();
                cell_iterator.Next(),cell++,flat_index++)
                {                          
                    constraint_multilinear_coordinates.Last().Set(TV::All_Ones_Vector()*0.5, cell);
                    constraint_node_positions.Last().Set(TV(), cell);
                    constraint_spring_constants.Last().Set(T(0), cell);
                    spring_id.Last().Set(0, cell);
                    spring_id_X.Last().Set(0, cell);
                    spring_id_Y.Last().Set(0, cell);
                    spring_id_Z.Last().Set(0, cell);
                    const T_INDEX& index = cell_iterator.Index();
                    int is_mesh = hybrid_block_cells(block)(flat_index);

                    if(constraint_block_bins(flat_index).m > block_constraint)
                        {
                            switch(is_mesh){
                            case -1:
                                // Do Nothing. Already cleared out
                                break;
                            case 0:
                                if(constraint_block_bins(flat_index)(block_constraint+1).x.is_grid){
                                    constraint_multilinear_coordinates.Last().Set(constraint_block_bins(flat_index)(block_constraint+1).y.getEmbeddedPoint().multilinear_coordinates(), cell);
                                    constraint_node_positions.Last().Set( BASE::grid.Node(constraint_block_bins(flat_index)(block_constraint+1).x.grid_index()), cell );
                                    if(constraint_block_bins(flat_index)(block_constraint+1).y.is_reference){
                                        int spring_reference_loc = constraint_block_bins(flat_index)(block_constraint+1).y.spring_coefficient_ptr;
                                        spring_id.Last().Set(spring_reference_loc-1, cell);
                                        spring_id_X.Last().Set( ((spring_reference_loc-1)*3)+0, cell);
                                        spring_id_Y.Last().Set( ((spring_reference_loc-1)*3)+1, cell);
                                        spring_id_Z.Last().Set( ((spring_reference_loc-1)*3)+2, cell);
                                    }
                                    else{
                                        constraint_spring_constants.Last().Set(constraint_block_bins(flat_index)(block_constraint+1).y.spring_coefficient, cell);
                                        spring_id.Last().Set(0, cell);
                                        spring_id_X.Last().Set(0, cell);
                                        spring_id_Y.Last().Set(1, cell);
                                        spring_id_Z.Last().Set(2, cell);
                                    }
                                }
                                break;
                            default:
                                if(!(constraint_block_bins(flat_index)(block_constraint+1).x.is_grid) && constraint_block_bins(flat_index)(block_constraint+1).x.mesh_index() == is_mesh){
                                    PHYSBAM_ASSERT( (constraint_block_bins(flat_index)(block_constraint+1).y.getEmbeddedPoint().type == CONSTRAINT_NODE<T,d>::MESH_FIXED) );
                                    constraint_multilinear_coordinates.Last().Set(constraint_block_bins(flat_index)(block_constraint+1).y.getEmbeddedPoint().multilinear_coordinates(), cell);
                                    constraint_node_positions.Last().Set( BASE::grid.Node(cell_indices_mesh(constraint_block_bins(flat_index)(block_constraint+1).x.mesh_index())), cell );
                                    if(constraint_block_bins(flat_index)(block_constraint+1).y.is_reference){
                                        int spring_reference_loc = constraint_block_bins(flat_index)(block_constraint+1).y.spring_coefficient_ptr;
                                        spring_id.Last().Set(spring_reference_loc-1, cell);
                                        spring_id_X.Last().Set( ((spring_reference_loc-1)*3)+0, cell);
                                        spring_id_Y.Last().Set( ((spring_reference_loc-1)*3)+1, cell);
                                        spring_id_Z.Last().Set( ((spring_reference_loc-1)*3)+2, cell); 
                                    }
                                    else{
                                        constraint_spring_constants.Last().Set(constraint_block_bins(flat_index)(block_constraint+1).y.spring_coefficient, cell);
                                        spring_id.Last().Set(0, cell);
                                        spring_id_X.Last().Set(0, cell);
                                        spring_id_Y.Last().Set(1, cell);
                                        spring_id_Z.Last().Set(2, cell);
                                    }
                                }
                                break;
                            }
                        }

                    switch(is_mesh){
                    case -1:
                        // Do nothing. Already Cleared Out
                        break;
                    case 0:
                        if( !(cell_type(index)==INTERIOR_CELL_TYPE || (/*first_order && */cell_type(index)==BOUNDARY_CELL_TYPE)) )
                            {
                                constraint_multilinear_coordinates.Last().Set(TV::All_Ones_Vector()*0.5, cell);
                                constraint_spring_constants.Last().Set(T(0), cell);
                                constraint_node_positions.Last().Set(TV(), cell);
                                spring_id.Last().Set(0, cell);
                                spring_id_X.Last().Set(0, cell);
                                spring_id_Y.Last().Set(1, cell);
                                spring_id_Z.Last().Set(2, cell);
                            }
                        break;
                    default:
                        if( !(cell_type_mesh(is_mesh)==INTERIOR_CELL_TYPE || (/*first_order && */cell_type_mesh(is_mesh)==BOUNDARY_CELL_TYPE)) )
                            {
                                constraint_multilinear_coordinates.Last().Set(TV::All_Ones_Vector()*0.5, cell);
                                constraint_spring_constants.Last().Set(T(0), cell);
                                constraint_node_positions.Last().Set(TV(), cell);
                                spring_id.Last().Set(0, cell);
                                spring_id_X.Last().Set(0, cell);
                                spring_id_Y.Last().Set(1, cell);
                                spring_id_Z.Last().Set(2, cell);
                            }
                        break;
                    }                       
                }        
        }
    }

#if 0
    for(RANGE_ITERATOR<d> cell_iterator(unpadded_cell_domain);
        cell_iterator.Valid();
        cell_iterator.Next()){
        const T_INDEX& index = cell_iterator.Index();
        
        if( c_counter(index) ){
            LOG::cout << "Missed Constraints: " << c_counter(index) << std::endl;
            LOG::cout << "total_blocked_constraints: " << total_blocked_constraints << std::endl;
            LOG::cout << "index: " << index << std::endl;
            
            LOG::cout << "All constraints..." << std::endl;
            for( int i = 1; i <= cell_constraints(index).m; i++){
                LOG::cout << cell_constraints(index)(i).y << std::endl;    
            }               
            
            for( int block = 1; block <= number_of_blocks; block++ ){
                const T_INDEX& base_index = block_bases(block);
                RANGE<T_INDEX> block_range( base_index, base_index+1);
                if( block_range.Lazy_Inside( index ) ){
                    LOG::cout << "Block " << block << std::endl;
                    LOG::cout << hybrid_block_cells(block) << std::endl;
                    
                    //for( int j = 1; j <= 8; j ++){
                    //    LOG::cout << "Found constraints (" << j << ")..." << std::endl;
                    //    for( int i = 1; i <= constraint_block_bins(j).m; i++){
                    //        LOG::cout << constraint_block_bins(j)(i).y << std::endl;    
                    //    }
                    //}
                }
            }            
        }
        PHYSBAM_ASSERT( c_counter(index) == 0 );
    }
#endif


    ARRAY<T>& collision_spring_constants = BASE::collision_spring_constants;
    ARRAY<TV>& collision_spring_locations = BASE::collision_spring_locations;
    BASE::specialized_data.InitializeConstraintData(constraint_base_offsets, constraint_base_lengths, constraint_multilinear_coordinates, constraint_spring_constants, constraint_node_positions,
                                                    spring_id, spring_id_X, spring_id_Y, spring_id_Z,
                                                    &(collision_spring_constants),
                                                    &(collision_spring_locations));

    LOG::cout << "Total blocked constraints: " << total_blocked_constraints << std::endl;

    PHYSBAM_ASSERT( BASE::collision_constraints.m == total_blocked_constraints );
#endif
}
//#####################################################################
// Function Initialize_Blocks_Muscles
//#####################################################################
template<class T_NONLINEAR_ELASTICITY> void HYBRID_NONLINEAR_ELASTICITY<T_NONLINEAR_ELASTICITY>::
Initialize_Blocks_Muscles(int number_of_partitions)
{
    LOG::SCOPE scope("HYBRID_NONLINEAR_ELASTICITY::Initialize_Blocks_Muscles()");
    const int DEFAULT_MUSCLE = 1;

    int block_muscle_offset = 1;
    ARRAY<int> muscle_base_offsets;
    ARRAY<int> muscle_base_lengths;
    ARRAY<BLOCKED_TYPE<int,8> > muscle_id;
    ARRAY<BLOCKED_TYPE<T,8> > muscle_density;
    ARRAY<BLOCKED_TYPE<TV,8> > muscle_fiber;

    muscle_base_offsets.Resize(number_of_blocks);
    muscle_base_lengths.Resize(number_of_blocks);

    for( int block = 1; block <= number_of_blocks; block++ ){
        const T_INDEX& base_index = block_bases(block);
       
        //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        //
        //              Setup Blocked Muscle Data
        //
        //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        
        int block_muscle_count = 0;
        
        int flat_index=1;
        VECTOR<int,8> block_cell_bins;
        
        for(RANGE_ITERATOR<d> cell_iterator(RANGE<T_INDEX>(base_index,base_index+1));
            cell_iterator.Valid();
            cell_iterator.Next(),flat_index++){
            const T_INDEX& index = cell_iterator.Index();
            int is_mesh = hybrid_block_cells(block)(flat_index);

            switch(is_mesh){
            case -1:
                break;
            case 0:
                block_muscle_count = max(cell_muscles(index).Size(), block_muscle_count);
                break;
            default:
                block_muscle_count = max(cell_muscles_mesh(is_mesh).Size(), block_muscle_count);
                break;
            }
        }               
                          
        muscle_base_offsets(block) = block_muscle_offset;
        muscle_base_lengths(block) = block_muscle_count;
        block_muscle_offset += block_muscle_count;

        for( int block_muscle = 0; block_muscle < block_muscle_count; block_muscle++){  
            muscle_id.Append(BLOCKED_TYPE<int,8>());
            muscle_density.Append(BLOCKED_TYPE<T,8>());
            muscle_fiber.Append(BLOCKED_TYPE<TV,8>());
                        
            int cell=0;
            int flat_index=1;
            for(RANGE_ITERATOR<d> cell_iterator(RANGE<T_INDEX>(base_index,base_index+1));
                cell_iterator.Valid();
                cell_iterator.Next(),cell++,flat_index++)
                {
                    const T_INDEX& index = cell_iterator.Index();
                    int is_mesh = hybrid_block_cells(block)(flat_index);
                    
                    muscle_id.Last().Set(DEFAULT_MUSCLE-1,cell);
                    muscle_fiber.Last().Set(TV(), cell); 
                    muscle_density.Last().Set(T(), cell); 

                    switch(is_mesh){
                    case -1:
                        break;
                    case 0:
                        if(cell_muscles(index).m > block_muscle){
                            muscle_id.Last().Set(cell_muscles(index)(block_muscle+1)-1, cell);
                            muscle_fiber.Last().Set(cell_fibers(index)(block_muscle+1), cell); 
                            muscle_density.Last().Set(cell_densities(index)(block_muscle+1), cell); 
                        }
                        break;
                    default:
                        if(cell_muscles_mesh(is_mesh).m > block_muscle){
                            muscle_id.Last().Set(cell_muscles_mesh(is_mesh)(block_muscle+1)-1, cell);
                            muscle_fiber.Last().Set(cell_fibers_mesh(is_mesh)(block_muscle+1), cell); 
                            muscle_density.Last().Set(cell_densities_mesh(is_mesh)(block_muscle+1), cell); 
                        }
                        break;
                    }

                }
        }

    }

    BASE::specialized_data.InitializeMuscleData(muscle_base_offsets, muscle_base_lengths, muscle_id,
                                                muscle_fiber, muscle_density);
}
//#####################################################################
// Function Initialize_Muscles
//#####################################################################
template<class T_NONLINEAR_ELASTICITY> void HYBRID_NONLINEAR_ELASTICITY<T_NONLINEAR_ELASTICITY>::
Initialize_Muscles()
{
    LOG::SCOPE scope("HYBRID_NONLINEAR_ELASTICITY::Initialize_Muscles()");
    BASE::Initialize_Muscles();
    cell_muscles_mesh.Resize(Number_Of_Mesh_Cells());
    cell_fibers_mesh.Resize(Number_Of_Mesh_Cells());
    cell_F_fibers_mesh.Resize(Number_Of_Mesh_Cells());
    cell_densities_mesh.Resize(Number_Of_Mesh_Cells());
    cell_c1_mesh.Resize(Number_Of_Mesh_Cells());
    cell_c2_mesh.Resize(Number_Of_Mesh_Cells());
}

//#####################################################################
// Function Update_Position_Based_State_Muscles_And_Constraints
//#####################################################################
template<class T_NONLINEAR_ELASTICITY> void HYBRID_NONLINEAR_ELASTICITY<T_NONLINEAR_ELASTICITY>::
Update_Position_Based_State_Muscles_And_Constraints(const T_STATE& state_u)
{
    if(Supports_Constraints())
        Update_Position_Based_State_Constraints_Mesh(BASE::View_Convert(state_u.x),state_u.p,View_Convert(state_u.x_mesh),state_u.p_mesh);
#ifndef USE_SPECIALIZED_KERNELS
    if(Supports_Muscles())
        Update_Position_Based_State_Muscles_Mesh(BASE::View_Convert(state_u.x),state_u.p,View_Convert(state_u.x_mesh),state_u.p_mesh);
#endif
}
//#####################################################################
// Function Update_Position_Based_State_Constraints_Mesh
//#####################################################################
template<class T_NONLINEAR_ELASTICITY> void HYBRID_NONLINEAR_ELASTICITY<T_NONLINEAR_ELASTICITY>::
Update_Position_Based_State_Constraints_Mesh(T_VECTOR_VARIABLE_VIEW_CONST u,T_SCALAR_VARIABLE_VIEW_CONST p,T_VECTOR_VARIABLE_MESH_VIEW_CONST u_mesh,T_SCALAR_VARIABLE_MESH_VIEW_CONST p_mesh)
{
    LOG::SCOPE scope("HYBRID_NONLINEAR_ELASTICITY::Update_Position_Based_State_Constraints_Mesh()");
    //Update_Collision_Constraints(u, u_mesh);
}
//#####################################################################
// Function Update_Position_Based_State_Muscles_Mesh
//#####################################################################
template<class T_NONLINEAR_ELASTICITY> void HYBRID_NONLINEAR_ELASTICITY<T_NONLINEAR_ELASTICITY>::
Update_Position_Based_State_Muscles_Mesh(T_VECTOR_VARIABLE_VIEW_CONST u,T_SCALAR_VARIABLE_VIEW_CONST p,T_VECTOR_VARIABLE_MESH_VIEW_CONST u_mesh,T_SCALAR_VARIABLE_MESH_VIEW_CONST p_mesh)
{
    LOG::SCOPE scope("HYBRID_NONLINEAR_ELASTICITY::Update_Position_Based_State_Muscles_Mesh()");
    // TODO: Need to think about this more as well. Muscles are always problematic
    PHYSBAM_NOT_IMPLEMENTED();
}  
/*
//#####################################################################
// Function Update_Collision_Constraints
//#####################################################################
template<class T_NONLINEAR_ELASTICITY> void HYBRID_NONLINEAR_ELASTICITY<T_NONLINEAR_ELASTICITY>::
Update_Collision_Constraints(T_VECTOR_VARIABLE_VIEW_CONST u, T_VECTOR_VARIABLE_MESH_VIEW_CONST u_mesh)
{
    if( !BASE::collision_shape )
        return;

    LOG::SCOPE scope("HYBRID_NONLINEAR_ELASTICITY:: Update_Collision_Constraints()");

   
    TV collision_center;
    collision_center(2) = -0.055f;
    //T collision_radius = 0.05f;
    TV collision_normal;
    collision_normal(2) = 1.0f; 

    for( int m=1; m<=BASE::collision_constraints.m; m++){
        CONSTRAINT_SEGMENT<T,d>& ic = BASE::collision_constraints(m);
        TV p;

        if(ic.endpoints[1].type == CONSTRAINT_NODE<T,d>::GRID_FIXED)
            p = Deformation_Grid(ic.endpoints[1].grid_index(),ic.endpoints[1].multilinear_coordinates(), u);
        else if(ic.endpoints[1].type == CONSTRAINT_NODE<T,d>::MESH_FIXED)
            p = Deformation_Mesh(ic.endpoints[1].mesh_index(),ic.endpoints[1].multilinear_coordinates(), u, u_mesh);
        else
            PHYSBAM_FATAL_ERROR();

        T depth = BASE::collision_shape->Phi( p );

        if( depth < 0 ) { // We are inside the collision shape!
            //LOG::cout << "COLLISION DETECTED!!!" << std::endl;
            //LOG::cout << "Point of Collision: " << p << std::endl;
            //LOG::cout << "New Constraint Point: " << BASE::collision_shape->ClosestPoint(p) << std::endl;
#if 1
            BASE::collision_spring_locations(ic.spring_coefficient_ptr) = BASE::collision_shape->ClosestPoint(p);
            BASE::collision_spring_constants(ic.spring_coefficient_ptr) = 1e4;
#else
            BASE::collision_spring_locations(ic.spring_coefficient_ptr) = TV();
            BASE::collision_spring_constants(ic.spring_coefficient_ptr) = 0.0;
#endif
        }
        else{
            BASE::collision_spring_locations(ic.spring_coefficient_ptr) = TV();
            BASE::collision_spring_constants(ic.spring_coefficient_ptr) = 0.0;
        }
    }
}
*/
//#####################################################################
// Function Add_Force_Muscles_And_Constraints
//#####################################################################
template<class T_NONLINEAR_ELASTICITY> void HYBRID_NONLINEAR_ELASTICITY<T_NONLINEAR_ELASTICITY>::
Add_Force_Muscles_And_Constraints(const T_STATE& state_u, T_STATE& state_f) const
{    
#ifndef USE_SPECIALIZED_KERNELS
    if(Supports_Muscles())
        Add_Force_Muscles_Mesh(BASE::View_Convert(state_u.x),state_u.p,View_Convert(state_u.x_mesh),state_u.p_mesh,state_f.x,state_f.p,state_f.x_mesh,state_f.p_mesh);
    if(Supports_Constraints()){
        Add_Force_Constraints_Mesh(BASE::dynamic_point_constraints,BASE::View_Convert(state_u.x),state_u.p,View_Convert(state_u.x_mesh),state_u.p_mesh,state_f.x,state_f.p,state_f.x_mesh,state_f.p_mesh);
        Add_Force_Constraints_Mesh(BASE::static_point_constraints,BASE::View_Convert(state_u.x),state_u.p,View_Convert(state_u.x_mesh),state_u.p_mesh,state_f.x,state_f.p,state_f.x_mesh,state_f.p_mesh);
        Add_Force_Constraints_Mesh(BASE::collision_constraints,BASE::View_Convert(state_u.x),state_u.p,View_Convert(state_u.x_mesh),state_u.p_mesh,state_f.x,state_f.p,state_f.x_mesh,state_f.p_mesh);
    }
#else
    if(Supports_Constraints()){
        Add_Force_Constraints_Mesh(BASE::dynamic_point_constraints,BASE::View_Convert(state_u.x),state_u.p,View_Convert(state_u.x_mesh),state_u.p_mesh,state_f.x,state_f.p,state_f.x_mesh,state_f.p_mesh);
        Add_Force_Constraints_Mesh(BASE::static_point_constraints,BASE::View_Convert(state_u.x),state_u.p,View_Convert(state_u.x_mesh),state_u.p_mesh,state_f.x,state_f.p,state_f.x_mesh,state_f.p_mesh);
        //Add_Force_Constraints_Mesh(BASE::collision_constraints,BASE::View_Convert(state_u.x),state_u.p,View_Convert(state_u.x_mesh),state_u.p_mesh,state_f.x,state_f.p,state_f.x_mesh,state_f.p_mesh);
    }
#endif
}
//#####################################################################
// Function Add_Force_Constraints_Mesh
//#####################################################################
template<class T_NONLINEAR_ELASTICITY> void HYBRID_NONLINEAR_ELASTICITY<T_NONLINEAR_ELASTICITY>::
Add_Force_Constraints_Mesh(const ARRAY<CONSTRAINT_SEGMENT<T,d> >& constraints, T_VECTOR_VARIABLE_VIEW_CONST u,T_SCALAR_VARIABLE_VIEW_CONST p,T_VECTOR_VARIABLE_MESH_VIEW_CONST u_mesh,T_SCALAR_VARIABLE_MESH_VIEW_CONST p_mesh, T_VECTOR_VARIABLE_VIEW f,T_SCALAR_VARIABLE_VIEW q,T_VECTOR_VARIABLE_MESH_VIEW f_mesh,T_SCALAR_VARIABLE_MESH_VIEW q_mesh) const
{
    LOG::SCOPE scope("HYBRID_NONLINEAR_ELASTICITY::Add_Force_Constraints_Mesh()");

    VECTOR<TV,2> X, G;  // G is the vector of endpoint forces
//    VECTOR<TV,2> computed_velocity;
        
        for(int c=1;c<=constraints.m;c++){
            
            // PHASE 1
            for(int e=1; e<=2; e++){
                    const CONSTRAINT_NODE<T,d>& endpoint = constraints(c).endpoints[e-1];
                    switch(endpoint.type){
                    case CONSTRAINT_NODE<T,d>::GRID_FIXED:
                        X(e) = Deformation_Grid(endpoint.grid_index(),endpoint.multilinear_coordinates(), u);
                        break;
                    case CONSTRAINT_NODE<T,d>::MESH_FIXED:
                        X(e) = Deformation_Mesh(endpoint.mesh_index(),endpoint.multilinear_coordinates(), u, u_mesh);
                        break;
                    case CONSTRAINT_NODE<T,d>::KINEMATIC:
                        X(e) = endpoint.spatial_location();
                        break;
                    case CONSTRAINT_NODE<T,d>::KINEMATIC_PTR:
                        X(e) = endpoint.spatial_location_ptr(BASE::collision_spring_locations);
                        break;
                    }
            }
            // PHASE 2
            TV force = X(1) - X(2);
            if( constraints(c).is_reference )
                force *= this->collision_spring_constants(constraints(c).spring_coefficient_ptr);
            else
                force *= constraints(c).spring_coefficient;
            G(1) = -force;  //Each endpoint undergoes equal and opposite force
            G(2) = force;
            
            // PHASE 3
            for(int e=1; e<=2; e++){
                    const CONSTRAINT_NODE<T,d>& endpoint = constraints(c).endpoints[e-1];
                    switch(endpoint.type){                       
                    case CONSTRAINT_NODE<T,d>::GRID_FIXED:
                        {
                            const T_INDEX& grid_index = endpoint.grid_index();
                            const T_STENCIL interpolation_stencil=Multilinear_Interpolation_Stencil_Grid(T_INDEX(),endpoint.multilinear_coordinates());
                            for(T_CONST_STENCIL_ITERATOR stencil_iterator(interpolation_stencil);stencil_iterator.Valid();stencil_iterator.Next())
                                for(int v=1;v<=d;v++)
                                    f(v)(stencil_iterator.Key()+grid_index)+=G(e)(v)*stencil_iterator.Data();}
                        break;
                    
                    case CONSTRAINT_NODE<T,d>::MESH_FIXED:
                        {
                            const T_INDEX& grid_index = cell_indices_mesh( endpoint.mesh_index() );
                            const T_STENCIL interpolation_stencil=Multilinear_Interpolation_Stencil_Mesh(endpoint.multilinear_coordinates());
                            int flat_index=1;
                            for(T_CONST_STENCIL_ITERATOR stencil_iterator(interpolation_stencil);stencil_iterator.Valid();stencil_iterator.Next(),flat_index++){
                                int mesh_node = cells_mesh(endpoint.mesh_index())(flat_index);
                                if( mesh_node == 0){
                                    for(int v=1;v<=d;v++)  f(v)(stencil_iterator.Key()+grid_index)+=G(e)(v)*stencil_iterator.Data();
                                }
                                else{
                                    for(int v=1;v<=d;v++)  f_mesh(v)(mesh_node)+=G(e)(v)*stencil_iterator.Data();                            
                                }
                            }
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
// Function Add_Force_Muscles_Mesh
//#####################################################################
template<class T_NONLINEAR_ELASTICITY> void HYBRID_NONLINEAR_ELASTICITY<T_NONLINEAR_ELASTICITY>::
Add_Force_Muscles_Mesh(T_VECTOR_VARIABLE_VIEW_CONST u,T_SCALAR_VARIABLE_VIEW_CONST p,T_VECTOR_VARIABLE_MESH_VIEW_CONST u_mesh,T_SCALAR_VARIABLE_MESH_VIEW_CONST p_mesh, T_VECTOR_VARIABLE_VIEW f,T_SCALAR_VARIABLE_VIEW q,T_VECTOR_VARIABLE_MESH_VIEW f_mesh,T_SCALAR_VARIABLE_MESH_VIEW q_mesh) const
{
    PHYSBAM_NOT_IMPLEMENTED();
}
//#####################################################################
// Function Add_Force_Differential_Muscles_And_Constraints
//#####################################################################
template<class T_NONLINEAR_ELASTICITY> void HYBRID_NONLINEAR_ELASTICITY<T_NONLINEAR_ELASTICITY>::
Add_Force_Differential_Muscles_And_Constraints(const T_STATE& state_du, T_STATE& state_df) const
{
    //LOG::SCOPE scope("HYBRID_NONLINEAR_ELASTICITY::Add_Force_Differential_Muscles_and_Constraints()");
#ifndef USE_SPECIALIZED_KERNELS
    if(Supports_Muscles())
        Add_Force_Differential_Muscles_Mesh(BASE::View_Convert(state_du.x),state_du.p,View_Convert(state_du.x_mesh),state_du.p_mesh,
                                            state_df.x,state_df.p,state_df.x_mesh,state_df.p_mesh);
    if(Supports_Constraints()){
        Add_Force_Differential_Constraints_Mesh(BASE::dynamic_point_constraints,
                                                BASE::View_Convert(state_du.x),state_du.p,View_Convert(state_du.x_mesh),state_du.p_mesh,
                                                state_df.x,state_df.p,state_df.x_mesh,state_df.p_mesh);
        Add_Force_Differential_Constraints_Mesh(BASE::static_point_constraints,
                                                BASE::View_Convert(state_du.x),state_du.p,View_Convert(state_du.x_mesh),state_du.p_mesh,
                                                state_df.x,state_df.p,state_df.x_mesh,state_df.p_mesh);
        Add_Force_Differential_Constraints_Mesh(BASE::collision_constraints,
                                                BASE::View_Convert(state_du.x),state_du.p,View_Convert(state_du.x_mesh),state_du.p_mesh,
                                                state_df.x,state_df.p,state_df.x_mesh,state_df.p_mesh);
    }
#else
    if(Supports_Constraints()){
        Add_Force_Differential_Constraints_Mesh(BASE::dynamic_point_constraints,
                                                BASE::View_Convert(state_du.x),state_du.p,View_Convert(state_du.x_mesh),state_du.p_mesh,
                                                state_df.x,state_df.p,state_df.x_mesh,state_df.p_mesh);
        Add_Force_Differential_Constraints_Mesh(BASE::static_point_constraints,
                                                BASE::View_Convert(state_du.x),state_du.p,View_Convert(state_du.x_mesh),state_du.p_mesh,
                                                state_df.x,state_df.p,state_df.x_mesh,state_df.p_mesh);
        //Add_Force_Differential_Constraints_Mesh(BASE::collision_constraints,
        //                                        BASE::View_Convert(state_du.x),state_du.p,View_Convert(state_du.x_mesh),state_du.p_mesh,
        //                                        state_df.x,state_df.p,state_df.x_mesh,state_df.p_mesh);
    }
#endif


}
//#####################################################################
// Function Add_Force_Differential_Constraints
//#####################################################################
template<class T_NONLINEAR_ELASTICITY> void HYBRID_NONLINEAR_ELASTICITY<T_NONLINEAR_ELASTICITY>::
Add_Force_Differential_Constraints_Mesh(const ARRAY<CONSTRAINT_SEGMENT<T,d> >& constraints,T_VECTOR_VARIABLE_VIEW_CONST du,T_SCALAR_VARIABLE_VIEW_CONST dp,T_VECTOR_VARIABLE_MESH_VIEW_CONST du_mesh, T_SCALAR_VARIABLE_MESH_VIEW_CONST dp_mesh,T_VECTOR_VARIABLE_VIEW df,T_SCALAR_VARIABLE_VIEW dq,T_VECTOR_VARIABLE_MESH_VIEW df_mesh,T_SCALAR_VARIABLE_MESH_VIEW dq_mesh) const
{
    //LOG::SCOPE scope("HYBRID_NONLINEAR_ELASTICITY::Add_Force_Differential_Constraints()");

    VECTOR<TV,2> dX,dG; // dG is the vector of endpoint forces
    
    for(int c=1;c<=constraints.m;c++){
        
        // PHASE 1 - Interpolate or collect displacements
        for(int e=1; e<=2; e++){
            const CONSTRAINT_NODE<T,d>& endpoint = constraints(c).endpoints[e-1];
            switch(endpoint.type){
            case CONSTRAINT_NODE<T,d>::GRID_FIXED:
                dX(e) = Displacement_Grid(endpoint.grid_index(),endpoint.multilinear_coordinates(),du);
                break;
            case CONSTRAINT_NODE<T,d>::MESH_FIXED:
                dX(e) = Displacement_Mesh(endpoint.mesh_index(),endpoint.multilinear_coordinates(),du,du_mesh);
                break;
            case CONSTRAINT_NODE<T,d>::KINEMATIC:
            case CONSTRAINT_NODE<T,d>::KINEMATIC_PTR:
                dX(e) = TV(); // This is differential, fixed_points evaluate to zero
                break;
                
            }}
        
        // PHASE 2
        TV dforce=dX(1)-dX(2);
        if(constraints(c).is_reference)
            dforce *= this->collision_spring_constants(constraints(c).spring_coefficient_ptr);
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
                    const T_INDEX& grid_index = endpoint.grid_index();
                    const T_STENCIL interpolation_stencil=Multilinear_Interpolation_Stencil_Grid(T_INDEX(),endpoint.multilinear_coordinates());
                    for(T_CONST_STENCIL_ITERATOR stencil_iterator(interpolation_stencil);stencil_iterator.Valid();stencil_iterator.Next())
                        for(int v=1;v<=d;v++)
                            df(v)(stencil_iterator.Key()+grid_index)+=dG(e)(v)*stencil_iterator.Data();
                }
                break;
            case CONSTRAINT_NODE<T,d>::KINEMATIC:
            case CONSTRAINT_NODE<T,d>::KINEMATIC_PTR:
                // Do nothing. Free space locations do not receive forces.
                break;                    
            case CONSTRAINT_NODE<T,d>::MESH_FIXED:
                {
                    const T_INDEX& grid_index = cell_indices_mesh( endpoint.mesh_index() );
                    const T_STENCIL interpolation_stencil=Multilinear_Interpolation_Stencil_Mesh(endpoint.multilinear_coordinates());
                    int flat_index=1;
                    for(T_CONST_STENCIL_ITERATOR stencil_iterator(interpolation_stencil);stencil_iterator.Valid();stencil_iterator.Next(),flat_index++){
                        int mesh_node = cells_mesh(endpoint.mesh_index())(flat_index);
                        if( mesh_node == 0){
                            for(int v=1;v<=d;v++)  df(v)(stencil_iterator.Key()+grid_index)+=dG(e)(v)*stencil_iterator.Data();
                        }
                        else{
                            for(int v=1;v<=d;v++)  df_mesh(v)(mesh_node)+=dG(e)(v)*stencil_iterator.Data();                            
                        }
                    }
                }
                break;
            }
        }
    }
}
//#####################################################################
// Function Add_Force_Differential_Muscles
//#####################################################################
template<class T_NONLINEAR_ELASTICITY> void HYBRID_NONLINEAR_ELASTICITY<T_NONLINEAR_ELASTICITY>::
Add_Force_Differential_Muscles_Mesh(T_VECTOR_VARIABLE_VIEW_CONST du,T_SCALAR_VARIABLE_VIEW_CONST dp,T_VECTOR_VARIABLE_MESH_VIEW_CONST du_mesh, T_SCALAR_VARIABLE_MESH_VIEW_CONST dp_mesh,T_VECTOR_VARIABLE_VIEW df,T_SCALAR_VARIABLE_VIEW dq,T_VECTOR_VARIABLE_MESH_VIEW df_mesh,T_SCALAR_VARIABLE_MESH_VIEW dq_mesh) const
{
    PHYSBAM_NOT_IMPLEMENTED();
}

//#####################################################################
// Function Cell_Has_Constraint
//#####################################################################
template<class T_NONLINEAR_ELASTICITY> bool HYBRID_NONLINEAR_ELASTICITY<T_NONLINEAR_ELASTICITY>::
Cell_Has_Constraint( const T_INDEX& index, int mesh_index ){
    if( mesh_index == 0 ){
        for( int constraint = 1; constraint <= BASE::collision_constraints.m; constraint++ ){
            CONSTRAINT_SEGMENT<T,d>& pc = BASE::collision_constraints(constraint);
            CONSTRAINT_NODE<T,d>& e1 = pc.endpoints[1];
            if( e1.type == CONSTRAINT_NODE<T,d>::GRID_FIXED ){
                const T_INDEX& cindex = e1.grid_index();
                if( index == cindex)
                    return true;
            }
        }  
    }
    else{
        for( int constraint = 1; constraint <= BASE::collision_constraints.m; constraint++ ){
            CONSTRAINT_SEGMENT<T,d>& pc = BASE::collision_constraints(constraint);
            CONSTRAINT_NODE<T,d>& e1 = pc.endpoints[1];
            if( e1.type == CONSTRAINT_NODE<T,d>::MESH_FIXED ){
                int mindex = e1.mesh_index();
                if( mesh_index == mindex)
                    return true;
            }
        }  
    }
    return false;
}



template class HYBRID_NONLINEAR_ELASTICITY<SKINNING_NONLINEAR_ELASTICITY<float,2,true,true> >;
template class HYBRID_NONLINEAR_ELASTICITY<SKINNING_NONLINEAR_ELASTICITY<float,2,true,false> >;
template class HYBRID_NONLINEAR_ELASTICITY<SKINNING_NONLINEAR_ELASTICITY<float,3,true,true> >;
template class HYBRID_NONLINEAR_ELASTICITY<SKINNING_NONLINEAR_ELASTICITY<float,3,true,false> >;
#ifndef USE_SPECIALIZED_KERNELS
template class HYBRID_NONLINEAR_ELASTICITY<SKINNING_NONLINEAR_ELASTICITY<double,2,true,true> >;
template class HYBRID_NONLINEAR_ELASTICITY<SKINNING_NONLINEAR_ELASTICITY<double,2,true,false> >;
template class HYBRID_NONLINEAR_ELASTICITY<SKINNING_NONLINEAR_ELASTICITY<double,3,true,true> >;
template class HYBRID_NONLINEAR_ELASTICITY<SKINNING_NONLINEAR_ELASTICITY<double,3,true,false> >;
#endif

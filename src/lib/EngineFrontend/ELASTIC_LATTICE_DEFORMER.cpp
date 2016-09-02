//#####################################################################
// Copyright 2013, Eftychios Sifakis, Nathan Mitchell
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#include <PhysBAM_Tools/Utilities/PROCESS_UTILITIES.h>
#include <PhysBAM_Tools/Random_Numbers/RANDOM_NUMBERS.h>

#include <EngineInterface/CONSTRAINTS.h>

#ifdef ENABLE_MPI_OFFLOAD
#include <EngineInterfaceMPI/ENGINE_INTERFACE_MPI.h>
#else
#include <EngineInterfaceLocal/ENGINE_INTERFACE_LOCAL.h>
#endif

#include "ELASTIC_LATTICE_DEFORMER.h"
#include <Common_Tools/Math_Tools/RANGE_ITERATOR.h>


using namespace PhysBAM;

extern bool boundary_moved;

ELASTIC_LATTICE_DEFORMER::ELASTIC_LATTICE_DEFORMER() : engine_created(false), solver_initialized(false), domain_initialized(false), collision_shape(NULL), discretization(NULL){
    pthread_mutex_init(&simulate_lock, 0);
    pthread_mutex_init(&geometry_buffer_lock, 0);
    working_vertices = &vertices;
    display_vertices = &vertices_B;
    working_stress = &stress;
    display_stress = &stress_B;
    working_strain = &strain;
    display_strain = &strain_B;
    simulation_queue = new PTHREAD_QUEUE(1);
    collision_proxy_set = false;
#ifdef ENABLE_MPI_OFFLOAD
    discretization = new ENGINE_INTERFACE_MPI();
#else 
    discretization = new ENGINE_INTERFACE_LOCAL();
#endif
}

ELASTIC_LATTICE_DEFORMER::~ELASTIC_LATTICE_DEFORMER(){
    pthread_mutex_destroy(&simulate_lock);
    pthread_mutex_destroy(&geometry_buffer_lock);
    delete simulation_queue;
    if( collision_shape)
        delete collision_shape;
    if( discretization )
        delete discretization;
}

//#####################################################################
// Add_Embedded_Point_To_Fixed_Point_Spring_Constraint
//#####################################################################
int ELASTIC_LATTICE_DEFORMER::
Add_Embedded_Point_To_Fixed_Point_Spring_Constraint(const T spring_coefficient,const TV& embedded_point_material_space_location,const TV& fixed_point_world_space_location)
{
    T_INDEX cell_index=fine_grid.Cell(embedded_point_material_space_location,0);
    PHYSBAM_ASSERT(unpadded_fine_domain.Lazy_Inside(cell_index));

    LOG::cout << embedded_point_material_space_location << " " << cell_index << std::endl;

    CONSTRAINT_SEGMENT<T,d> new_constraint;
    new_constraint.is_reference = false;

    if( ! voxmap(cell_index) ){
        LOG::cout << "Adding a dummy constraint!" << std::endl;
        new_constraint.endpoints[0].type=CONSTRAINT_NODE<T,d>::KINEMATIC;
        new_constraint.endpoints[1].type=CONSTRAINT_NODE<T,d>::KINEMATIC;
        new_constraint.endpoints[0].spatial_location()=TV();
        new_constraint.endpoints[1].spatial_location()=TV();
        return Add_Discretization_Constraint(new_constraint);
    }

    TV multilinear_coordinates=(embedded_point_material_space_location-fine_grid.Node(cell_index))/h;
    PHYSBAM_ASSERT(multilinear_coordinates.Min()>-1e-4 && multilinear_coordinates.Max()<1+1e-4);

    //LOG::cout << "Adding fine single point constraint." << std::endl;
    //LOG::cout << "Mesh index is " << fine_to_coarsemesh( cell_index ) << std::endl;

    new_constraint.endpoints[0].type=CONSTRAINT_NODE<T,d>::KINEMATIC;
    new_constraint.endpoints[0].spatial_location()=TV();
    new_constraint.endpoints[1].type=CONSTRAINT_NODE<T,d>::GRID_FIXED;
    new_constraint.endpoints[1].grid_index()=cell_index;
    new_constraint.endpoints[1].multilinear_coordinates() = multilinear_coordinates;
    new_constraint.spring_coefficient = spring_coefficient;

    return Add_Discretization_Constraint(new_constraint);
}



//#####################################################################
//  Add_Two_Embedded_Point_Spring_Constraint
//#####################################################################
int ELASTIC_LATTICE_DEFORMER::
Add_Two_Embedded_Point_Spring_Constraint(const T spring_coefficient,const TV& embedded_point_material_space_location1,const TV& embedded_point_material_space_location2)
{

    T_INDEX cell_index1=fine_grid.Cell(embedded_point_material_space_location1,0);
    PHYSBAM_ASSERT(unpadded_fine_domain.Lazy_Inside(cell_index1));
    PHYSBAM_ASSERT( voxmap(cell_index1) );
    TV multilinear_coordinates1=(embedded_point_material_space_location1-fine_grid.Node(cell_index1))/h;
    PHYSBAM_ASSERT(multilinear_coordinates1.Min()>-1e-4 && multilinear_coordinates1.Max()<1+1e-4);

    T_INDEX cell_index2=fine_grid.Cell(embedded_point_material_space_location2,0);
    PHYSBAM_ASSERT(unpadded_fine_domain.Lazy_Inside(cell_index2));
    PHYSBAM_ASSERT( voxmap(cell_index2) );
    TV multilinear_coordinates2=(embedded_point_material_space_location2-fine_grid.Node(cell_index2))/h;
    PHYSBAM_ASSERT(multilinear_coordinates2.Min()>-1e-4 && multilinear_coordinates2.Max()<1+1e-4);
    
    CONSTRAINT_SEGMENT<T,d> new_constraint;
    //LOG::cout << "Grid cell 1 " << cell_index1 << std::endl;
    //LOG::cout << "Grid cell 2 " << cell_index2 << std::endl;

    new_constraint.is_reference = false;
    new_constraint.endpoints[0].type=CONSTRAINT_NODE<T,d>::GRID_FIXED;
    new_constraint.endpoints[0].grid_index()=cell_index1;
    new_constraint.endpoints[0].multilinear_coordinates() = multilinear_coordinates1;
    new_constraint.endpoints[1].type=CONSTRAINT_NODE<T,d>::GRID_FIXED;
    new_constraint.endpoints[1].grid_index()=cell_index2;
    new_constraint.endpoints[1].multilinear_coordinates() = multilinear_coordinates2;
    new_constraint.spring_coefficient = spring_coefficient;
    
    return Add_Discretization_Constraint(new_constraint);

}

//#####################################################################
//  Update_Discretization_Constraints
//#####################################################################
int ELASTIC_LATTICE_DEFORMER::
Add_Discretization_Constraint(const CONSTRAINT_SEGMENT<T,d>& fine_constraint)
{
    if( fine_constraint.isSinglePoint() ){
        const T& spring_coefficient = fine_constraint.spring_coefficient;
        PHYSBAM_ASSERT( (fine_constraint.endpoints[1].type == CONSTRAINT_NODE<T,d>::GRID_FIXED) );
        const T_INDEX& cell_index = fine_constraint.endpoints[1].grid_index();
        const TV& multilinear_coordinates = fine_constraint.endpoints[1].multilinear_coordinates();
        PHYSBAM_ASSERT( (fine_constraint.endpoints[0].type == CONSTRAINT_NODE<T,d>::KINEMATIC) );
        const TV& world_space_location = fine_constraint.endpoints[0].spatial_location();

        TV embedded_location = (multilinear_coordinates * h) + fine_grid.Node(cell_index);
        int coarse_cid;

        {
            T_INDEX coarsecell_index=discretization->Cell(embedded_location);
            PHYSBAM_ASSERT(discretization->Unpadded_Cell_Domain().Lazy_Inside(coarsecell_index));
            TV multilinear_coordinates=(embedded_location-coarse_grid.Node(coarsecell_index))/coarse_h;
            PHYSBAM_ASSERT(multilinear_coordinates.Min()>-1e-4 && multilinear_coordinates.Max()<1+1e-4);
            
            coarse_cid=discretization->CreateNewConstraint(ENGINE_INTERFACE::DYNAMIC);
            CONSTRAINT_SEGMENT<T,d> new_constraint;
            discretization->GetConstraint(ENGINE_INTERFACE::DYNAMIC, coarse_cid, new_constraint);

            new_constraint.is_reference = false;
            
            if(fine_to_coarsemesh(cell_index) == 0 ){
                new_constraint.endpoints[0].type=CONSTRAINT_NODE<T,d>::KINEMATIC;
                new_constraint.endpoints[0].spatial_location()=TV();
                new_constraint.endpoints[1].type=CONSTRAINT_NODE<T,d>::GRID_FIXED;
                new_constraint.endpoints[1].grid_index()=coarsecell_index;
                new_constraint.endpoints[1].multilinear_coordinates() = multilinear_coordinates;
                new_constraint.spring_coefficient = spring_coefficient;           
            }
            else{
                //LOG::cout << "Adding a mesh single point constraint." << std::endl;
                int mesh_cell = fine_to_coarsemesh(cell_index);
                PHYSBAM_ASSERT( mesh_cell > 0 );
                new_constraint.endpoints[0].type=CONSTRAINT_NODE<T,d>::KINEMATIC;
                new_constraint.endpoints[0].spatial_location()=TV();
                new_constraint.endpoints[1].type=CONSTRAINT_NODE<T,d>::MESH_FIXED;
                new_constraint.endpoints[1].mesh_index()=mesh_cell;
                new_constraint.endpoints[1].multilinear_coordinates() = multilinear_coordinates;
                new_constraint.spring_coefficient = spring_coefficient; 
            }
            discretization->SetConstraint(ENGINE_INTERFACE::DYNAMIC, coarse_cid, new_constraint);
        }

        return coarse_cid;
    }

    if( fine_constraint.isDualPoint() ){
        const T& spring_coefficient = fine_constraint.spring_coefficient;
        PHYSBAM_ASSERT((fine_constraint.endpoints[0].type == CONSTRAINT_NODE<T,d>::GRID_FIXED));
        const T_INDEX& cell_index1 = fine_constraint.endpoints[0].grid_index();
        const TV& multilinear_coordinates1 = fine_constraint.endpoints[0].multilinear_coordinates();

        PHYSBAM_ASSERT((fine_constraint.endpoints[1].type == CONSTRAINT_NODE<T,d>::GRID_FIXED));
        const T_INDEX& cell_index2 = fine_constraint.endpoints[1].grid_index();
        const TV& multilinear_coordinates2 = fine_constraint.endpoints[1].multilinear_coordinates();

        TV embedded_location1 = (multilinear_coordinates1 * h) + fine_grid.Node(cell_index1);
        TV embedded_location2 = (multilinear_coordinates2 * h) + fine_grid.Node(cell_index2);
        int coarse_cid;

        //LOG::cout << "Grid cell 1 " << cell_index1 << std::endl;
        //LOG::cout << "Grid cell 2 " << cell_index2 << std::endl;

        int is_mesh_1 = fine_to_coarsemesh(cell_index1);
        int is_mesh_2 = fine_to_coarsemesh(cell_index2);
        
        //LOG::cout << "Mesh cell 1 " << is_mesh_1 << std::endl;
        //LOG::cout << "Mesh cell 2 " << is_mesh_2 << std::endl;


        {
            T_INDEX cell_index1=discretization->Cell(embedded_location1);
            T_INDEX cell_index2=discretization->Cell(embedded_location2);
            PHYSBAM_ASSERT(discretization->Unpadded_Cell_Domain().Lazy_Inside(cell_index1));
            PHYSBAM_ASSERT(discretization->Unpadded_Cell_Domain().Lazy_Inside(cell_index2));
            TV multilinear_coordinates_1=(embedded_location1-coarse_grid.Node(cell_index1))/coarse_h;
            TV multilinear_coordinates_2=(embedded_location2-coarse_grid.Node(cell_index2))/coarse_h;
            PHYSBAM_ASSERT(multilinear_coordinates_1.Min()>-1e-4 && multilinear_coordinates_1.Max()<1+1e-4);
            PHYSBAM_ASSERT(multilinear_coordinates_2.Min()>-1e-4 && multilinear_coordinates_2.Max()<1+1e-4);
            
            coarse_cid=discretization->CreateNewConstraint(ENGINE_INTERFACE::DYNAMIC);
            CONSTRAINT_SEGMENT<T,d> new_constraint;
            discretization->GetConstraint(ENGINE_INTERFACE::DYNAMIC, coarse_cid, new_constraint);
            new_constraint.is_reference = false;

            if( is_mesh_1 == 0 && is_mesh_2 == 0){
                new_constraint.endpoints[0].type=CONSTRAINT_NODE<T,d>::GRID_FIXED;
                new_constraint.endpoints[0].grid_index()=cell_index1;
                new_constraint.endpoints[0].multilinear_coordinates() = multilinear_coordinates_1;
                new_constraint.endpoints[1].type=CONSTRAINT_NODE<T,d>::GRID_FIXED;
                new_constraint.endpoints[1].grid_index()=cell_index2;
                new_constraint.endpoints[1].multilinear_coordinates() = multilinear_coordinates_2;
                new_constraint.spring_coefficient = spring_coefficient; 
            }
            else if(is_mesh_1 == 0 && is_mesh_2 > 0){
                int mesh_cell2 = is_mesh_2;
                PHYSBAM_ASSERT( mesh_cell2 > 0 );
                new_constraint.endpoints[0].type=CONSTRAINT_NODE<T,d>::GRID_FIXED;
                new_constraint.endpoints[0].grid_index()=cell_index1;
                new_constraint.endpoints[0].multilinear_coordinates() = multilinear_coordinates_1;
                new_constraint.endpoints[1].type=CONSTRAINT_NODE<T,d>::MESH_FIXED;
                new_constraint.endpoints[1].mesh_index()=mesh_cell2;
                new_constraint.endpoints[1].multilinear_coordinates() = multilinear_coordinates_2;
                new_constraint.spring_coefficient = spring_coefficient; 
            }
            else if(is_mesh_1 > 0 && is_mesh_2 == 0){
                int mesh_cell1 = is_mesh_1;
                PHYSBAM_ASSERT( mesh_cell1 > 0 );
                new_constraint.endpoints[0].type=CONSTRAINT_NODE<T,d>::MESH_FIXED;
                new_constraint.endpoints[0].mesh_index()=mesh_cell1;
                new_constraint.endpoints[0].multilinear_coordinates() = multilinear_coordinates_1;
                new_constraint.endpoints[1].type=CONSTRAINT_NODE<T,d>::GRID_FIXED;
                new_constraint.endpoints[1].grid_index()=cell_index2;
                new_constraint.endpoints[1].multilinear_coordinates() = multilinear_coordinates_2;
                new_constraint.spring_coefficient = spring_coefficient;
            }
            else{
                int mesh_cell1 = is_mesh_1;
                PHYSBAM_ASSERT( mesh_cell1 > 0 );
                int mesh_cell2 = is_mesh_2;
                PHYSBAM_ASSERT( mesh_cell2 > 0 );
                new_constraint.endpoints[0].type=CONSTRAINT_NODE<T,d>::MESH_FIXED;
                new_constraint.endpoints[0].mesh_index()=mesh_cell1;
                new_constraint.endpoints[0].multilinear_coordinates() = multilinear_coordinates_1;
                new_constraint.endpoints[1].type=CONSTRAINT_NODE<T,d>::MESH_FIXED;
                new_constraint.endpoints[1].mesh_index()=mesh_cell2;
                new_constraint.endpoints[1].multilinear_coordinates() = multilinear_coordinates_2;
                new_constraint.spring_coefficient = spring_coefficient;
            }
            discretization->SetConstraint(ENGINE_INTERFACE::DYNAMIC, coarse_cid, new_constraint);
        }
        return coarse_cid;
    }

    // If all else fails, this must be a dummy constraint.

    const TV& world_space_location = TV();
    const TV& world_space_velocity = TV();
    
    int coarse_cid=discretization->CreateNewConstraint(ENGINE_INTERFACE::DYNAMIC);
    CONSTRAINT_SEGMENT<T,d> new_constraint;
    discretization->GetConstraint(ENGINE_INTERFACE::DYNAMIC, coarse_cid, new_constraint);
    new_constraint.endpoints[0].type=CONSTRAINT_NODE<T,d>::KINEMATIC;
    new_constraint.endpoints[0].spatial_location()=world_space_location;
    new_constraint.endpoints[1].type=CONSTRAINT_NODE<T,d>::KINEMATIC;
    new_constraint.endpoints[1].spatial_location()=world_space_location;
    discretization->SetConstraint(ENGINE_INTERFACE::DYNAMIC, coarse_cid, new_constraint);

    return coarse_cid;
}


//#####################################################################
//  Update_Discretization_Constraints
//#####################################################################
void ELASTIC_LATTICE_DEFORMER::
Update_Discretization_Constraint(int cid)
{
#if 0
    CONSTRAINT_SEGMENT<T,d>& fine_constraint = fine_point_constraints(cid);
    if( fine_constraint.isSinglePoint() ){
        TV world_space_location = fine_constraint.endpoints[0].spatial_location();
        PHYSBAM_ASSERT( (fine_constraint.endpoints[0].type == CONSTRAINT_NODE<T,d>::KINEMATIC) );
        CONSTRAINT_SEGMENT<T,d> coarse_constraint;
        discretization->GetConstraint(ENGINE_INTERFACE::DYNAMIC, cid, coarse_constraint);
        PHYSBAM_ASSERT( (coarse_constraint.endpoints[0].type == CONSTRAINT_NODE<T,d>::KINEMATIC) );
        coarse_constraint.endpoints[0].spatial_location() = world_space_location;
        discretization->SetConstraint(ENGINE_INTERFACE::DYNAMIC, cid, coarse_constraint);
    }
    else{
        PHYSBAM_FATAL_ERROR( "Can not update other types of constraints!" );
    }
#endif
}


//#####################################################################
//  Update_Discretization_Constraints
//#####################################################################
int ELASTIC_LATTICE_DEFORMER::
Remove_Discretization_Constraint(int cid)
{
    int new_cid = cid;
    discretization->RemoveConstraint( ENGINE_INTERFACE::DYNAMIC, new_cid );  
    return new_cid;
}


#if 1
//#####################################################################
// Function Initialize_Collision_Points
//#####################################################################
void ELASTIC_LATTICE_DEFORMER::
Initialize_Collision_Points()
{
    LOG::SCOPE scope( "ELASTIC_LATTICE_DEFORMER::Initialize_Collision_Points()" );

    ARRAY<TV>* _vertices;
    ARRAY<T_INDEX> *_embedding_map;

    ARRAY<T> collision_constants;
    ARRAY<TV> collision_locations;

    discretization->GetCollisionConstants( collision_constants );
    discretization->GetCollisionLocations( collision_locations );

    collision_constants.Append(T());
    collision_locations.Append(TV());

    if( collision_proxy_set ){
        _vertices = &collision_proxy_vertices;
        _embedding_map=&collision_proxy_embedding_map;
    }
    else{
        _vertices = &vertices;
        _embedding_map= &embedding_map;      
    }

    LOG::cout << "Number of vertices to create proxies for: " << (*_vertices).m << std::endl;
    int constraints_created = 0;
    for( int i=1;i<=(*_vertices).m;i++){
        const T_INDEX& index = (*_embedding_map)(i);

        if( voxmap( index ) == 1){        
       
            //discretization->collision_constraints.Append(CONSTRAINT_SEGMENT<T,d>());
            int cid = discretization->CreateNewConstraint(discretization->COLLISION);
            constraints_created++;
            collision_constants.Append(T());
            collision_locations.Append(TV());
            int m = collision_constants.m; 
            CONSTRAINT_SEGMENT<T,d> pc;
            discretization->GetConstraint(discretization->COLLISION, cid, pc);

            pc.is_reference = true;
            pc.spring_coefficient_ptr = m;
            
            CONSTRAINT_NODE<T,d>& e0 = pc.endpoints[0];
            CONSTRAINT_NODE<T,d>& e1 = pc.endpoints[1];
            
            e0.type = CONSTRAINT_NODE<T,d>::KINEMATIC_PTR;
            e0._spatial_location_ptr = m;
            
            
            T_INDEX coarse_cell = ((index - 1) / refinement) + 1;
            TV weights = ( (*_vertices)(i) - coarse_grid.Node(coarse_cell))/coarse_h;       
            if(fine_to_coarsemesh( index ) == 0){
                e1.type = CONSTRAINT_NODE<T,d>::GRID_FIXED;
                e1.grid_index() = coarse_cell;
                e1.multilinear_coordinates() = weights;
            }
            else{
                e1.type = CONSTRAINT_NODE<T,d>::MESH_FIXED;
                e1.mesh_index() = fine_to_coarsemesh( index );
                e1.multilinear_coordinates() = weights;          
            }

            discretization->SetConstraint(discretization->COLLISION, cid, pc);
            
        }
    }     
    LOG::cout << "Created " << constraints_created << " collision proxies." << std::endl;

    discretization->SetCollisionConstants( collision_constants );
    discretization->SetCollisionLocations( collision_locations );
}
#endif

//#####################################################################
// Function Update_Muscles
//#####################################################################
void ELASTIC_LATTICE_DEFORMER::
Update_Muscles( const ARRAY<int, T_INDEX>& ids,  const ARRAY<int, int>& ids_mesh, 
                const ARRAY<T, T_INDEX>& density,  const ARRAY<T, int>& density_mesh,
                const ARRAY<TV, T_INDEX>& fiber,  const ARRAY<TV, int>& fiber_mesh,
                float maxstress)
{
    LOG::SCOPE scope("ELASTIC_LATTICE_DEFORMER::Update_Muscles()");


    if( muscle_id.Size().Sum() == 0 ){
        muscle_id.Resize( discretization->Unpadded_Cell_Domain() );
        muscle_density.Resize( discretization->Unpadded_Cell_Domain() );
        muscle_fiber.Resize( discretization->Unpadded_Cell_Domain() );
        muscle_id_mesh.Resize( discretization->Number_Of_Mesh_Cells() );
        muscle_density_mesh.Resize( discretization->Number_Of_Mesh_Cells() );
        muscle_fiber_mesh.Resize( discretization->Number_Of_Mesh_Cells() );
    }


    int id = max_muscle+1;

    for( RANGE_ITERATOR<d> iterator(discretization->Unpadded_Cell_Domain());iterator.Valid();iterator.Next()){
        const T_INDEX& index = iterator.Index();       
        if( ids(index) != -1 ){
            muscle_id(index).Append(id);
            muscle_density(index).Append(density(index));
            muscle_fiber(index).Append(fiber(index));
        }
    }

    int Number_Of_Mesh_Cells = discretization->Number_Of_Mesh_Cells();
    for( int cell = 1; cell < Number_Of_Mesh_Cells; cell++){
        if( ids_mesh(cell) != -1 ){
            muscle_id_mesh(cell).Append(id);
            muscle_density_mesh(cell).Append(density_mesh(cell));
            muscle_fiber_mesh(cell).Append(fiber_mesh(cell));
        }
    }

    max_muscle++;

    discretization->SetMuscleData( max_muscle,
                                   muscle_id, muscle_id_mesh,
                                   muscle_density, muscle_density_mesh,
                                   muscle_fiber, muscle_fiber_mesh);

    ARRAY<T, int> activations;
    activations.Resize( max_muscle );
    activations.Fill( 0.0 );
    discretization->SetMuscleActivations( activations );

    discretization->GetMuscleFiberMaxStress( activations );
    activations(id) = maxstress;
    discretization->SetMuscleFiberMaxStress( activations );
}

//#####################################################################
// Function Initialize_Solver
//#####################################################################
void ELASTIC_LATTICE_DEFORMER::
Initialize_Solver()
{
#if 1
    LOG::SCOPE scope("ELASTIC_LATTICE_DEFORMER::Initialize_Solver()");

    if( !discretization->EngineReady() )
        PHYSBAM_FATAL_ERROR( "Could not initialize solver. Deformer not available." );
    discretization->InitializeEngine();
    solver_initialized = true;
    frame = 0;
#endif
}

#define USE_CG

//#####################################################################
// Function Exact_Solve
//#####################################################################
float ELASTIC_LATTICE_DEFORMER::
Exact_Solve(const int krylov_iterations,const int newton_iterations,const T krylov_tolerance,const T newton_tolerance,const bool no_cut_cells)
{
#if 1
    float result;
    LOG::SCOPE scope("ELASTIC_LATTICE_DEFORMER::Exact_Solve()");
    if( !solver_initialized )
        PHYSBAM_FATAL_ERROR( "Can not perform exact solve until solver is initialized.");
    discretization->Exact_Solve(krylov_iterations,newton_iterations,krylov_tolerance,newton_tolerance,no_cut_cells, result);
    frame++;
    return result;
#endif
}
//#####################################################################
// Function Print_Volume_Change_Diagnostics
//#####################################################################
void ELASTIC_LATTICE_DEFORMER::
Print_Volume_Change_Diagnostics()
{
#if 0
    T least_volume_change=1.;
    T greatest_volume_change=0.;
    T accumulated_volume_change=0.;
    int number_of_cells=0;
    int p_m_01=0,p_m_05=0,p_m_1=0,p_m_2=0,p_m_more=0;
    const int vertices_per_cell=(d==2)?4:8;

    for(RANGE_ITERATOR<d> cell_iterator(discretization->unpadded_cell_domain);cell_iterator.Valid();cell_iterator.Next()){
        const T_INDEX& cell_index=cell_iterator.Index();

        if(discretization->cell_type(cell_index)==INTERIOR_CELL_TYPE || discretization->cell_type(cell_index)==BOUNDARY_CELL_TYPE){
            MATRIX_MXN<T> Ds;Ds.Resize(d,vertices_per_cell);

            {int vertex=1;
            for(RANGE_ITERATOR<d> iterator(RANGE<T_INDEX>(cell_index,cell_index+1));iterator.Valid();iterator.Next(),vertex++){
                const T_INDEX& index=iterator.Index();
                for(int v=1;v<=d;v++)
                    Ds(v,vertex)=u_internal.x(v)(index)+discretization->grid.Node(index)(v);}}

	    MATRIX<T,d> F=MATRIX<T,d>(Ds*discretization->G_One.Transposed());
	    T J=F.Determinant();
	    accumulated_volume_change+=J;
	    number_of_cells++;

	    if(J<least_volume_change) least_volume_change=J;
	    if(J>greatest_volume_change) greatest_volume_change=J;

	    if(J>=0.99 && J<=1.01) p_m_01++;
	    else if(J>=0.95 && J<=1.05) p_m_05++;
	    else if(J>=0.9 && J<=1.1) p_m_1++;
	    else if(J>=0.8 && J<=1.2) p_m_2++;
	    else if(J<0.8 || J>1.2) p_m_more++;
        }
    }

    LOG::cout<<"Least J = "<<least_volume_change<<std::endl;
    LOG::cout<<"Greatest J = "<<greatest_volume_change<<std::endl;
    LOG::cout<<"Average J   =  "<<accumulated_volume_change/(T)number_of_cells<<std::endl;
    LOG::cout<<"+/-.01        +/-.05        +/-.1        +/-.2        +/->.2"<<std::endl;
    LOG::cout<<p_m_01<<"           "<<p_m_05<<"             "<<p_m_1<<"           "<<p_m_2<<"            "<<p_m_more<<std::endl;
#endif
}
//#####################################################################
// Function Print_Force_Diagnostics
//#####################################################################
void ELASTIC_LATTICE_DEFORMER::
Print_Force_Diagnostics(const bool no_cut_cells)
{
#if 0
    HYBRID_CG_SYSTEM<T,d,T_DISCRETIZATION> krylov_system(*discretization);
    HYBRID_CG_VECTOR<T,d,T_DISCRETIZATION> krylov_b_debug(*discretization,b_X_debug,b_P_debug);

    LOG::SCOPE scope("Computing off-balance forces");
    for(int v=1;v<=d;v++) b_X_debug(v).Fill(T());
    b_P_debug.Fill(T());

    discretization->Update_Position_Based_State(u_internal, false);

    if(no_cut_cells)
        discretization->Add_Force(b_X_debug,b_P_debug,true);
    else{
        LOG::cout<<"USING CUT CELLS!"<<std::endl;
        discretization->Add_Force(b_X_debug,b_P_debug,false);}
    
    b_X+=discretization->f_external;
    T norm=krylov_system.Convergence_Norm(krylov_b_debug);
    LOG::cout<<"#############################################################################"<<std::endl;
    LOG::cout<<"         Norm of out-of-balance forces :                    "<<norm<<std::endl;
    LOG::cout<<"#############################################################################"<<std::endl;
#endif
}

//#####################################################################
// Function Print_Matrix_Condition
//#####################################################################
void ELASTIC_LATTICE_DEFORMER::
Print_Matrix_Condition(const bool no_cut_cells)
{
#if 0
    CG_SYSTEM<T_DISCRETIZATION> krylov_system(*discretization);
    CG_VECTOR<CG_POLICY<T_DISCRETIZATION> > krylov_b_debug(*discretization,b_debug);
    CG_VECTOR<CG_POLICY<T_DISCRETIZATION> > krylov_x(*discretization,x);

    LOG::SCOPE scope("Probing Matrix Condition Number");
    

    T norm = 0;

    LOG::cout << "Node passes: "<< std::endl;
    {
        RANGE<T_INDEX> domain = discretization->unpadded_node_domain;
        for( int w = 1; w <= d; w++){
            LOG::cout << w << " Channel" << std::endl;

            for( RANGE_ITERATOR<d> offset_iterator(RANGE<T_INDEX>(T_INDEX(), T_INDEX::All_Ones_Vector()));offset_iterator.Valid();offset_iterator.Next()){
                RANGE<T_INDEX> offset_range( offset_iterator.Index(), offset_iterator.Index());
                RANGE<T_INDEX> shift_domain = domain + offset_range;

                krylov_b_debug.Clear();
                krylov_x.Clear();

                LOG::cout << shift_domain << std::endl;
                for( RANGE_ITERATOR<d,2> n_iter(shift_domain);n_iter.Valid();n_iter.Next()){
                    if(discretization->node_is_active(n_iter.Index())){
                        x.x(w)(n_iter.Index()) = T(1.0);
                        //x.x(1)(n_iter.Index()) = T(1.0);
                        //x.x(2)(n_iter.Index()) = T(1.0);
                        //x.x(3)(n_iter.Index()) = T(1.0);
                    }
                }
                //krylov_system.Multiply( krylov_x, krylov_b_debug);  
                discretization->Add_Force_Differential( x, b_debug);
                T t_norm = b_debug.x(w).Maxabs();
                //T t_norm = krylov_system.Convergence_Norm(krylov_b_debug);
                LOG::cout<<"         Matrix Condition Number :                    "<<t_norm<<std::endl;
                norm = max( norm, t_norm );
            }
        }
    }

    LOG::cout << "Pressure pass: "<< std::endl;
    {
        RANGE<T_INDEX> domain = discretization->unpadded_cell_domain;
        krylov_b_debug.Clear();
        krylov_x.Clear();
        for( RANGE_ITERATOR<d> c_iter(domain);c_iter.Valid();c_iter.Next()){
            if(discretization->cell_type(c_iter.Index()) == INTERIOR_CELL_TYPE ){
                x.p(c_iter.Index()) = T(1.0);
            }
        }
        //krylov_system.Multiply( krylov_x, krylov_b_debug);                             
        discretization->Add_Force_Differential( x, b_debug);
        T t_norm = b_debug.p.Maxabs();
        LOG::cout<<"         Matrix Condition Number :                    "<<t_norm<<std::endl;
        norm = max( norm, t_norm );
    }


    LOG::cout<<"#############################################################################"<<std::endl;
    LOG::cout<<"         Matrix Condition Number :                    "<<norm<<std::endl;
    LOG::cout<<"#############################################################################"<<std::endl;
    krylov_b_debug.Clear();
    krylov_x.Clear();
#endif
}

//#####################################################################
// Function Material_Convergence_Rate
//#####################################################################
void ELASTIC_LATTICE_DEFORMER::
Material_Convergence_Rate(const bool no_cut_cells)
{
#if 0
    typedef ARRAY<T,int> T_SCALAR_VARIABLE;
    typedef VECTOR<T_SCALAR_VARIABLE,d> T_VECTOR_VARIABLE;
    typedef ARRAY_VIEW<T,T_INDEX> T_SCALAR_VARIABLE_VIEW;
    typedef VECTOR<T_SCALAR_VARIABLE_VIEW,d> T_VECTOR_VARIABLE_VIEW;
    typedef ARRAY_VIEW<T,int> T_SCALAR_VARIABLE_MESH_VIEW;
    typedef VECTOR<T_SCALAR_VARIABLE_MESH_VIEW,d> T_VECTOR_VARIABLE_MESH_VIEW;

    CG_SYSTEM<T_DISCRETIZATION> krylov_system(*discretization);
    CG_VECTOR<CG_POLICY<T_DISCRETIZATION> > krylov_b_debug(*discretization,b_debug);
    CG_VECTOR<CG_POLICY<T_DISCRETIZATION> > krylov_b2_debug(*discretization,b2_debug);
    CG_VECTOR<CG_POLICY<T_DISCRETIZATION> > krylov_x(*discretization,x);
    CG_VECTOR<CG_POLICY<T_DISCRETIZATION> > krylov_r(*discretization,r);
    CG_VECTOR<CG_POLICY<T_DISCRETIZATION> > krylov_k(*discretization,k);
    CG_VECTOR<CG_POLICY<T_DISCRETIZATION> > krylov_z(*discretization,z);

    T_VECTOR_VARIABLE du;
    for(int v=1;v<=d;v++) du(v).Resize((discretization->padded_node_domain.Edge_Lengths()+1).Product(),true,false);
    for(int v=1;v<=d;v++) du(v).Fill(T(0));
    T_VECTOR_VARIABLE_VIEW du_v(T_SCALAR_VARIABLE_VIEW(discretization->padded_node_domain,du.x.Get_Array_Pointer()),
                                T_SCALAR_VARIABLE_VIEW(discretization->padded_node_domain,du.y.Get_Array_Pointer()),
                                T_SCALAR_VARIABLE_VIEW(discretization->padded_node_domain,du.z.Get_Array_Pointer()));

    LOG::SCOPE scope("Probing Convergence Rate");
    
    T norm = 0;
    for( int i = -2; i >= -5; i--){
        T  alpha = pow(10,i);
        krylov_b_debug.Clear();
        krylov_b2_debug.Clear();
        krylov_r.Clear();
        krylov_k.Clear();
        krylov_z.Clear();
        
        for(RANGE_ITERATOR<d> iterator(discretization->unpadded_node_domain);iterator.Valid();iterator.Next()){
            for(int v=1;v<=d;v++)
            z.x(v) = u_internal.x(v);
        }
        
        discretization->Update_Position_Based_State(u_internal);
        discretization->Add_Force(u_internal, b_debug);
        
        RANDOM_NUMBERS<T> random_numbers;
        random_numbers.Set_Seed(1);
        
        const TV offset = random_numbers.Get_Direction<TV>();
        for(RANGE_ITERATOR<d> iterator(discretization->unpadded_node_domain);iterator.Valid();iterator.Next()){
            const T_INDEX& node_index=iterator.Index();
            if(discretization->node_is_active(node_index)){
                for(int v=1;v<=d;v++)
                    du_v(v)(node_index)=alpha*offset(v);            
            }
        }
        
        u_internal.x += du_v;
        r.x += du_v;
        
        //discretization->Add_Force_Differential( u_internal_mesh, k );
        discretization->Add_Force_Differential( r, k );
        
        discretization->Update_Position_Based_State(u_internal);
        discretization->Add_Force(u_internal, b2_debug);
        
        krylov_b2_debug -= krylov_b_debug;
        
        LOG::cout << "################    ALPHA: "<<alpha<<"#######################################"<<std::endl;
        LOG::cout << "################    DIFFERENTIAL NUMERICAL   ################################"<<std::endl;
        norm = b2_debug.x(1).Maxabs();
        LOG::cout << "       Max in x: " << norm << std::endl;
        norm = b2_debug.x(2).Maxabs();
        LOG::cout << "       Max in y: " << norm << std::endl;
        norm = b2_debug.x(3).Maxabs();
        LOG::cout << "       Max in z: " << norm << std::endl;
        LOG::cout << "################    DIFFERENTIAL ANALYTICAL  ################################"<<std::endl;
        norm = k.x(1).Maxabs();
        LOG::cout << "       Max in x: " << norm << std::endl;
        norm = k.x(2).Maxabs();
        LOG::cout << "       Max in y: " << norm << std::endl;
        norm = k.x(3).Maxabs();
        LOG::cout << "       Max in z: " << norm << std::endl;
        LOG::cout << "################    DIFFERENCE               ################################"<<std::endl;
        norm = abs(k.x(1).Maxabs() - b2_debug.x(1).Maxabs());
        LOG::cout << "       Abs Diff in x: " << norm << std::endl;
        norm = abs(k.x(2).Maxabs() - b2_debug.x(2).Maxabs());
        LOG::cout << "       Abs Diff in y: " << norm << std::endl;
        norm = abs(k.x(3).Maxabs() - b2_debug.x(3).Maxabs());
        LOG::cout << "       Abs Diff in z: " << norm << std::endl;
        LOG::cout << "#############################################################################"<<std::endl;
        
        for(RANGE_ITERATOR<d> iterator(discretization->unpadded_node_domain);iterator.Valid();iterator.Next()){
            for(int v=1;v<=d;v++)
            u_internal.x(v) = z.x(v);
        }
    }

    krylov_b_debug.Clear();
    krylov_b2_debug.Clear();
    krylov_r.Clear();
    krylov_k.Clear();
    krylov_z.Clear();
#endif
}
// TODO: Need to handle these eventually.
/*

//#####################################################################
// Function Multilinear_Interpolation
//#####################################################################
template<class T_ELASTICITY> void ELASTIC_LATTICE_DEFORMER<T_ELASTICITY>::
Multilinear_Interpolation(const T_INDEX& cell_index,const TV& weights)
{
    TV result;
    for(RANGE_ITERATOR<d> iterator(RANGE<T_INDEX>(T_INDEX(),T_INDEX::All_Ones_Vector()));iterator.Valid();iterator.Next()){
        const T_INDEX& offset=iterator.Index();
        const T_INDEX node_index=cell_index+offset;

        TV partial_result=grid.Node(node_index);
        for(int v=1;v<=d;v++) partial_result(v)+=u(v)(node_index);

        for(int v=1;v<=d;v++)
            if(offset(v))
                partial_result*=weights(v);
            else
                partial_result*=(1.-weights(v));

        result+=partial_result;
    }

    return result;
}
//#####################################################################
// Function Interpolate_Embedded_Object
//#####################################################################
template<class T_ELASTICITY> void ELASTIC_LATTICE_DEFORMER<T_ELASTICITY>::
Interpolate_Embedded_Object()
{
    if(embedded_particles)
        for(int p=1;p<=embedded_particles->array_collection->Size();p++)
            embedded_particles->X(p)=Multilinear_Interpolation(embedding_map(p),embedding_weights(p));
    return;
}

*/



//template class ELASTIC_LATTICE_DEFORMER< HYBRID_NONLINEAR_ELASTICITY<SKINNING_NONLINEAR_ELASTICITY<float,3,true,false> > >;
//template class ELASTIC_LATTICE_DEFORMER< SKINNING_NONLINEAR_ELASTICITY<float,3,true,false> >;

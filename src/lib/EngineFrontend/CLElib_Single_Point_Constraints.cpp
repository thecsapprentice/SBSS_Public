//#####################################################################
// Function Add_Single_Point_Constraint
//#####################################################################
int CLElib::
Add_Single_Point_Constraint(const int triangle_id, const std::array<float,3>& uv){
#ifndef ENABLE_LOG_MESSAGES
    LOG::SCOPE scope("CLElib::Add_Single_Point_Constraint");
#endif
    
    typedef VECTOR<T,d> TV;
    typedef VECTOR<int,d> T_INDEX;

    PHYSBAM_ASSERT(implementation);
    PRIMARY_DEFORMER& elastic_lattice_deformer=Implementation<PRIMARY_DEFORMER>();
    if(!elastic_lattice_deformer.engine_created)
        PHYSBAM_FATAL_ERROR("Add_Single_Point_Constraint() called without first creating a model");
    PRIMARY_ELASTICITY& discretization=elastic_lattice_deformer.Discretization();


    TV location;
    
    const VECTOR<int,3>& triangle_vertices=elastic_lattice_deformer.triangles(triangle_id+1);
    const ARRAY<T_INDEX>& embedding_map=elastic_lattice_deformer.embedding_map;
    const ARRAY<TV>& embedding_weights=elastic_lattice_deformer.embedding_weights;
    TRIANGLE_3D<T> material_triangle(
        elastic_lattice_deformer.fine_grid.Node(embedding_map(triangle_vertices(1)))+embedding_weights(triangle_vertices(1))*elastic_lattice_deformer.h,
        elastic_lattice_deformer.fine_grid.Node(embedding_map(triangle_vertices(2)))+embedding_weights(triangle_vertices(2))*elastic_lattice_deformer.h,
        elastic_lattice_deformer.fine_grid.Node(embedding_map(triangle_vertices(3)))+embedding_weights(triangle_vertices(3))*elastic_lattice_deformer.h);
    location=material_triangle.Point_From_Barycentric_Coordinates(VECTOR<T,3>(uv[0], uv[1], uv[2]));
    
    int spc_id = Add_Single_Point_Constraint( std::array<float,3>({{location(1), location(2), location(3) }}));
    elastic_lattice_deformer.single_point_constraint_map.at(spc_id).triangle = triangle_id;
    Move_Single_Point_Constraint(spc_id, std::array<float,3>({{location(1), location(2), location(3) }}) );
    return spc_id;
}

//#####################################################################
// Function Add_Single_Point_Constraint
//#####################################################################
int CLElib::
Add_Single_Point_Constraint(const std::array<float,3>& location){
#ifndef ENABLE_LOG_MESSAGES
    LOG::SCOPE scope("CLElib::Add_Single_Point_Constraint");
#endif
    
    typedef VECTOR<T,d> TV;
    typedef VECTOR<int,d> T_INDEX;

    PHYSBAM_ASSERT(implementation);
    PRIMARY_DEFORMER& elastic_lattice_deformer=Implementation<PRIMARY_DEFORMER>();
    if(!elastic_lattice_deformer.engine_created)
        PHYSBAM_FATAL_ERROR("Add_Single_Point_Constraint() called without first creating a model");
    PRIMARY_ELASTICITY& discretization=elastic_lattice_deformer.Discretization();

    
    TV embedded_point_material_space_location(location[0], location[1], location[2]);

    
    pthread_mutex_lock(&elastic_lattice_deformer.simulate_lock);
    
#ifdef GRID_IN_GRID
    int cid=elastic_lattice_deformer.Add_Embedded_Point_To_Fixed_Point_Spring_Constraint(elastic_lattice_deformer.hook_stiffness,embedded_point_material_space_location,TV());
#else
    int cid=discretization.Add_Embedded_Point_To_Fixed_Point_Spring_Constraint(elastic_lattice_deformer.hook_stiffness,embedded_point_material_space_location,TV(), false);
#endif

    PRIMARY_DEFORMER::SINGLE_POINT_CONSTRAINT spc_entry;
    spc_entry.cid = cid;
    spc_entry.triangle = -1;

    pthread_mutex_unlock(&elastic_lattice_deformer.simulate_lock);

    auto status = elastic_lattice_deformer.single_point_constraint_map.insert( PRIMARY_DEFORMER::SINGLE_POINT_CONSTRAINT_MAP::value_type(elastic_lattice_deformer.single_point_constraint_maxid++, spc_entry) );

    return elastic_lattice_deformer.single_point_constraint_maxid-1;    
}

//#####################################################################
// Function Move_Single_Point_Constraint
//#####################################################################
void CLElib::
Move_Single_Point_Constraint(const int hook_id,const std::array<float,3>& location){
#ifndef ENABLE_LOG_MESSAGES
    LOG::SCOPE scope("CLElib::Move_Single_Point_Constraint");
#endif

    
    typedef VECTOR<T,d> TV;
    typedef VECTOR<int,d> T_INDEX;

    PHYSBAM_ASSERT(implementation);
    PRIMARY_DEFORMER& elastic_lattice_deformer=Implementation<PRIMARY_DEFORMER>();
    if(!elastic_lattice_deformer.engine_created)
        PHYSBAM_FATAL_ERROR("Move_Single_Point_Constraint() called without first creating a model");
    PRIMARY_ELASTICITY& discretization=elastic_lattice_deformer.Discretization();
    
    TV space_location(location[0], location[1], location[2]);

    pthread_mutex_lock(&elastic_lattice_deformer.simulate_lock);
    int cid;
    try{
        PRIMARY_DEFORMER::SINGLE_POINT_CONSTRAINT& spc_entry = elastic_lattice_deformer.single_point_constraint_map.at( hook_id );
        cid=spc_entry.cid;
    }
    catch(  const std::out_of_range& e ){
        LOG::cout << "Hook ID " << hook_id << " is invalid." << std::endl;
        pthread_mutex_unlock(&elastic_lattice_deformer.simulate_lock);
        return;
    }

    CONSTRAINT_SEGMENT<T,d> constraint;
    discretization.GetConstraint(ENGINE_INTERFACE::DYNAMIC, cid, constraint);
    PHYSBAM_ASSERT( (constraint.endpoints[0].type == CONSTRAINT_NODE<T,d>::KINEMATIC) );
    constraint.endpoints[0].spatial_location() = space_location;
    discretization.SetConstraint(ENGINE_INTERFACE::DYNAMIC, cid, constraint);
        
    pthread_mutex_unlock(&elastic_lattice_deformer.simulate_lock);
}

//#####################################################################
// Function Delete_Single_Point_Constraint
//#####################################################################
void CLElib::
Delete_Single_Point_Constraint(const int hook_id){
#ifndef ENABLE_LOG_MESSAGES
    LOG::SCOPE scope("CLElib::Delete_Single_Point_Constraint");
#endif

    
    typedef VECTOR<T,d> TV;
    typedef VECTOR<int,d> T_INDEX;

    PHYSBAM_ASSERT(implementation);
    PRIMARY_DEFORMER& elastic_lattice_deformer=Implementation<PRIMARY_DEFORMER>();
    if(!elastic_lattice_deformer.engine_created)
        PHYSBAM_FATAL_ERROR("Delete_Single_Point_Constraint() called without first creating a model");
    PRIMARY_ELASTICITY& discretization=elastic_lattice_deformer.Discretization();
    int cid, new_cid;
    
    pthread_mutex_lock(&elastic_lattice_deformer.simulate_lock);
    {
        LOG::SCOPE scope( "Removing.");
        try{
            PRIMARY_DEFORMER::SINGLE_POINT_CONSTRAINT& spc_entry = elastic_lattice_deformer.single_point_constraint_map.at( hook_id );            
            cid=spc_entry.cid;
        }
        catch( const std::out_of_range& e ){
            LOG::cout << "Hook ID " << hook_id << " is invalid." << std::endl;
            pthread_mutex_unlock(&elastic_lattice_deformer.simulate_lock);
            return;
        }
        new_cid = elastic_lattice_deformer.Remove_Discretization_Constraint(cid);
        LOG::cout << "Done." << std::endl;
        LOG::cout << "Original cid: " << cid << std::endl;
        LOG::cout << "New cid: " << new_cid << std::endl;            
    }
    pthread_mutex_unlock(&elastic_lattice_deformer.simulate_lock);

    
    elastic_lattice_deformer.single_point_constraint_map.erase( hook_id );
    for( auto & iter : elastic_lattice_deformer.single_point_constraint_map )
        if( iter.second.cid == new_cid )
            iter.second.cid = cid;

    for( auto const& iter : elastic_lattice_deformer.single_point_constraint_map )
        LOG::cout << iter.first << " : " << iter.second.cid << std::endl;
}

//#####################################################################
// Function Get_Single_Point_Constraint_Position
//#####################################################################
void CLElib::
Get_Single_Point_Constraint_Position(const std::vector<int>& ids, std::vector< std::array<float,3> >& position){
#ifndef ENABLE_LOG_MESSAGES
    LOG::SCOPE scope("CLElib::Get_Single_Point_Constraint_Position");
#endif

    
    typedef VECTOR<T,d> TV;
    typedef VECTOR<int,d> T_INDEX;

    PHYSBAM_ASSERT(implementation);
    PRIMARY_DEFORMER& elastic_lattice_deformer=Implementation<PRIMARY_DEFORMER>();
    if(!elastic_lattice_deformer.engine_created)
        PHYSBAM_FATAL_ERROR("Delete_Single_Point_Constraint() called without first creating a model");
    PRIMARY_ELASTICITY& discretization=elastic_lattice_deformer.Discretization();

    position.clear();
    
    pthread_mutex_lock(&elastic_lattice_deformer.simulate_lock);
    for(auto const& iter : ids ){
        PRIMARY_DEFORMER::SINGLE_POINT_CONSTRAINT_MAP::iterator iter2 = elastic_lattice_deformer.single_point_constraint_map.find( iter );
        CONSTRAINT_SEGMENT<T,d> constraint;
        LOG::cout << "Loading constraint with cid " << iter2->second.cid << " from hook " << iter << std::endl;
        discretization.GetConstraint(ENGINE_INTERFACE::DYNAMIC, iter2->second.cid, constraint);
        PHYSBAM_ASSERT( (constraint.endpoints[0].type == CONSTRAINT_NODE<T,d>::KINEMATIC) );
        TV sl = constraint.endpoints[0].spatial_location();
        position.push_back( std::array<float,3>( {{ sl(1), sl(2), sl(3) }} ));        
    }   
    pthread_mutex_unlock(&elastic_lattice_deformer.simulate_lock);

}

//#####################################################################
// Function Get_Single_Point_Constraint_Triangles
//#####################################################################
void CLElib::
Get_Single_Point_Constraint_Triangles(const std::vector<int>& ids, std::vector< int >& triangles){
#ifndef ENABLE_LOG_MESSAGES
    LOG::SCOPE scope("CLElib::Get_Single_Point_Constraint_Triangles");
#endif

    
    typedef VECTOR<T,d> TV;
    typedef VECTOR<int,d> T_INDEX;

    PHYSBAM_ASSERT(implementation);
    PRIMARY_DEFORMER& elastic_lattice_deformer=Implementation<PRIMARY_DEFORMER>();
    if(!elastic_lattice_deformer.engine_created)
        PHYSBAM_FATAL_ERROR("Delete_Single_Point_Constraint() called without first creating a model");
    PRIMARY_ELASTICITY& discretization=elastic_lattice_deformer.Discretization();
 
    triangles.clear();
    
    pthread_mutex_lock(&elastic_lattice_deformer.simulate_lock);
    for(auto const& iter : ids ){
        PRIMARY_DEFORMER::SINGLE_POINT_CONSTRAINT_MAP::iterator iter2 = elastic_lattice_deformer.single_point_constraint_map.find( iter );
        int tri = iter2->second.triangle;
        triangles.push_back( tri );        
    }   
    pthread_mutex_unlock(&elastic_lattice_deformer.simulate_lock);

}

//#####################################################################
// Function Get_Active_Single_Point_Constraints
//#####################################################################
void CLElib::
Get_Active_Single_Point_Constraints(std::vector<int>& ids){
#ifndef ENABLE_LOG_MESSAGES
    LOG::SCOPE scope("CLElib::Get_Active_Single_Point_Constraints");
#endif

    
    typedef VECTOR<T,d> TV;
    typedef VECTOR<int,d> T_INDEX;

    PHYSBAM_ASSERT(implementation);
    PRIMARY_DEFORMER& elastic_lattice_deformer=Implementation<PRIMARY_DEFORMER>();
    if(!elastic_lattice_deformer.engine_created)
        PHYSBAM_FATAL_ERROR("Delete_Single_Point_Constraint() called without first creating a model");
    PRIMARY_ELASTICITY& discretization=elastic_lattice_deformer.Discretization();

    
    ids.clear();
    
    for( auto const& iter : elastic_lattice_deformer.single_point_constraint_map )
        ids.push_back( iter.first );

}

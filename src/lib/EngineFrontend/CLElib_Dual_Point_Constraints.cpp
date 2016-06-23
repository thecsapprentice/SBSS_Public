//#####################################################################
// Function Add_Dual_Point_Constraint
//#####################################################################
int CLElib::
Add_Dual_Point_Constraint(const int triangleA_id, const std::array<float,3>& uv_A, 
                          const int triangleB_id, const std::array<float,3>& uv_B){
#ifndef ENABLE_LOG_MESSAGES
    LOG::SCOPE scope("CLElib::Add_Dual_Point_Constraint");
#endif
    
    typedef VECTOR<T,d> TV;
    typedef VECTOR<int,d> T_INDEX;

    PHYSBAM_ASSERT(implementation);
    PRIMARY_DEFORMER& elastic_lattice_deformer=Implementation<PRIMARY_DEFORMER>();
    if(!elastic_lattice_deformer.engine_created)
        PHYSBAM_FATAL_ERROR("Add_Dual_Point_Constraint() called without first creating a model");
    PRIMARY_ELASTICITY& discretization=elastic_lattice_deformer.Discretization();


    TV locationA;
    TV locationB;
    
    {
        const VECTOR<int,3>& triangle_vertices=elastic_lattice_deformer.triangles(triangleA_id+1);
        const ARRAY<T_INDEX>& embedding_map=elastic_lattice_deformer.embedding_map;
        const ARRAY<TV>& embedding_weights=elastic_lattice_deformer.embedding_weights;
        TRIANGLE_3D<T> material_triangle(
                                         elastic_lattice_deformer.fine_grid.Node(embedding_map(triangle_vertices(1)))+embedding_weights(triangle_vertices(1))*elastic_lattice_deformer.h,
                                         elastic_lattice_deformer.fine_grid.Node(embedding_map(triangle_vertices(2)))+embedding_weights(triangle_vertices(2))*elastic_lattice_deformer.h,
                                         elastic_lattice_deformer.fine_grid.Node(embedding_map(triangle_vertices(3)))+embedding_weights(triangle_vertices(3))*elastic_lattice_deformer.h);
        locationA=material_triangle.Point_From_Barycentric_Coordinates(VECTOR<T,3>(uv_A[0], uv_A[1], uv_A[2]));
    }
    
    {
        const VECTOR<int,3>& triangle_vertices=elastic_lattice_deformer.triangles(triangleB_id+1);
        const ARRAY<T_INDEX>& embedding_map=elastic_lattice_deformer.embedding_map;
        const ARRAY<TV>& embedding_weights=elastic_lattice_deformer.embedding_weights;
        TRIANGLE_3D<T> material_triangle(
                                         elastic_lattice_deformer.fine_grid.Node(embedding_map(triangle_vertices(1)))+embedding_weights(triangle_vertices(1))*elastic_lattice_deformer.h,
                                         elastic_lattice_deformer.fine_grid.Node(embedding_map(triangle_vertices(2)))+embedding_weights(triangle_vertices(2))*elastic_lattice_deformer.h,
                                         elastic_lattice_deformer.fine_grid.Node(embedding_map(triangle_vertices(3)))+embedding_weights(triangle_vertices(3))*elastic_lattice_deformer.h);
        locationB=material_triangle.Point_From_Barycentric_Coordinates(VECTOR<T,3>(uv_B[0], uv_B[1], uv_B[2]));
    }
    
    int spc_id = Add_Dual_Point_Constraint( std::array<float,3>({{locationA(1), locationA(2), locationA(3) }}), std::array<float,3>({{locationB(1), locationB(2), locationB(3) }}) );
    elastic_lattice_deformer.dual_point_constraint_map.at(spc_id).triangleA = triangleA_id;
    elastic_lattice_deformer.dual_point_constraint_map.at(spc_id).triangleB = triangleB_id;
    elastic_lattice_deformer.dual_point_constraint_map.at(spc_id).triA_uv = uv_A;
    elastic_lattice_deformer.dual_point_constraint_map.at(spc_id).triB_uv = uv_B;
    return spc_id;
}

//#####################################################################
// Function Add_Dual_Point_Constraint
//#####################################################################
int CLElib::
Add_Dual_Point_Constraint(const std::array<float,3>& locationA, const std::array<float,3>& locationB){
#ifndef ENABLE_LOG_MESSAGES
    LOG::SCOPE scope("CLElib::Add_Dual_Point_Constraint");
#endif
    
    typedef VECTOR<T,d> TV;
    typedef VECTOR<int,d> T_INDEX;

    PHYSBAM_ASSERT(implementation);
    PRIMARY_DEFORMER& elastic_lattice_deformer=Implementation<PRIMARY_DEFORMER>();
    if(!elastic_lattice_deformer.engine_created)
        PHYSBAM_FATAL_ERROR("Add_Dual_Point_Constraint() called without first creating a model");
    PRIMARY_ELASTICITY& discretization=elastic_lattice_deformer.Discretization();

    
    TV embedded_point_material_space_locationA(locationA[0], locationA[1], locationA[2]);
    TV embedded_point_material_space_locationB(locationB[0], locationB[1], locationB[2]);

    
    pthread_mutex_lock(&elastic_lattice_deformer.simulate_lock);
    
#ifdef GRID_IN_GRID
    int cid=elastic_lattice_deformer.Add_Two_Embedded_Point_Spring_Constraint(elastic_lattice_deformer.hook_stiffness,embedded_point_material_space_locationA,embedded_point_material_space_locationB);
#else
    int cid=discretization.Add_Two_Embedded_Point_Spring_Constraint(elastic_lattice_deformer.hook_stiffness,embedded_point_material_space_locationA,embedded_point_material_space_locationB);
#endif

    PRIMARY_DEFORMER::DUAL_POINT_CONSTRAINT spc_entry;
    spc_entry.cid = cid;
    spc_entry.triangleA = -1;
    spc_entry.triangleB = -1;
    spc_entry.triA_uv = {{-1,-1,-1,}};
    spc_entry.triB_uv = {{-1,-1,-1,}};

    pthread_mutex_unlock(&elastic_lattice_deformer.simulate_lock);

    auto status = elastic_lattice_deformer.dual_point_constraint_map.insert( PRIMARY_DEFORMER::DUAL_POINT_CONSTRAINT_MAP::value_type(elastic_lattice_deformer.dual_point_constraint_maxid++, spc_entry) );

    return elastic_lattice_deformer.dual_point_constraint_maxid-1;    
}

//#####################################################################
// Function Delete_Dual_Point_Constraint
//#####################################################################
void CLElib::
Delete_Dual_Point_Constraint(const int suture_id){
#ifndef ENABLE_LOG_MESSAGES
    LOG::SCOPE scope("CLElib::Delete_Dual_Point_Constraint");
#endif

    
    typedef VECTOR<T,d> TV;
    typedef VECTOR<int,d> T_INDEX;

    PHYSBAM_ASSERT(implementation);
    PRIMARY_DEFORMER& elastic_lattice_deformer=Implementation<PRIMARY_DEFORMER>();
    if(!elastic_lattice_deformer.engine_created)
        PHYSBAM_FATAL_ERROR("Delete_Dual_Point_Constraint() called without first creating a model");
    PRIMARY_ELASTICITY& discretization=elastic_lattice_deformer.Discretization();
    int cid, new_cid;
    
    pthread_mutex_lock(&elastic_lattice_deformer.simulate_lock);
    {        
        LOG::SCOPE scope( "Removing.");
        try{
            PRIMARY_DEFORMER::DUAL_POINT_CONSTRAINT& spc_entry = elastic_lattice_deformer.dual_point_constraint_map.at( suture_id );
            cid=spc_entry.cid;
        }
        catch( const std::out_of_range& e ){
            LOG::cout << "Suture ID " << suture_id << " is invalid." << std::endl;
            pthread_mutex_unlock(&elastic_lattice_deformer.simulate_lock);
            return;
        }

        new_cid = elastic_lattice_deformer.Remove_Discretization_Constraint(cid);
        LOG::cout << "Done." << std::endl;
        LOG::cout << "Original cid: " << cid << std::endl;
        LOG::cout << "New cid: " << new_cid << std::endl;        
    }
    pthread_mutex_unlock(&elastic_lattice_deformer.simulate_lock);

    
    elastic_lattice_deformer.dual_point_constraint_map.erase( suture_id );
    for( auto & iter : elastic_lattice_deformer.dual_point_constraint_map )
        if( iter.second.cid == new_cid )
            iter.second.cid = cid;

    for( auto const& iter : elastic_lattice_deformer.dual_point_constraint_map )
        LOG::cout << iter.first << " : " << iter.second.cid << std::endl;
}

//#####################################################################
// Function Get_Dual_Point_Constraint_Position
//#####################################################################
void CLElib::
Get_Dual_Point_Constraint_UV(const std::vector<int>& ids, std::vector< std::array<float,6> >& UVs){
#ifndef ENABLE_LOG_MESSAGES
    LOG::SCOPE scope("CLElib::Get_Dual_Point_Constraint_Position");
#endif

    
    typedef VECTOR<T,d> TV;
    typedef VECTOR<int,d> T_INDEX;

    PHYSBAM_ASSERT(implementation);
    PRIMARY_DEFORMER& elastic_lattice_deformer=Implementation<PRIMARY_DEFORMER>();
    if(!elastic_lattice_deformer.engine_created)
        PHYSBAM_FATAL_ERROR("Delete_Dual_Point_Constraint() called without first creating a model");
    PRIMARY_ELASTICITY& discretization=elastic_lattice_deformer.Discretization();

    UVs.clear();
    
    pthread_mutex_lock(&elastic_lattice_deformer.simulate_lock);
    for(auto const& iter : ids ){
        PRIMARY_DEFORMER::DUAL_POINT_CONSTRAINT_MAP::iterator iter2 = elastic_lattice_deformer.dual_point_constraint_map.find( iter );
        std::array<float,3> triA_uv = iter2->second.triA_uv;
        std::array<float,3> triB_uv = iter2->second.triB_uv;
        UVs.push_back( std::array<float,6>( {{ triA_uv[0], triA_uv[1], triA_uv[2], triB_uv[0], triB_uv[1], triB_uv[2] }} ));        
    }   
    pthread_mutex_unlock(&elastic_lattice_deformer.simulate_lock);

}

//#####################################################################
// Function Get_Dual_Point_Constraint_Triangles
//#####################################################################
void CLElib::
Get_Dual_Point_Constraint_Triangles(const std::vector<int>& ids, std::vector< std::array<int,2> >& triangles){
#ifndef ENABLE_LOG_MESSAGES
    LOG::SCOPE scope("CLElib::Get_Dual_Point_Constraint_Triangles");
#endif

    
    typedef VECTOR<T,d> TV;
    typedef VECTOR<int,d> T_INDEX;

    PHYSBAM_ASSERT(implementation);
    PRIMARY_DEFORMER& elastic_lattice_deformer=Implementation<PRIMARY_DEFORMER>();
    if(!elastic_lattice_deformer.engine_created)
        PHYSBAM_FATAL_ERROR("Delete_Dual_Point_Constraint() called without first creating a model");
    PRIMARY_ELASTICITY& discretization=elastic_lattice_deformer.Discretization();
 
    triangles.clear();
    
    pthread_mutex_lock(&elastic_lattice_deformer.simulate_lock);
    for(auto const& iter : ids ){
        PRIMARY_DEFORMER::DUAL_POINT_CONSTRAINT_MAP::iterator iter2 = elastic_lattice_deformer.dual_point_constraint_map.find( iter );
        int triA = iter2->second.triangleA;
        int triB = iter2->second.triangleB;
        triangles.push_back( std::array<int,2>( {{ triA, triB }} ) );        
    }   
    pthread_mutex_unlock(&elastic_lattice_deformer.simulate_lock);

}

//#####################################################################
// Function Get_Active_Dual_Point_Constraints
//#####################################################################
void CLElib::
Get_Active_Dual_Point_Constraints(std::vector<int>& ids){
#ifndef ENABLE_LOG_MESSAGES
    LOG::SCOPE scope("CLElib::Get_Active_Dual_Point_Constraints");
#endif

    
    typedef VECTOR<T,d> TV;
    typedef VECTOR<int,d> T_INDEX;

    PHYSBAM_ASSERT(implementation);
    PRIMARY_DEFORMER& elastic_lattice_deformer=Implementation<PRIMARY_DEFORMER>();
    if(!elastic_lattice_deformer.engine_created)
        PHYSBAM_FATAL_ERROR("Delete_Dual_Point_Constraint() called without first creating a model");
    PRIMARY_ELASTICITY& discretization=elastic_lattice_deformer.Discretization();

    
    ids.clear();
    
    for( auto const& iter : elastic_lattice_deformer.dual_point_constraint_map )
        ids.push_back( iter.first );

}

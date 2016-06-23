//#####################################################################
// Copyright 2010, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class CLElib
//#####################################################################
#ifndef __CLElib_h__
#define __CLElib_h__

#include <cassert>
#include <vector>
#include <array>

class EMBEDDED_DEFORMER_BASE{
public:
    virtual ~EMBEDDED_DEFORMER_BASE() {}
};

namespace V_DATA_TYPE {
    enum {POSITION=1, STRESS=2, STRAIN=4};
};

class CLElib{
    EMBEDDED_DEFORMER_BASE* implementation;
    long frame_counter;
    int refinement;

public:
    template<class T_IMPLEMENTATION>
    T_IMPLEMENTATION& Implementation() const
    {assert(implementation);return *dynamic_cast<T_IMPLEMENTATION*>(implementation);}
    
    typedef float T;
    static const int d=3;
   

//#####################################################################    
    CLElib();
    ~CLElib();

    // These parameters must be set before the call to Create_Model, or after Destroy_Model has been called
    void Set_Youngs_Modulus(const float youngs_modulus);
    void Set_Poissons_Ratio(const float poisson_ratio);
    
    // The size of vertices & triangles is 3x the number of vertices and triangles, respectively
    // Numbering of triangle nodes is 0-based
    // dx is the cell size of the simulation grid
    void Create_Model(const std::vector<float>& vertices,const std::vector<int>& triangles,const float dx, const int refine);
    void Add_Static_Model(const std::vector<float>& vertices,const std::vector<int>& triangles);
    void Add_Muscle_Layer(const std::vector<float>& vertices,const std::vector<int>& triangles, const std::vector<float>& fiber_direction, const float& maxstress);
    void Set_Collision_Model(const std::vector<float>& vertices,const std::vector<int>& triangles);
    void Set_Collision_Proxy(const std::vector<float>& vertices,const std::vector<int>& triangles);
    void Destroy_Model();    

    void SetTextureParameters(const std::vector<int>& dimensions, const std::vector<float>& texcoords, const std::vector<int>& textriangles, const std::vector<int>& coord_map );
    void GetTextureStress(std::vector<float>& texdata );
    
    void Set_Fixed_Geometry(const std::vector<float>& vertices, const std::vector<int>& points, const std::vector<int>& segments, const std::vector<int>& triangles);
    void Set_Fixed_Points(const std::vector<float>& vertices);
    void Set_Fixed_Segments(const std::vector<float>& vertices, const std::vector<int>& segments);
    void Set_Fixed_Triangles(const std::vector<float>& vertices, const std::vector<int>& triangles);
    void Set_Fixed_Volume(const std::vector<float>& min_corner, const std::vector<float>& max_corner);

    // Once Finalize_Initialization has been called, the model is ready for simulation
    // This function can only be called once (until we've called Destroy_Model)
    // This function needs to be called, even if no triangles are fixed, before the model can be simulated
    // (in this case, if a quasistatic simulation is used, enough hooks should be set up to prevent rigid motion
    //  before attempting to call Advance_One_Time_Step())
    void Finalize_Initialization();

    // These parameters may be adjusted anytime
    void Set_Damping(const float rayleigh_coefficient);
    void Set_Density(const float density);
    void Set_Time_Step(const float dt);
    
    // These will only affect the stiffness of hooks/sutures added after calling the Set_???_Stiffness function
    void Set_Hook_Stiffness(const float hook_stiffness);
    void Set_Suture_Stiffness(const float suture_stiffness);

    // May be called anytime, but after the call to Set_Fixed_Triangles
    int Add_Single_Point_Constraint(const int triangle_id, const std::array<float,3>& uv);
    int Add_Single_Point_Constraint(const std::array<float,3>& location);
    void Move_Single_Point_Constraint(const int hook_id,const std::array<float,3>& location);
    void Delete_Single_Point_Constraint(const int hook_id);
    void Get_Single_Point_Constraint_Position(const std::vector<int>& ids, std::vector< std::array<float,3> >& position);
    void Get_Single_Point_Constraint_Triangles(const std::vector<int>& ids, std::vector< int >& triangles);
    void Get_Active_Single_Point_Constraints(std::vector<int>& ids);
    
/*
    int Add_Suture(const int triangle_id1,const float (&weights1)[2],const int triangle_id2,const float (&weights2)[2]);
    int Add_Suture(const int triangle_id1,const std::array<float,2> &weights1,
                   const int triangle_id2,const std::array<float,2> &weights2);
    int Add_Suture(const float(&location1)[3], const float(&location2)[3]);
    int Add_Suture(std::array<float,3>& location1, std::array<float,3>& location2);
    void Delete_Suture(const int suture_id);
*/
    int Add_Dual_Point_Constraint(const int triangleA_id, const std::array<float,3>& uv_A, 
                                  const int triangleB_id, const std::array<float,3>& uv_B );
    int Add_Dual_Point_Constraint(const std::array<float,3>& locationA, const std::array<float,3>& locationB );
    void Delete_Dual_Point_Constraint( const int suture_id);
    void Get_Dual_Point_Constraint_UV( const std::vector<int>& ids, std::vector<std::array<float,6> >& UVs);
    void Get_Dual_Point_Constraint_Triangles( const std::vector<int>& ids, std::vector<std::array<int,2> >& triangles);
    void Get_Active_Dual_Point_Constraints( std::vector<int>& ids );

    // May be called anytime, but after the call to Set_Fixed_Triangles
    void Advance_One_Time_Step();
    void WaitForSolve();

    // These mirror the original vertices array, given in Create_Model
    void Update_Fine_Displacement();

    /* DEPRECATED */ void Get_Vertices(std::vector<double>& vertices) const;
    /* DEPRECATED */ void Get_Strain(std::vector<double>& strain) const;
    /* DEPRECATED */ void Get_Stress(std::vector<double>& strain) const;
    /* DEPRECATED */ void Get_Vertices(std::vector<double>& vertices, int& sinceFrame) const;
    /* DEPRECATED */ void Get_Vertex_Data( std::vector<double>& vdata, int type, int &sinceFrame) const;
    void Get_Vertices(std::vector<float>& vertices) const;
    void Get_Strain(std::vector<float>& strain) const;
    void Get_Stress(std::vector<float>& strain) const;
    void Get_Vertices(std::vector<float>& vertices, int& sinceFrame) const;
    void Get_Vertex_Data( std::vector<float>& vdata, int type, int &sinceFrame) const;


    void Update_Embedded_Surfaces();
    void Generate_UFine_MeshMap();
    void UpdateDirichletCells(const std::vector<std::vector<int> >& cells);
    void Update_Collisions();


    // Debugging
    int debug_output;
    void WriteDebug(bool andDie=false);
    void Apply_Perturbation(const float perturb_amount);

 private:
    void Schedule_Exact_Solve(int krylov_iterations, int newton_iterations,  float krylov_tolerance, float newton_tolerance, long int& frame);
   
//#####################################################################    
}; 
#endif

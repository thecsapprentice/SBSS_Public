//#####################################################################
// Copyright 2010-2012, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ELASTIC_LATTICE_DEFORMER
//#####################################################################
#ifndef __ELASTIC_LATTICE_DEFORMER__
#define __ELASTIC_LATTICE_DEFORMER__

#include <pthread.h>

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>
#include <PhysBAM_Tools/Data_Structures/HASHTABLE.h>

#include <EngineInterface/ENGINE_INTERFACE.h>
#include <EngineInterface/CONSTRAINTS.h>

#include <Cutter/VOXELIZED_REGION_GENERATOR.h>

#include <Thread_Queueing/PTHREAD_QUEUE.h>

#include "CLElib.h"
#include "COLLISION_INTERFACE.h"
#include <Common/GENERIC_CELL.h>

namespace PhysBAM{

template<class T,int d,bool enable_constraints, bool enable_muscles> class SKINNING_NONLINEAR_ELASTICITY;
template<class T_NONLINEAR_ELASTICITY> class HYBRID_NONLINEAR_ELASTICITY;
template<class T_NONLINEAR_ELASTICITY> class HYBRID_NONLINEAR_ELASTICITY_STATE;
template<class T> class TRIANGULATED_SURFACE;

template<class TV> class GEOMETRY_PARTICLES;
template<class TV> class STRUCTURE;


class ELASTIC_LATTICE_DEFORMER:public EMBEDDED_DEFORMER_BASE
{
    typedef float T;
    static const int d=3;
    typedef VECTOR<T,d> TV;
    typedef VECTOR<int,d> T_INDEX;
    typedef bool T_FLAG;
    typedef ARRAY<T_FLAG, T_INDEX> T_FLAG_ARRAY;
    typedef ARRAY<T,T_INDEX> T_SCALAR_VARIABLE;
    typedef VECTOR<T_SCALAR_VARIABLE,d> T_VECTOR_VARIABLE;   
    typedef GRID<TV> T_GRID;

 public:
    //typedef HYBRID_NONLINEAR_ELASTICITY<SKINNING_NONLINEAR_ELASTICITY<T,d,true,false> > T_DISCRETIZATION;
    //typedef SKINNING_NONLINEAR_ELASTICITY<T,d,true,false> T_DISCRETIZATION;
    //T_DISCRETIZATION* discretization;

    ENGINE_INTERFACE* discretization;

 private:
    T youngs_modulus;
    T poissons_ratio;
    T rayleigh_coefficient;
    T density;
    T dt;
//    T attachment_stiffness;
    T hook_stiffness;
    T suture_stiffness;

    T time;
 public:
    ARRAY<TV> vertices;
    ARRAY<TV> vertices_B;
    ARRAY<TV> stress;
    ARRAY<TV> stress_B;
    ARRAY<TV> strain;
    ARRAY<TV> strain_B;
    ARRAY<VECTOR<int,d> > triangles;

    ARRAY<VECTOR<int,d> > render_triangles;
    ARRAY<TV> render_vertices;
    ARRAY<T_INDEX> render_embedding_map;
    ARRAY<TV> render_embedding_weights;

    ARRAY<TV>* working_vertices;
    ARRAY<TV>* display_vertices;
    ARRAY<TV>* working_stress;
    ARRAY<TV>* display_stress;
    ARRAY<TV>* working_strain;
    ARRAY<TV>* display_strain;

    RANGE< VECTOR<int,2> > texture_dimensions;
    ARRAY< T_INDEX, VECTOR<int,2> > pixel_embedding_map;
    ARRAY< TV, VECTOR<int,2> > pixel_embedding_weights;    
    
    ARRAY< ARRAY<TV> > static_model_vertices;
    ARRAY< ARRAY<VECTOR<int,d> > > static_model_triangles;

    ARRAY<TV> collision_proxy_vertices;
    ARRAY<VECTOR<int,d> > collision_proxy_triangles;
    ARRAY<T_INDEX> collision_proxy_embedding_map;
    bool collision_proxy_set;

    ARRAY<TV> attachment_X_in_local_frame;
    ARRAY<T_INDEX> attachment_embedding_map;
    ARRAY<TV> attachment_embedding_weights;
    ARRAY<int> attachment_bone_id;
    ARRAY<TRIANGULATED_SURFACE<T>*> bone_surfaces;

    int max_muscle;
    ARRAY< ARRAY<int>, T_INDEX> muscle_id;
    ARRAY< ARRAY<int>, int> muscle_id_mesh;
    ARRAY< ARRAY<T>, T_INDEX> muscle_density;
    ARRAY< ARRAY<T>, int> muscle_density_mesh;
    ARRAY< ARRAY<TV>, T_INDEX> muscle_fiber;
    ARRAY< ARRAY<TV>, int> muscle_fiber_mesh;
    

    ARRAY<TRIANGULATED_SURFACE<T>*> embedded_surfaces;
    ARRAY<ARRAY<T_INDEX> > embedded_surfaces_map;
    ARRAY<ARRAY<TV> > embedded_surfaces_weights;

    ARRAY<T_INDEX> embedding_map;
    ARRAY<TV> embedding_weights;

    ARRAY<CONSTRAINT_SEGMENT<T,d> > fine_point_constraints;
    ARRAY<int> constraint_index_of_hook_id;
    ARRAY<int> constraint_index_of_suture_id;
    ARRAY<int> fine_to_coarse_constraint;

    T_FLAG_ARRAY voxmap;
    T_FLAG_ARRAY voxmap_node;
    T_FLAG_ARRAY voxmap_dirichlet;
    T_VECTOR_VARIABLE u_fine;
    T_VECTOR_VARIABLE u_fine_stress;
    T_VECTOR_VARIABLE u_fine_strain;
    T_GRID coarse_grid;
    T_GRID fine_grid;
    RANGE<T_INDEX> unpadded_fine_domain;
    RANGE<T_INDEX> padded_fine_domain;
    RANGE<T_INDEX> unpadded_fine_node_domain;
    RANGE<T_INDEX> padded_fine_node_domain;
    T h;
    T coarse_h;
    int refinement;
    T_SCALAR_VARIABLE fine_to_coarsemesh;

    VOXELIZED_REGION_GENERATOR<T,d>* vrg;
    VECTOR<HASHTABLE<T_INDEX,GENERIC_CELL<T,d> >,2> cell_is_mesh_from;

    bool engine_created;
    bool solver_initialized;
    bool domain_initialized;

    int frame;

    PTHREAD_QUEUE* simulation_queue;
    pthread_mutex_t simulate_lock;
    mutable pthread_mutex_t geometry_buffer_lock;

 protected:
    COLLISION_INTERFACE<T,d>* collision_shape;
    

public:
    ELASTIC_LATTICE_DEFORMER();
    ~ELASTIC_LATTICE_DEFORMER();

    ENGINE_INTERFACE& Discretization()
    {return *discretization;}

    const ENGINE_INTERFACE& ConstDiscretization() const
    {return *discretization;}

    void Destroy_Engine(){
        simulation_queue->Wait();
        Discretization().DestroyEngine();
        engine_created = false;
        delete vrg;
        vrg=NULL;
        fine_point_constraints.Clean_Memory();
        constraint_index_of_hook_id.Clean_Memory();
        constraint_index_of_suture_id.Clean_Memory();
        fine_to_coarse_constraint.Clean_Memory();
    }

    void Initialize_Fine_Domain( T_INDEX cells, T dx, TV min_corner, int ratio,
                                 int lcellpad, int ucellpad, int lnodepad,int unodepad)
    {
        max_muscle = 0;
        refinement = ratio;
        h = dx;
        fine_grid.Initialize(cells+1,RANGE<TV>(min_corner,min_corner+TV(cells)*h) );
        unpadded_fine_domain = RANGE<T_INDEX>(fine_grid.Cell_Indices());
        padded_fine_domain = 
        RANGE<T_INDEX>(unpadded_fine_domain.min_corner-lcellpad*ratio,
                       unpadded_fine_domain.max_corner+ucellpad*ratio);
        unpadded_fine_node_domain = RANGE<T_INDEX>(fine_grid.Node_Indices());
        padded_fine_node_domain = 
        RANGE<T_INDEX>(unpadded_fine_node_domain.min_corner-lnodepad*ratio,
                       unpadded_fine_node_domain.max_corner+unodepad*ratio);
        voxmap.Resize(padded_fine_domain); voxmap.Fill(false);
        voxmap_node.Resize(padded_fine_node_domain); voxmap.Fill(false);
        voxmap_dirichlet.Resize(padded_fine_domain); voxmap.Fill(false);
        for(int v=1;v<=3;v++){ u_fine(v).Resize(padded_fine_domain); u_fine(v).Fill(T()); }
        for(int v=1;v<=3;v++){ u_fine_strain(v).Resize(padded_fine_domain); u_fine_strain(v).Fill(T()); }
        for(int v=1;v<=3;v++){ u_fine_stress(v).Resize(padded_fine_domain); u_fine_stress(v).Fill(T()); }
        fine_to_coarsemesh.Resize(unpadded_fine_domain, true, false, 0); fine_to_coarsemesh.Fill(0);
    }

    int Number_Of_Bones() const
    {return bone_surfaces.m;}
    
    const TRIANGULATED_SURFACE<T>& Bone_Surface(int bone_id) const
    {return *bone_surfaces(bone_id);}
    
    int Number_Of_Embedded_Surfaces() const
    {return embedded_surfaces.m;}
    
    const TRIANGULATED_SURFACE<T>& Embedded_Surface(int surface_id) const
    {return *embedded_surfaces(surface_id);}
    
    const ARRAY<TV>& Vertices() const
    {return vertices;}
    
    const ARRAY<VECTOR<int,d> >& Triangles() const
    {return triangles;}

    T Current_Time() const
    {return time;}
    
    int Add_Embedded_Point_To_Fixed_Point_Spring_Constraint(const T spring_coefficient,const TV& embedded_point_material_space_location,const TV& fixed_point_world_space_location);
    int Add_Two_Embedded_Point_Spring_Constraint(const T spring_coefficient,const TV& embedded_point_material_space_location1,const TV& embedded_point_material_space_location2);

    int Constraint_Count() {return fine_point_constraints.m; };
    void Initialize_Solver();

    // Pressure routines for CG
    float Exact_Solve(const int krylov_iterations,const int newton_iterations,const T krylov_tolerance,const T newton_tolerance,const bool no_cut_cells);
    void Print_Volume_Change_Diagnostics();
    void Print_Force_Diagnostics(const bool no_cut_cells);
    void Print_Matrix_Condition(const bool no_cut_cells);
    void Material_Convergence_Rate(const bool no_cut_cells);

//    TV Multilinear_Interpolation(const T_INDEX& cell_index,const TV& weights);  // We need these for embedded objects, but where do the objects live? Here or in Nonlinear?
//    void Interpolate_Embedded_Object();

    void Add_Discretization_Constraint(int cid);
    void Update_Discretization_Constraint(int cid);
    void Remove_Discretization_Constraint(int cid, int other_cid);
    void Initialize_Collision_Points();

    void SwapGeometryBuffers(){
        pthread_mutex_lock(&geometry_buffer_lock);
        ARRAY<TV>* temp;

        // Swap Positional Data
        temp = working_vertices;
        working_vertices = display_vertices;
        display_vertices = temp;

        // Swap Stress Data
        temp = working_stress;
        working_stress = display_stress;
        display_stress = temp;
     
        // Swap Strain Data
        temp = working_strain;
        working_strain = display_strain;
        display_strain = temp;

        pthread_mutex_unlock(&geometry_buffer_lock);
    }   
 
    void Set_Collision_Object(COLLISION_INTERFACE<T,d>* new_collision_object){
        if( collision_shape )
            delete collision_shape;
        
        collision_shape = new_collision_object;        
    };

    void Update_Muscles( const ARRAY<int, T_INDEX>& ids,  const ARRAY<int, int>& ids_mesh, 
                         const ARRAY<T, T_INDEX>& density,  const ARRAY<T, int>& density_mesh,
                         const ARRAY<TV, T_INDEX>& fiber,  const ARRAY<TV, int>& fiber_mesh,
                         float maxstress);
    

    friend class ::CLElib;
//#####################################################################    
//#####################################################################    
}; 
}
#endif

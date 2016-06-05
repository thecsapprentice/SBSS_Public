//#####################################################################
// Copyright 2015, Raj Setaluri.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class FRACTURE_HELPER
//#####################################################################
#ifndef __FRACTURE_HELPER__
#define __FRACTURE_HELPER__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TETRAHEDRALIZED_VOLUME.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <PhysBAM_Geometry/Topology/POLYGON_MESH.h>
#include <PhysBAM_Tools/Read_Write/Arrays/READ_WRITE_ARRAY.h>
#include <PhysBAM_Tools/Read_Write/Vectors/READ_WRITE_VECTOR.h>
#include <PhysBAM_Tools/Read_Write/Data_Structures/READ_WRITE_HASHTABLE.h>
#include <PhysBAM_Geometry/Read_Write/Topology/READ_WRITE_POLYGON_MESH.h>
#include <PhysBAM_Geometry/Read_Write/Topology/READ_WRITE_SIMPLEX_MESH.h>

namespace PhysBAM{
//#####################################################################
// Struct CUTTING_SIMPLEX_RAW
//#####################################################################
template<class T,int d>
struct CUTTING_SIMPLEX_RAW
{
    enum SIMPLEX_TYPE{GLOBAL_EMBEDDING_FACE,LOCAL_EMBEDDING_FACE,GLOBAL_CUT_FACE,LOCAL_CUT_FACE};
    SIMPLEX_TYPE type;
    int parent;
    int element_owner;
    int original_index;
    VECTOR<int,d> nodes;
    VECTOR<VECTOR<T,d>,d> weights;
    T abs_tol;
    VECTOR<VECTOR<T,d>,d+1> element_original_coordinates;
    VECTOR<VECTOR<T,d>,d> simplex_original_coordinates;
    VECTOR<bool,d> node_in_embedded_simplex;
};
template<class RW,class T,int d>
class Read_Write<CUTTING_SIMPLEX_RAW<T,d>,RW>
{
public:
    static void Read(std::istream& input,CUTTING_SIMPLEX_RAW<T,d>& object)
    {Read_Binary<RW>(input,object.type,object.parent,object.element_owner,object.nodes,object.node_in_embedded_simplex);Read_Binary<T>(input,object.weights,object.abs_tol,object.element_original_coordinates,object.simplex_original_coordinates);}
    static void Write(std::ostream& output,const CUTTING_SIMPLEX_RAW<T,d>& object)
    {Write_Binary<RW>(output,object.type,object.parent,object.element_owner,object.nodes,object.node_in_embedded_simplex);Write_Binary<T>(output,object.weights,object.abs_tol,object.element_original_coordinates,object.simplex_original_coordinates);}
};
//#####################################################################
// Struct CUTTING_SIMPLICES_RAW
//#####################################################################
template<class T,int d>
struct CUTTING_SIMPLICES_RAW
{
    ARRAY<CUTTING_SIMPLEX_RAW<T,d> > simplices;
    int index_for_last_old_cutting_simplex;
};
template<class RW,class T,int d>
class Read_Write<CUTTING_SIMPLICES_RAW<T,d>,RW>
{
public:
    static void Read(std::istream& input,CUTTING_SIMPLICES_RAW<T,d>& object)
    {Read_Binary<RW>(input,object.simplices,object.index_for_last_old_cutting_simplex);}
    static void Write(std::ostream& output,const CUTTING_SIMPLICES_RAW<T,d>& object)
    {Write_Binary<RW>(output,object.simplices,object.index_for_last_old_cutting_simplex);}
};
//#####################################################################
// Struct INTERSECTION_REGISTRY_RAW
//#####################################################################
template<class T,int d>
struct INTERSECTION_REGISTRY_RAW
{
    CUTTING_SIMPLICES_RAW<T,d>& cutting_simplices;
    ARRAY<ARRAY<int> > intersections_on_simplex;
    ARRAY<ARRAY<int> > simplices_on_intersection;
    ARRAY<ARRAY<VECTOR<T,d-1> > > simplex_weights_on_intersection;
    int index_for_last_old_intersection;

    INTERSECTION_REGISTRY_RAW(CUTTING_SIMPLICES_RAW<T,d>& cutting_simplices_input):cutting_simplices(cutting_simplices_input) {}
};
template<class RW,class T,int d>
class Read_Write<INTERSECTION_REGISTRY_RAW<T,d>,RW>
{
public:
    static void Read(std::istream& input,INTERSECTION_REGISTRY_RAW<T,d>& object)
    {Read_Binary<RW>(input,object.cutting_simplices,object.intersections_on_simplex,object.simplices_on_intersection,object.simplex_weights_on_intersection,object.index_for_last_old_intersection);}
    static void Write(std::ostream& output,const INTERSECTION_REGISTRY_RAW<T,d>& object)
    {Write_Binary<RW>(output,object.cutting_simplices,object.intersections_on_simplex,object.simplices_on_intersection,object.simplex_weights_on_intersection,object.index_for_last_old_intersection);}
};
//#####################################################################
// Struct CUTTING_PARTICLES_RAW
//#####################################################################
struct CUTTING_PARTICLES_RAW
{
    enum CUTTING_PARTICLE_ID_TYPE{TET_NODE_ID,INTERSECTION_ID,TET_NODE_AND_INTERSECTION_ID};
    ARRAY<int> tet_node_indices;
    ARRAY<int> intersection_indices;
    ARRAY<CUTTING_PARTICLE_ID_TYPE> particle_ids_types;
    HASHTABLE<int,int> intersection_to_particle_id;
    HASHTABLE<int,int> tet_node_to_particle_id;
};
template<class RW>
class Read_Write<CUTTING_PARTICLES_RAW,RW>
{
public:
    static void Read(std::istream& input,CUTTING_PARTICLES_RAW& object)
    {Read_Binary<RW>(input,object.tet_node_indices,object.intersection_indices,object.particle_ids_types,object.intersection_to_particle_id,object.tet_node_to_particle_id);}
    static void Write(std::ostream& output,const CUTTING_PARTICLES_RAW& object)
    {Write_Binary<RW>(output,object.tet_node_indices,object.intersection_indices,object.particle_ids_types,object.intersection_to_particle_id,object.tet_node_to_particle_id);}
};
//#####################################################################
// Struct CUTTING_POLYGON_RAW
//#####################################################################
struct CUTTING_POLYGON_RAW
{
    enum POLYGON_TYPE{FACE_BOUNDARY,FACE_INTERIOR,TRIANGLE_CLIPPED};
    int polygon_index; // indexes into polygon_mesh
    int simplex_owner; // indexes into cutting_simplices
    bool flipped;
    POLYGON_TYPE polygon_type;
};
template<class RW>
class Read_Write<CUTTING_POLYGON_RAW,RW>
{
public:
    static void Read(std::istream& input,CUTTING_POLYGON_RAW& object)
    {Read_Binary<RW>(input,object.polygon_index,object.simplex_owner,object.flipped,object.polygon_type);}
    static void Write(std::ostream& output,const CUTTING_POLYGON_RAW& object)
    {Write_Binary<RW>(output,object.polygon_index,object.simplex_owner,object.flipped,object.polygon_type);}
};
//#####################################################################
// Class FRACTURE_HELPER
//#####################################################################
template<class T>
class FRACTURE_HELPER
{
public:
    enum{d=3};
    typedef VECTOR<T,d> TV;
    typedef VECTOR<int,d> TV_INT;

public:
    // CURRENT STATE
    TETRAHEDRALIZED_VOLUME<T>* current_tetrahedralized_volume;

    // VALID ONLY FOR THE CURRENT CUT ITERATION
    TRIANGULATED_SURFACE<T>* cutting_triangulated_surface;
    //boost::scoped_ptr<TETRAHEDRALIZED_VOLUME<T> > next_tetrahedralized_volume;
    TRIANGLE_MESH final_duplicated_boundary_mesh;
    ARRAY<ARRAY<int> > final_parents_per_new_particle;
    ARRAY<ARRAY<T> > final_parent_weights_per_new_particle;
    ARRAY<int> old_particle_per_new_collapsed_particle;
    //ARRAY<ARRAY<int> > new_collapsed_tet_particle_per_old_tet_particle;

    // replaced at the end of current iteration
    CUTTING_SIMPLICES_RAW<T,3>* cutting_simplices; // TODO: try to make an instance
    INTERSECTION_REGISTRY_RAW<T,3>* intersection_registry; // TODO: try to make an instance
    POLYGON_MESH polygon_mesh;
    ARRAY<ARRAY<int> > simplices_per_current_tet;
    ARRAY<CUTTING_POLYGON_RAW> current_cutting_polygons;
    ARRAY<ARRAY<ARRAY<int> > > regions_per_tet;
    HASHTABLE<VECTOR<int,3>,VECTOR<int,3> > sorted_to_original_boundary_nodes;
    CUTTING_PARTICLES_RAW cutting_particles;

    // overwritten at end of current iteration and used in the next
    ARRAY<ARRAY<int> > cutting_polygons_per_cutting_simplex;

    // must be refreshed when first used in the current iteration
    //int intersection_counter;
    //ARRAY<int> last_old_simplex_index_per_current_tet;
    //ARRAY<ARRAY<int> > intersecting_simplices_per_simplex; // only records old-new and new-new intersections
    //ARRAY<ARRAY<int> > polygons_per_element;
    //HASHTABLE<PAIR<int,int>,int> hash_all_new_uncollapsed_particles; // <dup tet, orig particle> -> dup particle
    //ARRAY<ARRAY<int> > uncollapsed_new_particles_per_all_current_particle_ids;
    //ARRAY<ARRAY<int> > new_parents_per_new_particle;
    //ARRAY<ARRAY<T> > new_parent_weights_per_new_particle;
    //ARRAY<int> current_particle_id_per_uncollapsed_new_particle;
    //ARRAY<ARRAY<int> > new_tets_per_current;
    //ARRAY<int> current_particle_id_per_collapsed_new_particle;
    //int num_new_tets,num_new_particles;
    //UNION_FIND<> union_vertices;
    //ARRAY<int> new_particle_indices;
    //HASHTABLE<int,int> dup_tet_before_to_after_collapse;
    //HASHTABLE<int,int> dup_tet_after_to_before_collapse;
    // regions_per_otet(otet index)(dup tet local index)(region local index)(cutting polygon index)
    //ARRAY<ARRAY<ARRAY<ARRAY<int> > > > regions_per_otet;

    FRACTURE_HELPER()
    {current_tetrahedralized_volume=TETRAHEDRALIZED_VOLUME<T>::Create();
    cutting_triangulated_surface=TRIANGULATED_SURFACE<T>::Create();}

    ~FRACTURE_HELPER()
    {} // TODO: delete ptrs?

    void Read_From_File(STREAM_TYPE stream_type,const std::string prefix)
    {
        cutting_simplices=new CUTTING_SIMPLICES_RAW<T,d>();
        intersection_registry=new INTERSECTION_REGISTRY_RAW<T,d>(*cutting_simplices);
#if 0 // TODO: fix
        FILE_UTILITIES::Read_From_File(stream_type,prefix+"tet_volume",*current_tetrahedralized_volume);
#else
        Read_Tetrahedralized_Volume(stream_type,prefix+"tet_volume",*current_tetrahedralized_volume);
#endif
        FILE_UTILITIES::Read_From_File(stream_type,prefix+"cutting_simplices",*cutting_simplices);
        FILE_UTILITIES::Read_From_File(stream_type,prefix+"intersection_registry",*intersection_registry);
        FILE_UTILITIES::Read_From_File(stream_type,prefix+"simplices_per_tet",simplices_per_current_tet);
        FILE_UTILITIES::Read_From_File(stream_type,prefix+"cutting_polygons",current_cutting_polygons);
        FILE_UTILITIES::Read_From_File(stream_type,prefix+"regions_per_tet",regions_per_tet);
        FILE_UTILITIES::Read_From_File(stream_type,prefix+"boundary_nodes",sorted_to_original_boundary_nodes);
        FILE_UTILITIES::Read_From_File(stream_type,prefix+"particle_ids",cutting_particles);
        FILE_UTILITIES::Read_From_File(stream_type,prefix+"polygons_per_simplex",cutting_polygons_per_cutting_simplex);
        FILE_UTILITIES::Read_From_File(stream_type,prefix+"polygon_mesh",polygon_mesh);
        // for simulation (in addition to tet volume)
#if 0 // TODO: fix
        FILE_UTILITIES::Read_From_File(stream_type,prefix+"final_cutting_surface",*cutting_triangulated_surface);
#else
        Read_Triangulated_Surface(stream_type,prefix+"final_cutting_surface",*cutting_triangulated_surface);
#endif
        FILE_UTILITIES::Read_From_File(stream_type,prefix+"final_boundary_mesh",final_duplicated_boundary_mesh);
        FILE_UTILITIES::Read_From_File(stream_type,prefix+"final_parents",final_parents_per_new_particle);
        FILE_UTILITIES::Read_From_File(stream_type,prefix+"final_weights",final_parent_weights_per_new_particle);
        FILE_UTILITIES::Read_From_File(stream_type,prefix+"old_particle_per_new_particle",old_particle_per_new_collapsed_particle);
    }

private:
    void Read_Tetrahedralized_Volume(STREAM_TYPE stream_type,const std::string prefix,TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume)
    {
        ARRAY<TV> X;
        FILE_UTILITIES::Read_From_File(stream_type,prefix+"_X",X);
        ARRAY<VECTOR<int,4> > elements;
        FILE_UTILITIES::Read_From_File(stream_type,prefix+"_elements",elements);
        for(int p=1;p<=X.m;p++)
            tetrahedralized_volume.particles.X(tetrahedralized_volume.particles.array_collection->Add_Element())=X(p);
        for(int e=1;e<=elements.m;e++)
            tetrahedralized_volume.mesh.elements.Append(elements(e));
        tetrahedralized_volume.Update_Number_Nodes();
    }

    void Read_Triangulated_Surface(STREAM_TYPE stream_type,const std::string prefix,TRIANGULATED_SURFACE<T>& triangulated_surface)
    {
        ARRAY<TV> X;
        FILE_UTILITIES::Read_From_File(stream_type,prefix+"_X",X);
        ARRAY<VECTOR<int,3> > elements;
        FILE_UTILITIES::Read_From_File(stream_type,prefix+"_elements",elements);
        for(int p=1;p<=X.m;p++)
            triangulated_surface.particles.X(triangulated_surface.particles.array_collection->Add_Element())=X(p);
        for(int e=1;e<=elements.m;e++)
            triangulated_surface.mesh.elements.Append(elements(e));
        triangulated_surface.Update_Number_Nodes();
    }
};
}
#endif

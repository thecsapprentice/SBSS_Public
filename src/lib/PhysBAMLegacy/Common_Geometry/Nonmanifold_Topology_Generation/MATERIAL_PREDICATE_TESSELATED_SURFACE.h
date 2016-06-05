#ifndef __MATERIAL_PREDICATE_TESSELATED_SURFACE_H__
#define __MATERIAL_PREDICATE_TESSELATED_SURFACE_H__

#include <Common_Geometry/Nonmanifold_Topology_Generation/MATERIAL_PREDICATE.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TETRAHEDRALIZED_VOLUME.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_AREA.h>
#include <PhysBAM_Geometry/Basic_Geometry/TRIANGLE_2D.h>
#include <PhysBAM_Geometry/Basic_Geometry/TETRAHEDRON.h>
#include <PhysBAM_Geometry/Basic_Geometry/TRIANGLE_3D.h>
#include <PhysBAM_Geometry/Basic_Geometry/SEGMENT_2D.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/SEGMENTED_CURVE_2D.h>

namespace PhysBAM{

template<class T, int d> class MATERIAL_PREDICATE_TESSELATED_SURFACE_POLICY;

template<class T>
struct MATERIAL_PREDICATE_TESSELATED_SURFACE_POLICY<T,1>{
    typedef SEGMENTED_CURVE_2D<T> T_VOLUME;
    typedef SEGMENT_2D<T> T_ELEMENT;
};

template<class T>
struct MATERIAL_PREDICATE_TESSELATED_SURFACE_POLICY<T,2>{
    typedef TRIANGULATED_AREA<T> T_VOLUME;
    typedef SEGMENTED_CURVE_2D<T> T_BOUNDARY;
    typedef TRIANGLE_2D<T> T_ELEMENT;
    typedef SEGMENT_2D<T> T_BOUNDARY_ELEMENT;
};

template<class T>
struct MATERIAL_PREDICATE_TESSELATED_SURFACE_POLICY<T,3>{
    typedef TETRAHEDRALIZED_VOLUME<T> T_VOLUME;       
    typedef TRIANGULATED_SURFACE<T> T_BOUNDARY;
    typedef TETRAHEDRON<T> T_ELEMENT;
    typedef TRIANGLE_3D<T> T_BOUNDARY_ELEMENT;
};

template<class T, int d>
class MATERIAL_PREDICATE_TESSELATED_SURFACE: public MATERIAL_PREDICATE<T,d>{
public:
    typedef VECTOR<T,d> TV;
    typedef MATERIAL_PREDICATE<T,d> BASE;
    typedef typename BASE::T_CELL T_CELL;
    typedef typename BASE::T_MESH T_MESH;
    typedef typename MATERIAL_PREDICATE_TESSELATED_SURFACE_POLICY<T,d>::T_VOLUME T_VOLUME;
    typedef typename MATERIAL_PREDICATE_TESSELATED_SURFACE_POLICY<T,d>::T_BOUNDARY T_BOUNDARY;
    typedef typename MATERIAL_PREDICATE_TESSELATED_SURFACE_POLICY<T,d>::T_ELEMENT T_ELEMENT;
    typedef typename MATERIAL_PREDICATE_TESSELATED_SURFACE_POLICY<T,d>::T_BOUNDARY_ELEMENT T_BOUNDARY_ELEMENT;

    const T_VOLUME& current_tetrahedralized_volume;
    const ARRAY<VECTOR<bool,4> >& material_nodes_per_duplicate_tet;
    const ARRAY<int>& old_particle_per_new_collapsed_particle;
    const HASHTABLE<VECTOR<int,2>,ARRAY<VECTOR<int,2> > >& template_hex_node_to_template_tet_nodes;
    const ARRAY<int>& duplicate_tet_to_parent_tet;
    const ARRAY<ARRAY<int> >& template_hex_to_template_tets;
    const ARRAY<ARRAY<int> >& template_tet_to_duplicate_tets;
    const ARRAY<T>& duplicate_tet_nodal_distances;

    ARRAY<ARRAY<int> > material_fragments;
    ARRAY<HASHTABLE<int> > material_fragments_hash;
    HASHTABLE<int,int> duplicate_cell_to_root_cell;

    MATERIAL_PREDICATE_TESSELATED_SURFACE(const T_VOLUME& current_tetrahedralized_volume_input,const ARRAY<VECTOR<bool,4> >& material_nodes_per_duplicate_tet_input,const ARRAY<int>& old_particle_per_new_collapsed_particle_input,const HASHTABLE<VECTOR<int,2>,ARRAY<VECTOR<int,2> > >& template_hex_node_to_template_tet_nodes_input,const ARRAY<int>& duplicate_tet_to_parent_tet_input,const ARRAY<ARRAY<int> >& template_hex_to_template_tets_input,const ARRAY<ARRAY<int> >& template_tet_to_duplicate_tets_input,const ARRAY<T>& duplicate_tet_nodal_distances_input);
    ~MATERIAL_PREDICATE_TESSELATED_SURFACE();

    virtual void MaterialFragments(const T_MESH& mesh,const HASHTABLE<PAIR<VECTOR<int,2>,int> >& linkage_list,ARRAY<T_CELL>& sub_cells,HASHTABLE<int,int>& subcell_to_root_cell);
    virtual bool IsMaterialContinuous(const int axis,const int root_cell_high, const int sub_cell_high,const int root_cell_low,const int sub_cell_low) const;
    virtual void MergeSubcellMaterial(const UNION_FIND<int>& merge_map);
    virtual void ComputeNodalDistances(const ARRAY<TV>& nodal_positions,int subcell,ARRAY<T>& distances) const;
    virtual bool InsideMaterial(const int subcell_index,const int node_index) const;
    virtual int Node_Material_Representative(const int subcell_index,const int node_index) const;
    virtual int Subcell_Material_Representative(const int subcell_index) const;
};
}
#endif

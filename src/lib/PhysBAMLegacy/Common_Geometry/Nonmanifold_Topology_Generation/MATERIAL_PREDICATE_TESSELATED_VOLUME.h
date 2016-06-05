#ifndef __MATERIAL_PREDICATE_TESSELATED_VOLUME_H__
#define __MATERIAL_PREDICATE_TESSELATED_VOLUME_H__

#include <Common_Geometry/Nonmanifold_Topology_Generation/CUTTER.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TETRAHEDRALIZED_VOLUME.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_AREA.h>
#include <PhysBAM_Geometry/Basic_Geometry/TRIANGLE_2D.h>
#include <PhysBAM_Geometry/Basic_Geometry/TETRAHEDRON.h>
#include <PhysBAM_Geometry/Basic_Geometry/TRIANGLE_3D.h>
#include <PhysBAM_Geometry/Basic_Geometry/SEGMENT_2D.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/SEGMENTED_CURVE_2D.h>

namespace PhysBAM{

    template<class T, int d> class MATERIAL_PREDICATE_TESSELATED_VOLUME_POLICY;
   
    template<class T>
    class MATERIAL_PREDICATE_TESSELATED_VOLUME_POLICY<T,1>{
    public:
        typedef SEGMENTED_CURVE_2D<T> T_VOLUME;
        typedef SEGMENT_2D<T> T_ELEMENT;
    };
 
    template<class T>
    class MATERIAL_PREDICATE_TESSELATED_VOLUME_POLICY<T,2>{
    public:
        typedef TRIANGULATED_AREA<T> T_VOLUME;
        typedef SEGMENTED_CURVE_2D<T> T_BOUNDARY;
        typedef TRIANGLE_2D<T> T_ELEMENT;
        typedef SEGMENT_2D<T> T_BOUNDARY_ELEMENT;

    };

    template<class T>
    class MATERIAL_PREDICATE_TESSELATED_VOLUME_POLICY<T,3>{
    public:
        typedef TETRAHEDRALIZED_VOLUME<T> T_VOLUME;       
        typedef TRIANGULATED_SURFACE<T> T_BOUNDARY;
        typedef TETRAHEDRON<T> T_ELEMENT;
        typedef TRIANGLE_3D<T> T_BOUNDARY_ELEMENT;

    };



    template<class T, int d>
    class MATERIAL_PREDICATE_TESSELATED_VOLUME : public MATERIAL_PREDICATE<T,d> {
    public:
        typedef VECTOR<T,d> TV;
        typedef MATERIAL_PREDICATE<T,d> BASE;
        typedef typename BASE::T_CELL T_CELL;
        typedef typename BASE::T_MESH T_MESH;
        typedef typename MATERIAL_PREDICATE_TESSELATED_VOLUME_POLICY<T,d>::T_VOLUME T_VOLUME;
        typedef typename MATERIAL_PREDICATE_TESSELATED_VOLUME_POLICY<T,d>::T_BOUNDARY T_BOUNDARY;
        typedef typename MATERIAL_PREDICATE_TESSELATED_VOLUME_POLICY<T,d>::T_ELEMENT T_ELEMENT;
        typedef typename MATERIAL_PREDICATE_TESSELATED_VOLUME_POLICY<T,d>::T_BOUNDARY_ELEMENT T_BOUNDARY_ELEMENT;

        MATERIAL_PREDICATE_TESSELATED_VOLUME(T_VOLUME& volume );
        
        virtual void MaterialFragments(const T_MESH& mesh,
                                       const HASHTABLE<PAIR<VECTOR<int,2>,int> >& linkage_list,
                                       ARRAY<T_CELL>& sub_cells, 
                                       HASHTABLE<int, int>& subcell_to_root_cell);

        virtual bool IsMaterialContinuous(const int axis,
                                          const int root_cell_high, const int sub_cell_high,
                                          const int root_cell_low,  const int sub_cell_low) const;




        virtual void MergeSubcellMaterial( const UNION_FIND<int>& merge_map );

        virtual void ComputeNodalDistances( const ARRAY<TV>& nodal_positions, int subcell, ARRAY<T>& distances) const;
        virtual bool InsideMaterial(const int subcell_index,const int node_index) const;
        virtual int Node_Material_Representative(const int subcell_index,const int node_index) const;
        virtual int Subcell_Material_Representative(const int subcell_index) const;

        T_VOLUME& volume;
        T_BOUNDARY* boundary;
        BOX_HIERARCHY<TV> aabb_hierarchy;
        //BOX_HIERARCHY<TV> aabb_boundary_hierarchy;

        ARRAY< ARRAY<int> > material_fragments;
        ARRAY< VECTOR<ARRAY<int>, T_MESH::vertices_per_cell > > corner_fragments;
        ARRAY< ARRAY<int> > boundary_fragments;
    };



}

#endif

#ifndef __MATERIAL_PREDICATE_VOXELIZED_VOLUME_H__
#define __MATERIAL_PREDICATE_VOXELIZED_VOLUME_H__

#include <Common_Geometry/Nonmanifold_Topology_Generation/CUTTER.h>

namespace PhysBAM{

    template<class T, int d> class MATERIAL_PREDICATE_VOXELIZED_VOLUME_POLICY;
   
    template<class T>
    class MATERIAL_PREDICATE_VOXELIZED_VOLUME_POLICY<T,1>{
    public:
    };
 
    template<class T>
    class MATERIAL_PREDICATE_VOXELIZED_VOLUME_POLICY<T,2>{
    public:
    };

    template<class T>
    class MATERIAL_PREDICATE_VOXELIZED_VOLUME_POLICY<T,3>{
    public:
    };
    

    template<class T, int d>
    class MATERIAL_PREDICATE_VOXELIZED_VOLUME : public MATERIAL_PREDICATE<T,d> {
    public:
        typedef VECTOR<T,d> TV;
        typedef VECTOR<int,d> T_INDEX;
        typedef MATERIAL_PREDICATE<T,d> BASE;
        typedef typename BASE::T_CELL T_CELL;
        typedef typename BASE::T_MESH T_MESH;
        typedef ARRAY< bool, T_INDEX > VOXMAP;

        MATERIAL_PREDICATE_VOXELIZED_VOLUME(const int refinement_input,
                                            const ARRAY< VOXMAP > voxmap_input );
        
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

        const VOXMAP& SubcellVoxMaterialFragment(const int subcell_index ) const;

    private:
        const int refinement;
        const ARRAY< VOXMAP > voxmap;
        const RANGE<T_INDEX> voxmap_domain;

        ARRAY< VOXMAP > material_fragments;
        ARRAY< VECTOR< T_INDEX , T_MESH::vertices_per_cell > > corner_fragments;

        int VoxelIndexToFlat( const T_INDEX voxel_index ) const ;
        T_INDEX FlatIndexToVoxel( const int flat_index ) const ;

    };



}

#endif

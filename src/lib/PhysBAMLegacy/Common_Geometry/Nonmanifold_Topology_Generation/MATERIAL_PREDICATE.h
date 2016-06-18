#ifndef __MATERIAL_PREDICATE_H__
#define __MATERIAL_PREDICATE_H__

#include <Common_Geometry/Topology/REGULAR_HYPERCUBE_MESH.h>
#include <PhysBAM_Tools/Data_Structures/UNION_FIND.h>
#include <PhysBAM_Tools/Data_Structures/HASHTABLE.h>
#include <PhysBAM_Tools/Data_Structures/PAIR.h>
#include <stdint.h>

#include <memory>
#if __cplusplus <= 199711L
// Fall back to using old auto_ptr
#define MP_AUTO_PTR_TYPE std::auto_ptr
#else
// Use the new C++11 unique_ptr
#define MP_AUTO_PTR_TYPE std::unique_ptr
#endif


namespace PhysBAM
{

template<class T,int d>
class MATERIAL_PREDICATE{        
public:
    typedef VECTOR<T,d> TV;
    typedef REGULAR_HYPERCUBE_MESH<T,d> T_MESH;
    typedef typename REGULAR_HYPERCUBE_MESH<T,d>::T_CELL T_CELL;

    virtual std::string Name() {return std::string("GENERIC MATERIAL PREDICATE"); };
    virtual void MaterialFragments(const T_MESH& mesh,const HASHTABLE<PAIR<VECTOR<int,2>,int> >& linkage_list,ARRAY<T_CELL>& sub_cells,HASHTABLE<int, int>& subcell_to_root_cell)=0;
    // Axis - the face/edge axis (x,y,z...)
    // High Cells - cell whose interesting face is on the positive axis side
    // Low Cells - cell whose interesting face is on the negative axis side
    //   Note: This is the SAME face, just different perspectives of it.
    virtual bool IsMaterialContinuous(const int axis,const int root_cell_high,const int sub_cell_high,const int root_cell_low,const int sub_cell_low) const=0;
    virtual void MergeSubcellMaterial(const UNION_FIND<int>& merge_map)=0;
    virtual void ComputeNodalDistances(const ARRAY<TV>& nodal_positions,int subcell,ARRAY<T>& distances) const=0;
    virtual bool InsideMaterial(const int subcell_index,const int node_index) const=0;
    virtual int Node_Material_Representative(const int subcell_index,const int node_index) const=0;
    virtual int Subcell_Material_Representative(const int subcell_index) const=0;
};
}
#endif

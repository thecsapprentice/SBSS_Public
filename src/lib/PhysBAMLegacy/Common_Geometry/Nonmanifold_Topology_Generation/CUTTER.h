#ifndef __CUTTER_H__
#define __CUTTER_H__

#include <Common_Geometry/Topology/REGULAR_HYPERCUBE_MESH.h> // Input
#include <Nonmanifold_Implicit_Objects/NONMANIFOLD_LEVELSET_MESH.h> // Output
#include <Common_Tools/Arrays/HYBRID_ARRAY.h>
#include <Common_Tools/Arrays/HYBRID_ARRAY_ITERATOR.h>
#include <Common_Geometry/Nonmanifold_Topology_Generation/MATERIAL_PREDICATE.h>
#include <Common_Geometry/Nonmanifold_Topology_Generation/CUTTER_STRATEGY.h>

//#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Data_Structures/UNION_FIND.h>
#include <PhysBAM_Tools/Data_Structures/HASHTABLE.h>
#include <PhysBAM_Tools/Data_Structures/HASHTABLE_ITERATOR.h>
#include <stdint.h>


namespace PhysBAM
{


    template< class T, int d >
    class CUTTER_REGIONS
    {
    public:
        typedef VECTOR<T,d> TV;
        typedef REGULAR_HYPERCUBE_MESH<T,d> T_MESH;
        typedef typename REGULAR_HYPERCUBE_MESH<T,d>::T_CELL T_CELL;
        typedef typename REGULAR_HYPERCUBE_MESH<T,d>::T_FACE T_FACE;

        typedef ARRAY<T_CELL> T_REGION;
        typedef ARRAY<int> T_REGION_IDS;
        typedef ARRAY<T_REGION> T_REGION_ARRAY; 
        typedef ARRAY<T_REGION_IDS> T_REGION_IDS_ARRAY; 

        //typedef VECTOR<int,d> T_INDEX;
        //typedef ARRAY<T_INDEX> T_GRID_REGION;
        //typedef ARRAY<T_GRID_REGION> T_GRID_REGION_ARRAY;
        
        T_REGION_ARRAY regions;
        T_REGION_IDS_ARRAY region_ids;
        //T_GRID_REGION_ARRAY grid_regions;
        int max_index;
        ARRAY<TV> vertices;
    };
    

    template<class T, int d>
    class CUTTER
    {

        typedef VECTOR<T,d> TV;
        typedef VECTOR<int,d> T_INDEX;
        typedef REGULAR_HYPERCUBE_MESH<T,d> T_MESH;
        typedef typename REGULAR_HYPERCUBE_MESH<T,d>::T_CELL T_CELL;
        typedef typename REGULAR_HYPERCUBE_MESH<T,d>::T_FACE T_FACE;

    public:
    
        CUTTER(const REGULAR_HYPERCUBE_MESH<T,d>& mesh);

        void Set_Material_Predicate(  MATERIAL_PREDICATE<T,d>* mp_instance_input) {
            mp_instance = mp_instance_input;
        }

        void Set_Strategy(  CUTTER_STRATEGY<T,d>* strategy_input) {
            strategy = strategy_input;
        }

        void Generate();
        void Clear_Linkages(){
            linkage_list.Clean_Memory();
        }

        void Insert_Linkage(int A, int B, int Cell){
            if( ! linkage_list.Contains( PAIR<VECTOR<int,2>,int>( VECTOR<int,2>(A,B), Cell ) ) )
                linkage_list.Insert( PAIR<VECTOR<int,2>,int>( VECTOR<int,2>(A,B), Cell ));
        }

        const CUTTER_REGIONS<T,d>* GetRegionData() const {
            return &regions;
        }

        const ARRAY<T_CELL>& GetSubCells() const {
            return subcells;
        }

        ARRAY<int> SubCells_From_Root(int root) const {
            return root_cell_to_subcell(root);
        }

        int Root_From_Subcell(int subcell) const {
            return subcell_to_root_cell.Get(subcell);
        }

        const HASHTABLE<int,int>& Subcell_To_Root_Cell() const;
        const ARRAY<ARRAY<int> >& Root_Cell_To_Subcell() const;
          
    private:
        const T_MESH& mesh;

        CUTTER_REGIONS<T,d> regions;

        UNION_FIND<int>  region_sets;
        UNION_FIND<int>  equivalence_sets;
        HASHTABLE<PAIR<VECTOR<int,2>,int> > linkage_list;

        ARRAY<T_CELL> subcells;
        HASHTABLE<int, int> subcell_to_root_cell;
        ARRAY< ARRAY<int> > root_cell_to_subcell;
        ARRAY< ARRAY< int > > subcell_incident_list; // For each subcell, provide a list of incident, touching, cells
        ARRAY< VECTOR<ARRAY< int >, T_MESH::faces_per_cell > > subcell_face_neighbors;
        // For each subcell, provide a list of neighbors along each face

        ARRAY< ARRAY< int > > node_to_subcell; // For each mesh vertex, list cells it belongs to.

        CUTTER_STRATEGY<T,d>* strategy;
        MATERIAL_PREDICATE<T,d>* mp_instance;
    };
}


#endif

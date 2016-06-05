#ifndef __CUTTER_STRATEGY_H__
#define __CUTTER_STRATEGY_H__

#include "MATERIAL_PREDICATE.h"
#include <Common_Geometry/Topology/REGULAR_HYPERCUBE_MESH.h>

namespace PhysBAM{

    template<class T, int d> 
    class CUTTER_STRATEGY {
    public:
        typedef REGULAR_HYPERCUBE_MESH<T,d> T_MESH;
        typedef typename REGULAR_HYPERCUBE_MESH<T,d>::T_CELL T_CELL;
        typedef typename REGULAR_HYPERCUBE_MESH<T,d>::T_FACE T_FACE;

        CUTTER_STRATEGY() : mp_instance(NULL), restartRequired(false) {};
        virtual ~CUTTER_STRATEGY() {};

        void Set_Material_Predicate( MATERIAL_PREDICATE<T,d>* mp_instance_input ) { 
            mp_instance = mp_instance_input;
        }

        void Set_Relationship_Data( const REGULAR_HYPERCUBE_MESH<T,d>& mesh_in,
                                    ARRAY<T_CELL>& subcells_in,
                                    HASHTABLE<PAIR<VECTOR<int,2>,int> >& linkage_list_in,
                                    HASHTABLE<int, int>& subcell_to_root_cell_in,
                                    ARRAY< ARRAY<int> >& root_cell_to_subcell_in,
                                    ARRAY< ARRAY< int > >& subcell_incident_list_in,
                                    ARRAY< VECTOR<ARRAY< int >, T_MESH::faces_per_cell > >& subcell_face_neighbors_in,
                                    ARRAY< ARRAY< int > >& node_to_subcell_in){
            
            mesh = &mesh_in;
            subcells = &subcells_in;
            linkage_list = &linkage_list_in;
            subcell_to_root_cell = &subcell_to_root_cell_in;
            root_cell_to_subcell = &root_cell_to_subcell_in;
            subcell_incident_list = &subcell_incident_list_in;
            subcell_face_neighbors = &subcell_face_neighbors_in;
            node_to_subcell = &node_to_subcell_in;
        }
                                    
                                    

        void ClearRestart() {
            restartRequired = false;
        }

        bool RestartRequired() {
            return restartRequired;
        }
        
        virtual void Reset();
        virtual void Generate_Subcells();
        virtual void Generate_Initial_Equivalence_Classes(  ) {} ;
        virtual void Collapse_Equivalence_Classes(  ) {} ;
        virtual void Renumber_And_Build_Mesh( );
        virtual void Collect_Regions( ) ;
                
    protected:
        MATERIAL_PREDICATE<T,d>* mp_instance;
        bool restartRequired;

        // External Data
        const REGULAR_HYPERCUBE_MESH<T,d>* mesh;
        ARRAY<T_CELL>* subcells;
        HASHTABLE<PAIR<VECTOR<int,2>,int> >* linkage_list;
        HASHTABLE<int, int>* subcell_to_root_cell;
        ARRAY< ARRAY<int> >* root_cell_to_subcell;
        ARRAY< ARRAY< int > >* subcell_incident_list;
        ARRAY< VECTOR<ARRAY< int >, T_MESH::faces_per_cell > >* subcell_face_neighbors;
        ARRAY< ARRAY< int > >* node_to_subcell;

        // Internal Data
        UNION_FIND<int>  equivalence_sets;

    private:
        void BuildIncidentLists();
    };

}

#endif



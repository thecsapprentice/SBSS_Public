#ifndef __CUTTER_NONMANIFOLD_LEVELSET_STRATEGY_H__
#define __CUTTER_NONMANIFOLD_LEVELSET_STRATEGY_H__

#include "CUTTER_STRATEGY.h"
#include <PhysBAM_Tools/Data_Structures/TRIPLE.h>
#include <PhysBAM_Tools/Data_Structures/HASHTABLE.h>
#include <Nonmanifold_Implicit_Objects/NONMANIFOLD_LEVELSET_MESH.h>

namespace PhysBAM{

    template<class T, int d> 
    class CUTTER_NONMANIFOLD_LEVELSET_STRATEGY : public CUTTER_STRATEGY<T,d> {
    public:
        typedef typename CUTTER_STRATEGY<T,d>::T_MESH T_MESH;
        typedef typename CUTTER_STRATEGY<T,d>::T_CELL T_CELL;
        typedef typename CUTTER_STRATEGY<T,d>::T_FACE T_FACE;


        struct TRANSITION_FACE{
            uint32_t type;
            ARRAY<TRIPLE<int,int,int> > left_cells,right_cells;     // entries are ordered in the same order as the bit mask in type (lower bits being lower index)
        };

        CUTTER_NONMANIFOLD_LEVELSET_STRATEGY() {};
        virtual ~CUTTER_NONMANIFOLD_LEVELSET_STRATEGY() {};
        
        virtual void Reset();
        virtual void Generate_Subcells();
        virtual void Generate_Initial_Equivalence_Classes(  );
        virtual void Collapse_Equivalence_Classes(  );       
        virtual void Renumber_And_Build_Mesh( );
        virtual void Collect_Regions( ) ;

        void Initialize_Nonmanifold_Levelset( NONMANIFOLD_LEVELSET_MESH<T,d>& nm_mesh, 
                                              HYBRID_ARRAY<T,d>& phi,
                                              HYBRID_ARRAY<bool,d>& done );
        const HASHTABLE<int,int>& Nmcell_To_Subcell() const;
        const HASHTABLE<int,int>& Subcell_To_Nmcell() const;

    private:
        using CUTTER_STRATEGY<T,d>::mp_instance;
        using CUTTER_STRATEGY<T,d>::restartRequired;

        // External Data
        using CUTTER_STRATEGY<T,d>::mesh;
        using CUTTER_STRATEGY<T,d>::subcells;
        using CUTTER_STRATEGY<T,d>::linkage_list;
        using CUTTER_STRATEGY<T,d>::subcell_to_root_cell;
        using CUTTER_STRATEGY<T,d>::root_cell_to_subcell;
        using CUTTER_STRATEGY<T,d>::subcell_incident_list;
        using CUTTER_STRATEGY<T,d>::subcell_face_neighbors;
        using CUTTER_STRATEGY<T,d>::node_to_subcell;

        // Internal Data
        using CUTTER_STRATEGY<T,d>::equivalence_sets;

        HASHTABLE<int,int> nmcell_to_subcell;
        HASHTABLE<int,int> subcell_to_nmcell;
        ARRAY<TRANSITION_FACE> transition_faces;
        HASHTABLE<uint32_t> allowable_transition_types;
        typedef enum { LEFT=0x1, RIGHT=0x2 } TRANSITION_CORNERS_2D;        
        typedef enum { LOWER_LEFT=0x1, UPPER_LEFT=0x2, LOWER_RIGHT=0x4, UPPER_RIGHT=0x8 } TRANSITION_CORNERS_3D;        

        int BuildTransitionID( int left1, int left2, int left3, int left4,
                               int right1, int right2, int right3, int right4);
        void Set_Allowable_Transition_Types();
        
    };

}

#endif

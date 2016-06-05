#ifndef __CUTTER_BASIC_STRATEGY_H__
#define __CUTTER_BASIC_STRATEGY_H__

#include "CUTTER_STRATEGY.h"
#include <PhysBAM_Tools/Data_Structures/TRIPLE.h>
#include <PhysBAM_Tools/Data_Structures/HASHTABLE.h>

namespace PhysBAM{

    template<class T, int d> 
    class CUTTER_BASIC_STRATEGY : public CUTTER_STRATEGY<T,d> {
    public:
        typedef typename CUTTER_STRATEGY<T,d>::T_MESH T_MESH;
        typedef typename CUTTER_STRATEGY<T,d>::T_CELL T_CELL;
        typedef typename CUTTER_STRATEGY<T,d>::T_FACE T_FACE;

        CUTTER_BASIC_STRATEGY() {};
        virtual ~CUTTER_BASIC_STRATEGY() {};
        
        virtual void Reset();
        virtual void Generate_Subcells();
        virtual void Generate_Initial_Equivalence_Classes(  );
        virtual void Collapse_Equivalence_Classes(  );       
        virtual void Renumber_And_Build_Mesh( );
        virtual void Collect_Regions( ) ;

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
    };

}

#endif

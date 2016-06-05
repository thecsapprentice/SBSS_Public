#include <PhysBAM_Tools/Arrays_Computations/SORT.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Data_Structures/QUEUE.h>
#include <PhysBAM_Tools/Data_Structures/PAIR.h>
#include <PhysBAM_Tools/Data_Structures/TRIPLE.h>
#include <PhysBAM_Tools/Data_Structures/OPERATION_HASH.h>

#include <Common_Geometry/Nonmanifold_Topology_Generation/CUTTER.h>

using namespace PhysBAM;

template<class T,int d> CUTTER<T,d>::
CUTTER(const REGULAR_HYPERCUBE_MESH<T,d>& mesh_input)
    :mesh(mesh_input),mp_instance(NULL)
{
}

template< class T, int d > void
CUTTER<T,d>::Generate() {

    LOG::SCOPE scope("Generating Doman Cuts");

    strategy->Set_Material_Predicate( mp_instance );
    strategy->Set_Relationship_Data( mesh, 
                                     subcells,
                                     linkage_list,
                                     subcell_to_root_cell,
                                     root_cell_to_subcell,
                                     subcell_incident_list,
                                     subcell_face_neighbors,
                                     node_to_subcell );

    strategy->Reset();
    
    while( strategy->RestartRequired() ){
        LOG::SCOPE scope("Restarting...");
        strategy->ClearRestart();
        {
            LOG::SCOPE scope( "Generating Subcells");
            strategy->Generate_Subcells();
        }
        if( strategy->RestartRequired() )
            continue;
        {
            LOG::SCOPE scope( "Determine Initial Equivalence Classes");
            strategy->Generate_Initial_Equivalence_Classes();
        }
        if( strategy->RestartRequired() )
            continue;
        {
            LOG::SCOPE scope( "Collapse Equivalence Classes");
            strategy->Collapse_Equivalence_Classes();
        }
    }
    {
        LOG::SCOPE scope( "Renumber and Build Mesh");
        strategy->Renumber_And_Build_Mesh();
    }
    {
        LOG::SCOPE scope( "Forming Regions" );
        strategy->Collect_Regions();
    }
}

//#####################################################################
// Subcell_To_Root_Cell
//#####################################################################
template<class T,int d> const HASHTABLE<int,int>& CUTTER<T,d>::
Subcell_To_Root_Cell() const
{
    return subcell_to_root_cell;
}
//#####################################################################
// Root_Cell_To_Subcell
//#####################################################################
template<class T,int d> const ARRAY<ARRAY<int> >& CUTTER<T,d>::
Root_Cell_To_Subcell() const
{
    return root_cell_to_subcell;
}


template class CUTTER<float, 2>;
template class CUTTER<float, 3>;
template class CUTTER<double, 2>;
template class CUTTER<double, 3>;



#include "CUTTER_BASIC_STRATEGY.h"
#include <PhysBAM_Tools/Arrays_Computations/SORT.h>

using namespace PhysBAM;

template<class T, int d> void
CUTTER_BASIC_STRATEGY<T,d>::Reset() {
    CUTTER_STRATEGY<T,d>::Reset();
}

template<class T, int d> void
CUTTER_BASIC_STRATEGY<T,d>::Generate_Subcells(){
    CUTTER_STRATEGY<T,d>::Generate_Subcells();
}

template<class T, int d> void
CUTTER_BASIC_STRATEGY<T,d>::Generate_Initial_Equivalence_Classes(){

    int last_progress = 0;
    
    for( int root_cell = 1; root_cell <= (*mesh).elements.m; root_cell++ )
        {
            if( (int)(((T)(root_cell) / (*mesh).elements.m) * 100) % 10 == 0  &&
                (int)(((T)(root_cell) / (*mesh).elements.m) * 100) != last_progress){
                last_progress = (int)(((T)(root_cell) / (*mesh).elements.m) * 100);
                LOG::cout << "Progress: "<< (int)(((T)(root_cell) / (*mesh).elements.m) * 100) << "%" << std::endl;
            }
            
            //const T_INDEX& index = mesh.elements(root_cell).index;
            PHYSBAM_ASSERT( (*mesh).neighbors );
            const VECTOR<PAIR<int,int>, d>& cell_neighbors = (*mesh).neighbors->operator()(root_cell);
            
            for( int i = 1; i <= d; i++)
                {
                    int root_cell_neighbor = cell_neighbors(i).x; // Take the first neighbor of axis
                    
                    if( root_cell_neighbor == -1 ) // No Neighbor on this side;
                        continue;
                    
                    int last_sub_cell = 0;
                    for( int t=1; t <= (*root_cell_to_subcell)(root_cell).m; t++)
                        {
                            last_sub_cell=(*root_cell_to_subcell)(root_cell)(t);
                            int last_sub_cell_neighbor = 0;
                            for( int r=1; r <= (*root_cell_to_subcell)(root_cell_neighbor).m ; r++)
                                {
                                    last_sub_cell_neighbor=(*root_cell_to_subcell)(root_cell_neighbor)(r);
                                    const T_CELL& cellA = (*subcells)( last_sub_cell );
                                    const T_CELL& cellB = (*subcells)( last_sub_cell_neighbor );
                                       
                                    bool isContinous = mp_instance->IsMaterialContinuous(i,
                                                                                         root_cell,
                                                                                         last_sub_cell,
                                                                                         root_cell_neighbor,
                                                                                         last_sub_cell_neighbor);
                                                                        
                                    // This is the trivial condition
                                    if(isContinous){T_FACE faceA,faceB;
                                        faceA=cellA.face(i,0);
                                        faceB=cellB.face(i,1);
                                        for(int j=1;j<=T_MESH::vertices_per_face;j++) equivalence_sets.Union(faceA(j),faceB(j));
                                        //region_sets.Union(last_sub_cell,last_sub_cell_neighbor);
                                    }
                                }
                        }
                }
        }

}

template<class T, int d> void
CUTTER_BASIC_STRATEGY<T,d>::Collapse_Equivalence_Classes(  ){

}

template<class T, int d> void
CUTTER_BASIC_STRATEGY<T,d>::Renumber_And_Build_Mesh( ){

    CUTTER_STRATEGY<T,d>::Renumber_And_Build_Mesh();
}

template<class T, int d> void 
CUTTER_BASIC_STRATEGY<T,d>::Collect_Regions( ){

}


template class CUTTER_BASIC_STRATEGY<float,2>;
template class CUTTER_BASIC_STRATEGY<float,3>;
template class CUTTER_BASIC_STRATEGY<double,2>;
template class CUTTER_BASIC_STRATEGY<double,3>;

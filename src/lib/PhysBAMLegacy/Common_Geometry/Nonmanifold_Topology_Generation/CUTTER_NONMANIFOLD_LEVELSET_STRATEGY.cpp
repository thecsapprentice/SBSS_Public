
#include "CUTTER_NONMANIFOLD_LEVELSET_STRATEGY.h"
#include <PhysBAM_Tools/Arrays_Computations/SORT.h>

using namespace PhysBAM;

template<class T, int d> void
CUTTER_NONMANIFOLD_LEVELSET_STRATEGY<T,d>::Reset() {
    CUTTER_STRATEGY<T,d>::Reset();
    
    // 
    transition_faces.Clean_Memory();
}

template<class T, int d> void
CUTTER_NONMANIFOLD_LEVELSET_STRATEGY<T,d>::Generate_Subcells(){
    CUTTER_STRATEGY<T,d>::Generate_Subcells();
    transition_faces.Clean_Memory();
}


template<class T,int d> int
CUTTER_NONMANIFOLD_LEVELSET_STRATEGY<T,d>::
BuildTransitionID( int left1, int left2, int left3, int left4,
                   int right1, int right2, int right3, int right4)
{

    // Checking order
    PHYSBAM_ASSERT( (left1==0 || left1 > left2) );
    PHYSBAM_ASSERT( ((left1==0 && left2==0) || left2 > left3) );
    PHYSBAM_ASSERT( ((left1==0 && left2==0 && left3==0) || left3 > left4) );
    PHYSBAM_ASSERT( (right1==0 || right1 > right2) );
    PHYSBAM_ASSERT( ((right1==0 && right2==0) || right2 > right3) );
    PHYSBAM_ASSERT( ((right1==0 && right2==0 && right3==0) || right3 > right4) );
 
    // Checking Overlaps
    PHYSBAM_ASSERT( (left1 | left2 | left3 | left4) == 
                    (left1 ^ left2 ^ left3 ^ left4) );

    PHYSBAM_ASSERT( (right1 | right2 | right3 | right4) == 
                    (right1 ^ right2 ^ right3 ^ right4) );
   
    // Checking Completeness
    PHYSBAM_ASSERT( (left1 | left2 | left3 | left4) == 
                    (right1 | right2 | right3 | right4) );


    // Build ID:
    //
    //  31|                                                 |0 
    //    left1 left2 left3 left4 right1 right2 right3 right4
    //

    uint32_t ID = 0;
    ID |= left1 << 28;
    ID |= left2 << 24;
    ID |= left3 << 20;
    ID |= left4 << 16;

    ID |= right1 << 12;
    ID |= right2 << 8;
    ID |= right3 << 4;
    ID |= right4 << 0;
    return ID;
}


template<class T,int d> void
CUTTER_NONMANIFOLD_LEVELSET_STRATEGY<T,d>::
Set_Allowable_Transition_Types()
{

#define INSERT_TRANSITION( L1, L2, L3, L4, R1, R2, R3, R4 ) \
    allowable_transition_types.Insert(BuildTransitionID(L1, L2, L3, L4, R1, R2, R3, R4 ) ); \
    allowable_transition_types.Insert(BuildTransitionID(R1, R2, R3, R4, L1, L2, L3, L4 ) );

    if(d==2){   // 2D
        INSERT_TRANSITION( 0, 0, 0, LEFT|RIGHT, 0, 0, RIGHT, LEFT );
    }
    else{       // 3D

        ////// Two Points Filled
        
        // Single Edge
        INSERT_TRANSITION( 0, 0, 0, LOWER_LEFT|UPPER_LEFT, 0, 0, UPPER_LEFT, LOWER_LEFT );
        INSERT_TRANSITION( 0, 0, 0, LOWER_LEFT|LOWER_RIGHT, 0, 0, LOWER_RIGHT, LOWER_LEFT );
        INSERT_TRANSITION( 0, 0, 0, UPPER_RIGHT|LOWER_RIGHT, 0, 0, UPPER_RIGHT, LOWER_RIGHT );
        INSERT_TRANSITION( 0, 0, 0, UPPER_LEFT|UPPER_RIGHT, 0, 0, UPPER_RIGHT, UPPER_LEFT );
        
        // Opposite Corners
        INSERT_TRANSITION( 0, 0, 0, LOWER_LEFT|UPPER_RIGHT, 0, 0, UPPER_RIGHT, LOWER_LEFT );
        INSERT_TRANSITION( 0, 0, 0, UPPER_LEFT|LOWER_RIGHT, 0, 0, LOWER_RIGHT, UPPER_LEFT );

        ////// Three Points Filled

        // Point And Edge
        INSERT_TRANSITION( 0, 0, 0, LOWER_LEFT|UPPER_LEFT|LOWER_RIGHT, 0, 0, LOWER_RIGHT|LOWER_LEFT, UPPER_LEFT );
        INSERT_TRANSITION( 0, 0, 0, LOWER_LEFT|UPPER_LEFT|LOWER_RIGHT, 0, 0, LOWER_RIGHT, LOWER_LEFT|UPPER_LEFT );

        INSERT_TRANSITION( 0, 0, 0, LOWER_LEFT|UPPER_RIGHT|LOWER_RIGHT, 0, 0, UPPER_RIGHT, LOWER_RIGHT|LOWER_LEFT );
        INSERT_TRANSITION( 0, 0, 0, LOWER_LEFT|UPPER_RIGHT|LOWER_RIGHT, 0, 0, UPPER_RIGHT|LOWER_RIGHT, LOWER_LEFT );

        INSERT_TRANSITION( 0, 0, 0, UPPER_LEFT|UPPER_RIGHT|LOWER_RIGHT, 0, 0, UPPER_RIGHT|UPPER_LEFT, LOWER_RIGHT );
        INSERT_TRANSITION( 0, 0, 0, UPPER_LEFT|UPPER_RIGHT|LOWER_RIGHT, 0, 0, UPPER_RIGHT|LOWER_RIGHT, UPPER_LEFT );

        INSERT_TRANSITION( 0, 0, 0, UPPER_LEFT|UPPER_RIGHT|LOWER_LEFT, 0, 0, UPPER_RIGHT, UPPER_LEFT|LOWER_LEFT );
        INSERT_TRANSITION( 0, 0, 0, UPPER_LEFT|UPPER_RIGHT|LOWER_LEFT, 0, 0, UPPER_RIGHT|UPPER_LEFT, LOWER_LEFT );
        
        // Three Points
        INSERT_TRANSITION( 0, 0, 0, LOWER_LEFT|UPPER_LEFT|LOWER_RIGHT, 0,  LOWER_RIGHT, UPPER_LEFT, LOWER_LEFT );
        INSERT_TRANSITION( 0, 0, 0, LOWER_LEFT|UPPER_RIGHT|LOWER_RIGHT, 0,  UPPER_RIGHT, LOWER_RIGHT, LOWER_LEFT );
        INSERT_TRANSITION( 0, 0, 0, UPPER_LEFT|UPPER_RIGHT|LOWER_RIGHT, 0,  UPPER_RIGHT, LOWER_RIGHT, UPPER_LEFT );
        INSERT_TRANSITION( 0, 0, 0, UPPER_LEFT|UPPER_RIGHT|LOWER_LEFT, 0,  UPPER_RIGHT, UPPER_LEFT, LOWER_LEFT );

        ////// All Points Filled

        // Dual Edges
        INSERT_TRANSITION( 0, 0, 0, LOWER_LEFT|UPPER_LEFT|LOWER_RIGHT|UPPER_RIGHT, 0, 0, UPPER_LEFT|UPPER_RIGHT, LOWER_LEFT|LOWER_RIGHT );
        INSERT_TRANSITION( 0, 0, 0, LOWER_LEFT|UPPER_LEFT|LOWER_RIGHT|UPPER_RIGHT, 0, 0, UPPER_RIGHT|LOWER_RIGHT, UPPER_LEFT|LOWER_LEFT );

        // Point And Two Edges
        INSERT_TRANSITION( 0, 0, 0, LOWER_LEFT|UPPER_LEFT|LOWER_RIGHT|UPPER_RIGHT, 0, 0, UPPER_LEFT|UPPER_RIGHT|LOWER_RIGHT,  LOWER_LEFT );
        INSERT_TRANSITION( 0, 0, 0, LOWER_LEFT|UPPER_LEFT|LOWER_RIGHT|UPPER_RIGHT, 0, 0, LOWER_LEFT|UPPER_RIGHT|LOWER_RIGHT,  UPPER_LEFT );
        INSERT_TRANSITION( 0, 0, 0, LOWER_LEFT|UPPER_LEFT|LOWER_RIGHT|UPPER_RIGHT, 0, 0, UPPER_LEFT|UPPER_RIGHT|LOWER_LEFT,  LOWER_RIGHT );
        INSERT_TRANSITION( 0, 0, 0, LOWER_LEFT|UPPER_LEFT|LOWER_RIGHT|UPPER_RIGHT, 0, 0, UPPER_RIGHT, UPPER_LEFT|LOWER_RIGHT|LOWER_LEFT );

        // Two Points And Edge
        INSERT_TRANSITION( 0, 0, 0, LOWER_LEFT|UPPER_LEFT|LOWER_RIGHT|UPPER_RIGHT, 0, LOWER_RIGHT | UPPER_RIGHT, UPPER_LEFT, LOWER_LEFT );
        INSERT_TRANSITION( 0, 0, 0, LOWER_LEFT|UPPER_LEFT|LOWER_RIGHT|UPPER_RIGHT, 0, UPPER_LEFT | UPPER_RIGHT, LOWER_RIGHT, LOWER_LEFT );
        INSERT_TRANSITION( 0, 0, 0, LOWER_LEFT|UPPER_LEFT|LOWER_RIGHT|UPPER_RIGHT, 0, UPPER_RIGHT, LOWER_RIGHT, UPPER_LEFT | LOWER_LEFT );
        INSERT_TRANSITION( 0, 0, 0, LOWER_LEFT|UPPER_LEFT|LOWER_RIGHT|UPPER_RIGHT, 0, UPPER_RIGHT, LOWER_RIGHT|LOWER_LEFT, UPPER_LEFT );

        // Four Points
        INSERT_TRANSITION( 0, 0, 0, LOWER_LEFT|UPPER_LEFT|LOWER_RIGHT|UPPER_RIGHT, UPPER_RIGHT, LOWER_RIGHT, UPPER_LEFT, LOWER_LEFT );
    }
}


template<class T, int d> void
CUTTER_NONMANIFOLD_LEVELSET_STRATEGY<T,d>::Generate_Initial_Equivalence_Classes(){

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
CUTTER_NONMANIFOLD_LEVELSET_STRATEGY<T,d>::Collapse_Equivalence_Classes(  ){

        {
            LOG::SCOPE scope("STEP 3: Merging pairs of duplicate cells all of whose vertices have been identified");

            for(int root_cell=1;root_cell<=(*mesh).elements.m;root_cell++){
                for(int t1=1;t1<(*root_cell_to_subcell)(root_cell).m;++t1){
                    for(int t2=t1+1;t2<=(*root_cell_to_subcell)(root_cell).m;++t2){
                        int subcell1=(*root_cell_to_subcell)(root_cell)(t1),subcell2=(*root_cell_to_subcell)(root_cell)(t2);
                        const T_CELL& cellA=(*subcells)(subcell1);const T_CELL& cellB=(*subcells)(subcell2);

                        bool all_vertices_identified=true;
                        for(int j=1;j<=T_MESH::vertices_per_cell;++j)
                            if(equivalence_sets.Find(cellA.vertices(j))!=equivalence_sets.Find(cellB.vertices(j))){
                                all_vertices_identified=false;break;}

                        if(all_vertices_identified){
                            int cell1_fragment_representative=mp_instance->Subcell_Material_Representative(subcell1);
                            PHYSBAM_ASSERT(cell1_fragment_representative != -1);

                            int cell2_fragment_representative=mp_instance->Subcell_Material_Representative(subcell2);
                            PHYSBAM_ASSERT(cell2_fragment_representative != -1);

                            if( cell1_fragment_representative < cell2_fragment_representative )
                                (*linkage_list).Insert(PAIR<VECTOR<int,2>,int>(VECTOR<int,2>(cell1_fragment_representative,cell2_fragment_representative),root_cell));
                            else
                                (*linkage_list).Insert(PAIR<VECTOR<int,2>,int>(VECTOR<int,2>(cell2_fragment_representative,cell1_fragment_representative),root_cell));
                            
                            LOG::cout<<"Appending linkage list with simplices "<<min(cell1_fragment_representative,cell2_fragment_representative)<<" and "<<max(cell1_fragment_representative,cell2_fragment_representative)<<" for cell "<<root_cell<<" with index "<<(*mesh).elements(root_cell).index<<std::endl;
                            restartRequired=true;}}}}

            if(restartRequired) return;
        }


        {
            LOG::SCOPE("STEP 4: Merging pairs of duplicate cells which have an identified vertex with negative sign");

            for(int root_cell=1;root_cell<=(*mesh).elements.m;++root_cell){
                for(int t1=1;t1<(*root_cell_to_subcell)(root_cell).m;++t1){
                    for(int t2=t1+1;t2<=(*root_cell_to_subcell)(root_cell).m;++t2){
                        int subcell1=(*root_cell_to_subcell)(root_cell)(t1),subcell2=(*root_cell_to_subcell)(root_cell)(t2);
                        const T_CELL& cellA=(*subcells)(subcell1);const T_CELL& cellB=(*subcells)(subcell2);

                        for(int j=1;j<=T_MESH::vertices_per_cell;++j){
                            if(equivalence_sets.Find(cellA.vertices(j))==equivalence_sets.Find(cellB.vertices(j))){
                                bool vertex_subcell1_inside_material=mp_instance->InsideMaterial(subcell1,j);
                                bool vertex_subcell2_inside_material=mp_instance->InsideMaterial(subcell2,j);
                                if(vertex_subcell1_inside_material && vertex_subcell2_inside_material){
                                    int cell1_fragment_representative=mp_instance->Node_Material_Representative(subcell1,j);
                                    int cell2_fragment_representative=mp_instance->Node_Material_Representative(subcell2,j);
                                    PHYSBAM_ASSERT(!(cell1_fragment_representative == cell2_fragment_representative));

                                    if( cell1_fragment_representative < cell2_fragment_representative )
                                        (*linkage_list).Insert(PAIR<VECTOR<int,2>,int>(VECTOR<int,2>(cell1_fragment_representative,cell2_fragment_representative),root_cell));
                                    else
                                        (*linkage_list).Insert(PAIR<VECTOR<int,2>,int>(VECTOR<int,2>(cell2_fragment_representative,cell1_fragment_representative),root_cell));

                                    LOG::cout<<"Appending linkage list with simplices "<<min(cell1_fragment_representative,cell2_fragment_representative)<<" and "<<max(cell1_fragment_representative,cell2_fragment_representative)<<" for cell "<<root_cell<<" with index "<<(*mesh).elements(root_cell).index<<std::endl;
                                    restartRequired=true;
                                    
                                    break; // Don't need to continue looking once the cell is merged
                                }}}}}}

            if(restartRequired) return;
        }


        {
            LOG::SCOPE("STEP 5: Generating transition faces or merging duplicate cells to avoid exotic transitions");


            // store all translated faces in a hashtable
            HASHTABLE<VECTOR<int,T_MESH::vertices_per_face>,ARRAY<TRIPLE<int,int,int> > > identified_faces;
            for(int root_cell=1;root_cell<=(*mesh).elements.m;++root_cell){
                for(int t=1;t<=(*root_cell_to_subcell)(root_cell).m;++t){
                    int subcell=(*root_cell_to_subcell)(root_cell)(t);
                    for(int axis=1;axis<=d;++axis) for(int side=0;side<=1;++side){
                        const VECTOR<int,T_MESH::vertices_per_face> face_indices=(*subcells)(subcell).face(axis,side);
                        VECTOR<int,T_MESH::vertices_per_face> translated_face;
                        for(int k=1;k<=T_MESH::vertices_per_face;++k) translated_face(k)=equivalence_sets.Find(face_indices(k));
                        ARRAY<TRIPLE<int,int,int> > mapped_cells=identified_faces.Get_Or_Insert(translated_face,ARRAY<TRIPLE<int,int,int> >());
                        mapped_cells.Append(TRIPLE<int,int,int>(root_cell,t,2*(axis-1)+side));
                        identified_faces.Set(translated_face,mapped_cells);}}}

            for(HASHTABLE_ITERATOR<VECTOR<int,T_MESH::vertices_per_face>,ARRAY<TRIPLE<int,int,int> > > iterator(identified_faces);iterator.Valid();iterator.Next()){
                const ARRAY<TRIPLE<int,int,int> >& all_cells=iterator.Data();

                UNION_FIND<int> partitions;partitions.Initialize(all_cells.m);
                for(int i=1;i<all_cells.m;++i) for(int j=i+1;j<=all_cells.m;++j){
                    int subcell1=(*root_cell_to_subcell)(all_cells(i).x)(all_cells(i).y),subcell2=(*root_cell_to_subcell)(all_cells(j).x)(all_cells(j).y);
                    int axis=1;
                    if(all_cells(i).z>=2 && all_cells(i).z<=3) axis=2;
                    else if(all_cells(i).z>=4 && all_cells(i).z<=5) axis=3;

                    const VECTOR<int,T_MESH::vertices_per_face> subcell1_face_indices=(*subcells)(subcell1).face_indices(axis,all_cells(i).z%2);
                    const VECTOR<int,T_MESH::vertices_per_face> subcell2_face_indices=(*subcells)(subcell2).face_indices(axis,all_cells(j).z%2);

                    for(int t=1;t<=T_MESH::vertices_per_face;++t){
                        bool vertex_subcell1_inside_material=mp_instance->InsideMaterial(subcell1,subcell1_face_indices(t));
                        bool vertex_subcell2_inside_material=mp_instance->InsideMaterial(subcell2,subcell2_face_indices(t));
                        if(vertex_subcell1_inside_material && vertex_subcell2_inside_material && partitions.Find(i)!=partitions.Find(j)) partitions.Union(i,j);}}

                // collect partitions
                HASHTABLE<int,ARRAY<TRIPLE<int,int,int> > > left_subcells,right_subcells;;
                for(int i=1;i<=all_cells.m;i++){
                    ARRAY<TRIPLE<int,int,int> > lc=left_subcells.Get_Or_Insert(partitions.Find(i),ARRAY<TRIPLE<int,int,int> >()),rc=right_subcells.Get_Or_Insert(partitions.Find(i),ARRAY<TRIPLE<int,int,int> >());
                    if(all_cells(i).z==0 || all_cells(i).z==2 || all_cells(i).z==4) rc.Append(TRIPLE<int,int,int>(all_cells(i).x,all_cells(i).y,all_cells(i).z));
                    else lc.Append(TRIPLE<int,int,int>(all_cells(i).x,all_cells(i).y,all_cells(i).z));
                    left_subcells.Set(partitions.Find(i),lc);right_subcells.Set(partitions.Find(i),rc);}

                for(HASHTABLE_ITERATOR<int,ARRAY<TRIPLE<int,int,int> > > partition_iterator(left_subcells);partition_iterator.Valid();partition_iterator.Next()){
                    const ARRAY<TRIPLE<int,int,int> >& lefties=partition_iterator.Data();const ARRAY<TRIPLE<int,int,int> > righties=right_subcells.Get(partition_iterator.Key());
                    bool disallow_transition=false;

                    if(lefties.m>1 || righties.m>1){
                        if(lefties.m>=4 || righties.m>=4) disallow_transition=true;
                        uint32_t type=0;
                        ARRAY<uint32_t> left_types,right_types;
                        ARRAY<int> left_indices,right_indices;
                        // assemble left types
                        for(int i=1;i<=lefties.m;++i){uint32_t current_type=0;
                            int subcell=(*root_cell_to_subcell)(lefties(i).x)(lefties(i).y),axis=1;
                            if(lefties(i).z>=2 && lefties(i).z<=3) axis=2;
                            else if(lefties(i).z>=4 && lefties(i).z<=5) axis=3;
                            const VECTOR<int,T_MESH::vertices_per_face> subcell_face_indices=(*subcells)(subcell).face_indices(axis,lefties(i).z%2);

                            for(int t=1;t<=T_MESH::vertices_per_face;++t) if(mp_instance->InsideMaterial(subcell,subcell_face_indices(t)))
                                current_type=current_type|(1<<(t-1));
                            left_types.Append(current_type);left_indices.Append(i);}

                        // assemble right types
                        for(int i=1;i<=righties.m;++i){uint32_t current_type=0;
                            int subcell=(*root_cell_to_subcell)(righties(i).x)(righties(i).y),axis=1;
                            if(righties(i).z>=2 && righties(i).z<=3) axis=2;
                            else if(righties(i).z>=4 && righties(i).z<=5) axis=3;
                            const VECTOR<int,T_MESH::vertices_per_face> subcell_face_indices=(*subcells)(subcell).face_indices(axis,righties(i).z%2);

                            for(int t=1;t<=T_MESH::vertices_per_face;++t) if(mp_instance->InsideMaterial(subcell,subcell_face_indices(t)))
                                current_type=current_type|(1<<(t-1));
                            right_types.Append(current_type);right_indices.Append(i);}

                        // sort left and right indices using indirect comparison with types
                        Sort(left_indices,Indirect_Comparison(left_types));Sort(right_indices,Indirect_Comparison(right_types));
                        // sort left and right types
                        Sort(left_types);Sort(right_types);     // smallest comes first and so on...
                        for(int i=left_types.m+1;i<=4;++i) left_types.Append((uint32_t)0);
                        for(int i=right_types.m+1;i<=4;++i) right_types.Append((uint32_t)0);

                        LOG::cout<<"Left Types: "<<left_types<<std::endl;
                        LOG::cout<<"Right Types: "<<right_types<<std::endl;

                        int shift=0;
                        for(int i=1;i<=right_types.m;++i){type=type|(right_types(i)<<shift);shift+=4;}
                        for(int i=1;i<=left_types.m;++i){type=type|(left_types(i)<<shift);shift+=4;}
                        LOG::cout<<"Type: "<<type<<std::endl;

                        if(!allowable_transition_types.Contains(type)) disallow_transition=true;
                        else{LOG::cout<<"Allowable transition face detected!!!"<<std::endl;
                            TRANSITION_FACE transition;transition.type=type;
                            LOG::cout<<"Left cells:";
                            for(int i=1;i<=lefties.m;++i){LOG::cout<<" "<<(*subcells)((*root_cell_to_subcell)(lefties(i).x)(lefties(i).y)).index;
                                transition.left_cells.Append(lefties(left_indices(i)));}
                            LOG::cout<<", Right cells:";
                            for(int i=1;i<=righties.m;++i){LOG::cout<<" "<<(*subcells)((*root_cell_to_subcell)(righties(i).x)(righties(i).y)).index;
                                transition.right_cells.Append(righties(right_indices(i)));}
                            LOG::cout<<std::endl;
                            transition_faces.Append(transition);}}

                    if(disallow_transition){
                        // collapse all left cells
                        for(int i=1;i<lefties.m;++i) for(int j=i+1;j<=lefties.m;++j){
                            int subcell1=(*root_cell_to_subcell)(lefties(i).x)(lefties(i).y),subcell2=(*root_cell_to_subcell)(lefties(j).x)(lefties(j).y);

                            int cell1_fragment_representative=mp_instance->Subcell_Material_Representative(subcell1);
                            PHYSBAM_ASSERT(cell1_fragment_representative != -1);
                            int cell2_fragment_representative=mp_instance->Subcell_Material_Representative(subcell2);
                            PHYSBAM_ASSERT(cell2_fragment_representative != -1);

                            if( cell1_fragment_representative < cell2_fragment_representative )
                                (*linkage_list).Insert(PAIR<VECTOR<int,2>,int>(VECTOR<int,2>(cell1_fragment_representative,cell2_fragment_representative),lefties(i).x));
                            else
                                (*linkage_list).Insert(PAIR<VECTOR<int,2>,int>(VECTOR<int,2>(cell2_fragment_representative,cell1_fragment_representative),lefties(i).x));

                            LOG::cout<<"Appending linkage list with simplices "<<min(cell1_fragment_representative,cell2_fragment_representative)<<" and "<<max(cell1_fragment_representative,cell2_fragment_representative)<<" for cell "<<lefties(i).x<<" with index "<<(*mesh).elements(lefties(i).x).index<<std::endl;
                            restartRequired=true;}

                        // collapse all right cells
                        for(int i=1;i<righties.m;++i) for(int j=i+1;j<=righties.m;++j){
                            int subcell1=(*root_cell_to_subcell)(righties(i).x)(righties(i).y),subcell2=(*root_cell_to_subcell)(righties(j).x)(righties(j).y);

                            int cell1_fragment_representative=mp_instance->Subcell_Material_Representative(subcell1);
                            PHYSBAM_ASSERT(cell1_fragment_representative != -1);
                            int cell2_fragment_representative=mp_instance->Subcell_Material_Representative(subcell2);
                            PHYSBAM_ASSERT(cell2_fragment_representative != -1);

                            if( cell1_fragment_representative < cell2_fragment_representative )
                                (*linkage_list).Insert(PAIR<VECTOR<int,2>,int>(VECTOR<int,2>(cell1_fragment_representative,cell2_fragment_representative),righties(i).x));
                            else
                                (*linkage_list).Insert(PAIR<VECTOR<int,2>,int>(VECTOR<int,2>(cell2_fragment_representative,cell1_fragment_representative),righties(i).x));

                            LOG::cout<<"Appending linkage list with simplices "<<min(cell1_fragment_representative,cell2_fragment_representative)<<" and "<<max(cell1_fragment_representative,cell2_fragment_representative)<<" for cell "<<righties(i).x<<" with index "<<(*mesh).elements(righties(i).x).index<<std::endl;

                            restartRequired=true;}}}}

            if(restartRequired) return;
        }
}

template<class T, int d> void
CUTTER_NONMANIFOLD_LEVELSET_STRATEGY<T,d>::Renumber_And_Build_Mesh( ){

    {
        LOG::SCOPE("STEP 6: Recompute equivalence_sets and region_sets identifying all positive nodes and all negative nodes");

        UNION_FIND<int> old_equivalence_sets=equivalence_sets;
        //UNION_FIND<int> old_region_sets=region_sets;

        equivalence_sets.Clear_Connectivity();
        equivalence_sets.Initialize((*subcells).m*T_MESH::vertices_per_cell);
        //region_sets.Clear_Connectivity();
        //region_sets.Initialize((*subcells).m);

        for(int root_cell=1;root_cell<=(*mesh).elements.m;root_cell++){
            const VECTOR<PAIR<int,int>,d>& cell_neighbors=(*mesh).neighbors->operator()(root_cell);
            for(int i=1;i<=d;i++){int root_cell_neighbor=cell_neighbors(i).x;   // Take the first neighbor of axis
                if(root_cell_neighbor==-1) continue;                            // No Neighbor on this side;
                for(int t=1;t<=(*root_cell_to_subcell)(root_cell).m;t++){
                    int last_sub_cell=(*root_cell_to_subcell)(root_cell)(t);
                    for(int r=1;r<=(*root_cell_to_subcell)(root_cell_neighbor).m;r++){
                        int last_sub_cell_neighbor=(*root_cell_to_subcell)(root_cell_neighbor)(r);
                        const T_CELL& cellA=(*subcells)(last_sub_cell);
                        const T_CELL& cellB=(*subcells)(last_sub_cell_neighbor);
                        T_FACE faceA=cellA.face(i,0),faceB=cellB.face(i,1);
                        const VECTOR<int,T_MESH::vertices_per_face> face_indicesA=cellA.face_indices(i,0),face_indicesB=cellB.face_indices(i,1);

                        //bool merge_regions=(old_region_sets.Find(last_sub_cell)==old_region_sets.Find(last_sub_cell_neighbor));
                        for(int j=1;j<=T_MESH::vertices_per_face;j++){
                            if(old_equivalence_sets.Find(faceA(j))==old_equivalence_sets.Find(faceB(j))){
                                bool vertex1_inside_material=mp_instance->InsideMaterial(last_sub_cell,face_indicesA(j));
                                bool vertex2_inside_material=mp_instance->InsideMaterial(last_sub_cell_neighbor,face_indicesB(j));
                                if(vertex1_inside_material==vertex2_inside_material)
                                    equivalence_sets.Union(faceA(j),faceB(j));
                                //else 
                                //    merge_regions=false;
                            }
                            //else
                            //    merge_regions=false;
                        }
                        //if(merge_regions) region_sets.Union(last_sub_cell,last_sub_cell_neighbor);
                    }
                }
            }
        }
    }   

    CUTTER_STRATEGY<T,d>::Renumber_And_Build_Mesh();
}

template<class T, int d> void 
CUTTER_NONMANIFOLD_LEVELSET_STRATEGY<T,d>::Collect_Regions( ){

}

//#####################################################################
// Nmcell_To_Subcell
//#####################################################################
template<class T,int d> const HASHTABLE<int,int>&
CUTTER_NONMANIFOLD_LEVELSET_STRATEGY<T,d>::
Nmcell_To_Subcell() const
{
    return nmcell_to_subcell;
}
//#####################################################################
// Subcell_To_Nmcell
//#####################################################################
template<class T,int d> const HASHTABLE<int,int>&
CUTTER_NONMANIFOLD_LEVELSET_STRATEGY<T,d>::
Subcell_To_Nmcell() const
{
    return subcell_to_nmcell;
}
//#####################################################################

template<class T, int d> void
CUTTER_NONMANIFOLD_LEVELSET_STRATEGY<T,d>::
Initialize_Nonmanifold_Levelset( NONMANIFOLD_LEVELSET_MESH<T,d>& nm_mesh, 
                                 HYBRID_ARRAY<T,d>& phi,
                                 HYBRID_ARRAY<bool,d>& done )
{
    LOG::SCOPE scope("Build Nonmanifold_Levelset_Mesh From Cutter...");
    typedef VECTOR<int,d> T_INDEX;
    typedef VECTOR<T,d> TV;

    //###################################################################
    //          BUILD NONMANIFOLD_LEVELSET_MESH FROM CUTTER
    //###################################################################
    
    nm_mesh.cells.Remove_All();
    nm_mesh.cell_is_grid.Resize( 0, true, false);
    nm_mesh.node_domain = RANGE<T_INDEX>();
    nm_mesh.node_mesh_count = 0;

    typedef PAIR<T_INDEX,int> HINDEX;
    ARRAY<RANGE<TV> > leaf_boxes;
    HASHTABLE<T_INDEX,ARRAY<int> > node_instances;
    HASHTABLE<int,HINDEX> node_map;
    nmcell_to_subcell.Clean_Memory();
    subcell_to_nmcell.Clean_Memory();

    //const CUTTER_REGIONS<T,d>* regions = GetRegionData();
    for( int subcell = 1; subcell <= (*subcells).m; subcell++){
        //for( int r = 1; r <= regions->regions.m; r++){
        //const typename CUTTER_REGIONS<T,d>::T_REGION& region = regions->regions(r);
        //for( int c = 1; c <= region.m; c++ ){
        const T_CELL& cell = (*subcells)(subcell);
        const T_INDEX& cell_index = cell.index;
        typename NONMANIFOLD_LEVELSET_MESH<T,d>::T_CELL mesh_cell;
        mesh_cell.cell_index = cell_index;
        int root_cell = (*subcell_to_root_cell).Get(subcell);
        int index=1;
        bool is_grid = true;
        RANGE<TV> cell_box;
        for(RANGE_ITERATOR<d> node_iterator(RANGE<T_INDEX>(cell_index,cell_index+1));node_iterator.Valid();node_iterator.Next()){
            if( nm_mesh.node_domain.Empty() )
                nm_mesh.node_domain = RANGE<T_INDEX>(node_iterator.Index(), node_iterator.Index());
            else
                nm_mesh.node_domain.Enlarge_To_Include_Point( node_iterator.Index() );

            if( cell_box.Empty() )
                cell_box = RANGE<TV>( (*mesh).Node(root_cell,index), (*mesh).Node(root_cell,index) );
            else
                cell_box.Enlarge_To_Include_Point( (*mesh).Node(root_cell,index) );
            
            ARRAY<int> ni = node_instances.Get_Default(node_iterator.Index(), ARRAY<int>() );
            if( !ni.Find(cell.vertices(index)) ){
                HINDEX nm_index;
                if( ni.m >= 1 ){
                    //PHYSBAM_FATAL_ERROR();
                    nm_index.x = T_INDEX();
                    nm_index.y = nm_mesh.node_mesh_count+1;
                    nm_mesh.node_mesh_count+=1;
                }
                else
                    {
                        nm_index.x = node_iterator.Index();
                        nm_index.y = 0;
                    }
                //LOG::cout << "Inserting vertex " << cell.vertices(index) << " with hindex " << nm_index << std::endl;
                node_map.Insert(cell.vertices(index), nm_index); 
                ni.Append( cell.vertices(index) );
            }
            node_instances.Set( node_iterator.Index(), ni );
            mesh_cell.nodes(index) = node_map.Get( cell.vertices(index) );
            nm_mesh.node_locations.Set( mesh_cell.nodes(index), (*mesh).Node( root_cell, index ));
            is_grid = is_grid && mesh_cell.nodes(index).y == 0;
            index++;
        }
        int nm_cell_id = nm_mesh.cells.Append( mesh_cell );
        nmcell_to_subcell.Set( nm_cell_id, subcell );
        subcell_to_nmcell.Set( subcell, nm_cell_id );
        nm_mesh.cell_is_grid.Append( is_grid );
        leaf_boxes.Append( cell_box );
    }
    
    nm_mesh.aabb_hierarchy.Set_Leaf_Boxes(leaf_boxes, true);
    nm_mesh.Update_Node_Neighbors();
    nm_mesh.initialized = true;

    // initialize transition faces
    for(int i=1;i<=transition_faces.m;++i){
        typename NONMANIFOLD_LEVELSET_MESH<T,d>::TRANSITION_FACE transition;
        transition.type=transition_faces(i).type;
        for(int j=1;j<=transition_faces(i).left_cells.m;++j){
            int subcell=(*root_cell_to_subcell)(transition_faces(i).left_cells(j).x)(transition_faces(i).left_cells(j).y);
            transition.left_cells.Append(PAIR<int,int>(subcell_to_nmcell.Get(subcell),transition_faces(i).left_cells(j).z+1));}
        for(int j=1;j<=transition_faces(i).right_cells.m;++j){
            int subcell=(*root_cell_to_subcell)(transition_faces(i).right_cells(j).x)(transition_faces(i).right_cells(j).y);
            transition.right_cells.Append(PAIR<int,int>(subcell_to_nmcell.Get(subcell),transition_faces(i).right_cells(j).z+1));}
        nm_mesh.transition_faces.Append(transition);}

    LOG::cout<<nm_mesh.transition_faces.m<<" transition faces initialized!"<<std::endl;

    nm_mesh.Compute_Acceleration_Structures_For_Transition_Faces();

    // compute cell neighbors
    for(int nm_cell=1;nm_cell<=nm_mesh.cells.m;++nm_cell){int subcell=nmcell_to_subcell.Get(nm_cell);
        for(int i=1;i<=T_MESH::faces_per_cell;++i){
            //if(nm_mesh.cells(nm_cell).transition(i).x!=0)
            //    nm_mesh.cells(nm_cell).cell_neighbors(i)=PAIR<bool,int>(true,nm_mesh.cells(nm_cell).transition(i).x);
            //else{
            const ARRAY<int>& face_neighbors=(*subcell_face_neighbors)(subcell)(i);
            //PHYSBAM_ASSERT(face_neighbors.m==1 || face_neighbors.m==0);
/*
            if( !(face_neighbors.m==1 || face_neighbors.m==0) ){
                LOG::cout << "FATAL ERROR: More than one face neighbor detected for a non-transition face!!!!" << std::endl;
                LOG::cout << "  Error Detected for Subcell Cell " << subcell << ", Face " << i << std::endl;
                const T_CELL& cell = (*subcells)(subcell);
                LOG::cout << "  Index: " << cell.index << std::endl;
                LOG::cout << "  Root: " << (*subcell_to_root_cell).Get(subcell) << std::endl;
                LOG::cout << "  Neighbors: " << std::endl;
                for( int n = 1; n <= face_neighbors.m; n++)
                    LOG::cout << "    " << face_neighbors(n) << std::endl;
            }
*/

            //if(face_neighbors.m==1) nm_mesh.cells(nm_cell).cell_neighbors(i)=PAIR<bool,int>(false,subcell_to_nmcell.Get(face_neighbors(1)));
            //else nm_mesh.cells(nm_cell).cell_neighbors(i)=PAIR<bool,int>(false,0);
            nm_mesh.cells(nm_cell).cell_neighbors(i) = face_neighbors;
        }}

    phi = HYBRID_ARRAY<T,d>( nm_mesh.Node_Domain(), nm_mesh.Node_Mesh_Counts() );
    phi.Fill(FLT_MAX);
    done = HYBRID_ARRAY<bool,d>(nm_mesh.Node_Domain(), nm_mesh.Node_Mesh_Counts());
    done.Fill(false);
    
    {
        LOG::SCOPE scope("Computing Nodal Distances near interface...");
        for( int nmc = 1; nmc <= nm_mesh.cells.m; nmc++){
            int subcell_id = nmcell_to_subcell.Get( nmc );
            int root_id = (*subcell_to_root_cell).Get(subcell_id);
            //LOG::cout << "Computing distances for cell " << nm_mesh.cells(nmc).cell_index << std::endl;
            ARRAY<TV> nodal_positions;
            ARRAY<T> distances;
            //LOG::cout << "Pre: " << std::endl;
            for( int i = 1; i <= T_MESH::vertices_per_cell; i++){            
                nodal_positions.Append( (*mesh).Node( root_id, i ) );
                distances.Append( phi(nm_mesh.cells(nmc).nodes(i)) );            
                //LOG::cout <<"\t"<< nodal_positions.Last() << "   " << distances.Last() << std::endl;
                
            }
            mp_instance->ComputeNodalDistances( nodal_positions, subcell_id, distances );
            //LOG::cout << "Post: " << std::endl;
            for( int i = 1; i <= T_MESH::vertices_per_cell; i++){
                //LOG::cout << "\t" << distances(i) << std::endl;
                phi(nm_mesh.cells(nmc).nodes(i)) = distances(i);
                if( abs(phi(nm_mesh.cells(nmc).nodes(i))) <= 2*(*mesh).dx ){
                    done( nm_mesh.cells(nmc).nodes(i)) = true;                
                }          
                else{
                    //intt sign = phi(nm_mesh.cells(nmc).nodes(i)) > 0 ? 1 : -1;
                    //if( phi(nm_mesh.cells(nmc).nodes(i)) != 0.0f )
                    //    phi(nm_mesh.cells(nmc).nodes(i)) = FLT_MAX * sign;
                }
            }       
        }
    }

    {
        LOG::SCOPE scope("Running Consistency Checks...");
        for( int nmc = 1; nmc <= nm_mesh.cells.m; nmc++){
            int scell = nmcell_to_subcell.Get( nmc );
            int root_id = (*subcell_to_root_cell).Get(scell);
            T_CELL& subcell = (*subcells)(scell);
            //LOG::cout << "Checking subcell " << scell << std::endl;
            //LOG::cout << subcell.index << std::endl;
            for( int f = 1; f <= T_MESH::faces_per_cell; f++){
                int axis = ((f-1) / 2)+1;
                bool direction = bool( ((f-1) % 2) );
                VECTOR<int,T_MESH::vertices_per_face> face_indices = subcell.face_indices( axis, direction );
                
                // If we have neighbors in this direction, nevermind
                if( (*subcell_face_neighbors)(scell)(f).m != 0 ) {
                    //LOG::cout << "We have neighbors along face " << f << std::endl;
                    continue;
                }
                
                //LOG::cout << "We DO NOT have neighbors along face " << f << std::endl;
                
                ARRAY<TV> nodal_positions;
                ARRAY<T> distances;
                for( int i = 1; i <= T_MESH::vertices_per_cell; i++){            
                    nodal_positions.Append( (*mesh).Node( root_id, i ) );
                    distances.Append( phi(nm_mesh.cells(nmc).nodes(i)) );
                }

                int positive=0, negative=0;
                for( int v = 1; v <= T_MESH::vertices_per_face; v++){
                    //LOG::cout << "Sign at position "<< face_indices(v) << "  " << distances(face_indices(v)) << std::endl;
                    if( sign(distances(face_indices(v))) < 0 ) negative ++;
                    else positive ++;
                }
                if( positive && negative ){
                    LOG::cout<<"For cell index: "<<nm_mesh.cells(nmc).cell_index<<" ";
                    LOG::cout << "Interface collides with mesh boundary in cell " << nmc << " (root " << root_id << ", index " << subcell.index << " ) along face " << f << std::endl;}
                
            }
        }
    }

}


template class CUTTER_NONMANIFOLD_LEVELSET_STRATEGY<float,2>;
template class CUTTER_NONMANIFOLD_LEVELSET_STRATEGY<float,3>;
template class CUTTER_NONMANIFOLD_LEVELSET_STRATEGY<double,2>;
template class CUTTER_NONMANIFOLD_LEVELSET_STRATEGY<double,3>;

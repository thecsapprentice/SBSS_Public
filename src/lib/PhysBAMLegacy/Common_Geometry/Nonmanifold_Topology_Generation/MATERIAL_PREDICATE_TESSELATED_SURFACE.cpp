#include <Common_Geometry/Nonmanifold_Topology_Generation/MATERIAL_PREDICATE_TESSELATED_SURFACE.h>
#include <PhysBAM_Tools/Data_Structures/UNION_FIND.h>

using namespace PhysBAM;

template<class T,int d>
MATERIAL_PREDICATE_TESSELATED_SURFACE<T,d>::MATERIAL_PREDICATE_TESSELATED_SURFACE(const T_VOLUME& current_tetrahedralized_volume_input,const ARRAY<VECTOR<bool,4> >& material_nodes_per_duplicate_tet_input,const ARRAY<int>& old_particle_per_new_collapsed_particle_input,const HASHTABLE<VECTOR<int,2>,ARRAY<VECTOR<int,2> > >& template_hex_node_to_template_tet_nodes_input,const ARRAY<int>& duplicate_tet_to_parent_tet_input,const ARRAY<ARRAY<int> >& template_hex_to_template_tets_input,const ARRAY<ARRAY<int> >& template_tet_to_duplicate_tets_input,const ARRAY<T>& duplicate_tet_nodal_distances_input)
    :current_tetrahedralized_volume(current_tetrahedralized_volume_input),material_nodes_per_duplicate_tet(material_nodes_per_duplicate_tet_input),old_particle_per_new_collapsed_particle(old_particle_per_new_collapsed_particle_input),template_hex_node_to_template_tet_nodes(template_hex_node_to_template_tet_nodes_input),duplicate_tet_to_parent_tet(duplicate_tet_to_parent_tet_input),template_hex_to_template_tets(template_hex_to_template_tets_input),template_tet_to_duplicate_tets(template_tet_to_duplicate_tets_input),duplicate_tet_nodal_distances(duplicate_tet_nodal_distances_input)
{
}

template<class T,int d>
MATERIAL_PREDICATE_TESSELATED_SURFACE<T,d>::~MATERIAL_PREDICATE_TESSELATED_SURFACE()
{
}

template<class T,int d> void
MATERIAL_PREDICATE_TESSELATED_SURFACE<T,d>::MaterialFragments(const T_MESH& mesh,const HASHTABLE<PAIR<VECTOR<int,2>,int> >& linkage_list,ARRAY<T_CELL>& sub_cells,HASHTABLE<int,int>& subcell_to_root_cell)
{
    int vertex_counter=1;
    sub_cells.Clean_Memory();
    subcell_to_root_cell.Clean_Memory();
    material_fragments.Clean_Memory();
    duplicate_cell_to_root_cell.Clean_Memory();

    int last_progress=0;
    for(int template_hex_index=1;template_hex_index<=mesh.elements.m;template_hex_index++){
        ARRAY<int> duplicate_tets; // all duplicate tets of this template hex
        const ARRAY<int>& template_tets=template_hex_to_template_tets(template_hex_index);
        for(int i=1;i<=template_tets.m;i++){
            const ARRAY<int>& duplicate_tets_of_template_tet=template_tet_to_duplicate_tets(template_tets(i));
            duplicate_tets.Append_Elements(duplicate_tets_of_template_tet);}

        UNION_FIND<int> duplicate_tet_partitions; // partitioning all tets into connected fragments
        duplicate_tet_partitions.Initialize(duplicate_tets.m);
        for(int i=1;i<duplicate_tets.m;i++) for(int j=i+1;j<=duplicate_tets.m;j++){
            const VECTOR<int,2> indices=VECTOR<int,2>(duplicate_tets(i),duplicate_tets(j)).Sorted();
            if(current_tetrahedralized_volume.mesh.Face_Neighbors(duplicate_tets(i),duplicate_tets(j)) || linkage_list.Contains(PAIR<VECTOR<int,2>,int>(indices,template_hex_index)))
               duplicate_tet_partitions.Union(i,j);}
        HASHTABLE<int,ARRAY<int> > components;
        for(int i=1;i<=duplicate_tets.m;i++){
            ARRAY<int> c=components.Get_Or_Insert(duplicate_tet_partitions.Find(i),ARRAY<int>());
            c.Append(duplicate_tets(i));
            components.Set(duplicate_tet_partitions.Find(i),c);}

        int num_components=0;
        for(HASHTABLE_ITERATOR<int,ARRAY<int> > iterator(components);iterator.Valid();iterator.Next()){
            num_components++;
            T_CELL new_subcell;
            new_subcell.index=mesh.elements(template_hex_index).index;
            VECTOR<ARRAY<int>,T_MESH::vertices_per_cell> corner_fragments_i;
            for(int v=1;v<=T_MESH::vertices_per_cell;v++)
                new_subcell.vertices(v)=vertex_counter++;
            const int s=sub_cells.Append(new_subcell);
            subcell_to_root_cell.Set(s,template_hex_index);
            material_fragments.Append(iterator.Data());}
        if( (int)(((T)(last_progress+1) / mesh.elements.m) * 100) % 5 == 0  &&
            (int)(((T)(last_progress) / mesh.elements.m) * 100) % 5 != 0 ){
            LOG::cout << "Progress: "<< (int)(((T)(last_progress+1) / mesh.elements.m) * 100) << "%" << std::endl;}
        last_progress++;}
    duplicate_cell_to_root_cell=subcell_to_root_cell;

    material_fragments_hash.Exact_Resize(material_fragments.m);
    for(int i=1;i<=material_fragments.m;i++){
        HASHTABLE<int> m_fragments;
        for(int j=1;j<=material_fragments(i).m;j++)
            m_fragments.Insert(material_fragments(i)(j));
        material_fragments_hash(i)=m_fragments;}
}

template<class T,int d> bool
MATERIAL_PREDICATE_TESSELATED_SURFACE<T,d>::IsMaterialContinuous(const int axis,const int root_cell_high, const int sub_cell_high,const int root_cell_low,const int sub_cell_low) const
{
    for(int i=1;i<=material_fragments(sub_cell_high).m;i++)
        for(int j=1;j<=material_fragments(sub_cell_low).m;j++)
            if(current_tetrahedralized_volume.mesh.Face_Neighbors(material_fragments(sub_cell_high)(i),material_fragments(sub_cell_low)(j))) return true;
    return false;
}

template<class T,int d> void MATERIAL_PREDICATE_TESSELATED_SURFACE<T,d>::
MergeSubcellMaterial(const UNION_FIND<int>& merge_map)
{
    for( int i = 1; i <= material_fragments.m; i++ ){
        ARRAY<int> master_fragment_list = material_fragments( merge_map.Find(i) );
        ARRAY<int> sub_fragment_list = material_fragments( i );
        for( int j = 1; j <= sub_fragment_list.m; j++)
            master_fragment_list.Append( sub_fragment_list(j) );
        master_fragment_list.Prune_Duplicates();
        material_fragments( merge_map.Find(i) ) = master_fragment_list;
    }
    for( int i = 1; i <= material_fragments.m; i++ ){
        material_fragments(i) = material_fragments( merge_map.Find(i) );
    }
}

template<class T,int d> void MATERIAL_PREDICATE_TESSELATED_SURFACE<T,d>::
ComputeNodalDistances(const ARRAY<TV>& nodal_positions,int subcell,ARRAY<T>& distances) const
{
    for(int node=1;node<=nodal_positions.m;++node){
        const ARRAY<VECTOR<int,2> >& template_tet_nodes=template_hex_node_to_template_tet_nodes.Get(VECTOR<int,2>(duplicate_cell_to_root_cell.Get(subcell),node));
        for(int i=1;i<=template_tet_nodes.m;++i){
            int template_tet=template_tet_nodes(i).x,template_node=template_tet_nodes(i).y;
            const ARRAY<int>& duplicate_tets=template_tet_to_duplicate_tets(template_tet);
            for(int j=1;j<=duplicate_tets.m;++j) if(material_fragments_hash(subcell).Contains(duplicate_tets(j))){int duplicate_tet=duplicate_tets(j);
                const VECTOR<int,4>& element=current_tetrahedralized_volume.mesh.elements(duplicate_tet);
                int duplicate_tet_node_index=0;
                for(int t=1;t<=4;++t) if(old_particle_per_new_collapsed_particle(element(t))==template_node){duplicate_tet_node_index=t;break;}
                PHYSBAM_ASSERT(duplicate_tet_node_index!=0);
                T dist=duplicate_tet_nodal_distances(element(duplicate_tet_node_index));
                if(material_nodes_per_duplicate_tet(duplicate_tet)(duplicate_tet_node_index)) dist*=(T)-1.;
                distances(node)=min(distances(node),dist);}}}
}

template<class T,int d> bool MATERIAL_PREDICATE_TESSELATED_SURFACE<T,d>::
InsideMaterial(const int subcell_index,const int node_index) const
{
    const ARRAY<VECTOR<int,2> >& template_tet_nodes=template_hex_node_to_template_tet_nodes.Get(VECTOR<int,2>(duplicate_cell_to_root_cell.Get(subcell_index),node_index));
    for(int i=1;i<=template_tet_nodes.m;++i){
        int template_tet=template_tet_nodes(i).x,template_node=template_tet_nodes(i).y;
        const ARRAY<int>& duplicate_tets=template_tet_to_duplicate_tets(template_tet);
        for(int j=1;j<=duplicate_tets.m;++j) if(material_fragments_hash(subcell_index).Contains(duplicate_tets(j))){int duplicate_tet=duplicate_tets(j);
            const VECTOR<int,4>& element=current_tetrahedralized_volume.mesh.elements(duplicate_tet);
            int duplicate_tet_node_index=0;
            for(int t=1;t<=4;++t) if(old_particle_per_new_collapsed_particle(element(t))==template_node){duplicate_tet_node_index=t;break;}
            PHYSBAM_ASSERT(duplicate_tet_node_index!=0);
            if(material_nodes_per_duplicate_tet(duplicate_tet)(duplicate_tet_node_index)) return true;}}
    return false;
}

template<class T,int d> int MATERIAL_PREDICATE_TESSELATED_SURFACE<T,d>::
Node_Material_Representative(const int subcell_index,const int node_index) const
{
    return Subcell_Material_Representative(subcell_index);
}

template<class T,int d> int MATERIAL_PREDICATE_TESSELATED_SURFACE<T,d>::
Subcell_Material_Representative(const int subcell_index) const
{
    if(material_fragments(subcell_index).m==0) return -1;
    return material_fragments(subcell_index)(1);
}

template class MATERIAL_PREDICATE_TESSELATED_SURFACE<float,3>;
template class MATERIAL_PREDICATE_TESSELATED_SURFACE<double,3>;

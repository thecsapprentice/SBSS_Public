//#####################################################################
// Copyright 2014, Raj Setaluri, Nathan Mitchell.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class NONMANIFOLD_DUALCONTOUR_2D
//#####################################################################
#include <PhysBAM_Tools/Math_Tools/sign.h>
#include <PhysBAM_Geometry/Basic_Geometry/BOX.h>
#include <Nonmanifold_Implicit_Objects/NONMANIFOLD_DUALCONTOUR_2D.h>
#include <Nonmanifold_Implicit_Objects/NONMANIFOLD_LEVELSET_MESH.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/SEGMENTED_CURVE_2D.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_AREA.h>
#include <Common_Tools/Arrays/HYBRID_ARRAY_ITERATOR.h>
#include <PhysBAM_Tools/Arrays_Computations/SUMMATIONS.h>

using namespace PhysBAM;
//#####################################################################
// Function Generate_Topology
//#####################################################################
template<class T> void NONMANIFOLD_DUALCONTOUR_2D<T>::
Generate_Topology()
{
#if 0
    levelset.mesh.Update_Node_Neighbors();
    topology.Preallocate(5000);
    vertices.Resize(levelset.phi.Domain_Indices(),levelset.phi.Flat_Size());
    vertices.Fill(0);
    int index=1;
    for(HYBRID_ARRAY_ITERATOR<2> iter(levelset.mesh.Node_Domain(),levelset.mesh.Node_Mesh_Counts());iter.Valid();iter.Next()){
        const T_INDEX& node_index=iter.Index();
        typename NONMANIFOLD_LEVELSET_MESH<T,2>::AXIS_ALIGNED_NEIGHBORS neighbors=levelset.mesh.Node_Neighbors(node_index);
        for(int i=1;i<=neighbors(1).m;i++){
            const T_INDEX& neighbor_index=neighbors(1)(i);
            const T phi1=levelset.phi(node_index),phi2=levelset.phi(neighbor_index);
            const int index_i=node_to_int.Get_Or_Insert(node_index,index);
            if(index_i==index){int_to_node.Append(node_index);index++;}
            for(int j=1;j<=neighbors(4).m;j++){
                const T_INDEX& other_neighbor_index=neighbors(4)(j);
                const int index_j=node_to_int.Get_Or_Insert(other_neighbor_index,index);
                if(index_j==index){int_to_node.Append(other_neighbor_index);index++;}
                if(phi1<=contour_value&&phi2>contour_value) topology.Append(VECTOR<int,2>(index_j,index_i));
                if(phi1>contour_value&&phi2<=contour_value) topology.Append(VECTOR<int,2>(index_i,index_j));}}
        for(int i=1;i<=neighbors(3).m;i++){
            const T_INDEX& neighbor_index=neighbors(3)(i);
            const T phi1=levelset.phi(node_index),phi2=levelset.phi(neighbor_index);
            const int index_i=node_to_int.Get_Or_Insert(node_index,index);
            if(index_i==index){int_to_node.Append(node_index);index++;}
            for(int j=1;j<=neighbors(2).m;j++){
                const T_INDEX& other_neighbor_index=neighbors(2)(j);
                const int index_j=node_to_int.Get_Or_Insert(other_neighbor_index,index);
                if(index_j==index){int_to_node.Append(other_neighbor_index);index++;}
                if(phi1<=contour_value&&phi2>contour_value) topology.Append(VECTOR<int,2>(index_j,index_i));
                if(phi1>contour_value&&phi2<=contour_value) topology.Append(VECTOR<int,2>(index_i,index_j));}}}
    topology.Compact();
    for(int t=1;t<=topology.m;t++){
        int i,j;topology(t).Get(i,j);
        vertices(int_to_node(i))=vertices(int_to_node(j))=1;}
#endif
}
//#####################################################################
// Function Generate_Vertices
//#####################################################################
template<class T> void NONMANIFOLD_DUALCONTOUR_2D<T>::
Generate_Vertices()
{
#if 0
    const int count=vertices.Sum();
    geometry.Resize(count);
    normals.Resize(count);
    levelset.Compute_Normals();
    int vertex=0;
    for(int i=1;i<=int_to_node.m;i++){
        const T_INDEX& index=int_to_node(i);
        if(!vertices(index)) continue;
        vertices(index)=++vertex;
        const TV base_X=levelset.mesh.Node(index);
        TV weights=TV();//(T).5*TV::All_Ones_Vector(); -- TODO: figure out
        T dx=levelset.mesh.dx;

        // debug
        //geometry(vertex)=base_X+weights*dx;
        //continue;
        // debug

        TV weights_guess;
        TV normal=levelset.Normal(index,weights);
        T phi=levelset.Phi(index,weights);
        T delta_phi;
        T min_dX=dx;
        bool failed_in_positive_normal=false;
        int iterations=0;
        if(is_distance_field) while(abs(phi-contour_value)>1e-5*min_dX && (iterations++)<10){
            weights-=(phi-contour_value)*normal;
            clamp(weights,TV(),TV::All_Ones_Vector());
            phi=levelset.Phi(index,weights);
            normal=levelset.Normal(index,weights);}
        else while(abs(phi-contour_value)>1e-5*min_dX && (iterations++)<10){
            delta_phi=levelset.Phi(index,weights+normal*dx)-phi;
            if(abs(delta_phi)>1e-8*min_dX) weights_guess=weights-dx*((phi-contour_value)/delta_phi)*normal; // Newton-Raphson
            if(abs(levelset.Phi(index,weights_guess)-contour_value)>abs(phi-contour_value)){
                if(!failed_in_positive_normal){failed_in_positive_normal=true;normal=-normal;} // Look away from a shock
                else break;}
            else{
                failed_in_positive_normal=false;
                weights=weights_guess;
                phi=levelset.Phi(index,weights);
                normal=levelset.Normal(index,weights);}}
        geometry(vertex)=base_X+weights*dx;normals(vertex)=normal;}
#endif
}
//#####################################################################
// Function Ensure_Vertices_In_Correct_Cells
//#####################################################################
template<class T> void NONMANIFOLD_DUALCONTOUR_2D<T>::
Ensure_Vertices_In_Correct_Cells()
{
    // TODO: fix
    return;
    PHYSBAM_NOT_IMPLEMENTED();
#if 0
    int vertex=0;
    TV_INT i;
    for(i.x=1;i.x<grid.counts.x;i.x++) for(i.y=1;i.y<grid.counts.y;i.y++) if(vertices(i)){
        ++vertex;TV_INT v=grid.Cell(geometry(vertex),0);
        if(i!=v){
            TV cell_center=grid.Center(i);TV offset=(T).5*grid.dX;
            geometry(vertex)=BOX<TV>(cell_center-offset,cell_center+offset).Surface(geometry(vertex));}}
#endif 
}
//#####################################################################
// Function Get_Segmented_Curve
//#####################################################################
template<class T> SEGMENTED_CURVE_2D<T>* NONMANIFOLD_DUALCONTOUR_2D<T>::
Get_Segmented_Curve()
{
    SEGMENTED_CURVE_2D<T>* curve=SEGMENTED_CURVE_2D<T>::Create();
#if 0
    curve->particles.array_collection->Add_Elements(geometry.m);
    curve->particles.X=geometry;
    curve->mesh.number_nodes=geometry.m;
    curve->mesh.elements.Exact_Resize(topology.m);
    for(int t=1;t<=topology.m;t++){
        int i,j;topology(t).Get(i,j);i=vertices(int_to_node(i));j=vertices(int_to_node(j));
        curve->mesh.elements(t).Set(i,j);}
    curve->Update_Segment_List();
#endif
    return curve;
}
//#####################################################################
// Function Get_Triangulated_Surface
//#####################################################################
template<class T> TRIANGULATED_AREA<T>* NONMANIFOLD_DUALCONTOUR_2D<T>::
Get_Triangulated_Area(const int sign)
{
    PHYSBAM_NOT_IMPLEMENTED();
    TRIANGULATED_AREA<T>* area=TRIANGULATED_AREA<T>::Create();
#if 0
    GEOMETRY_PARTICLES<TV>& particles=area->particles;
    particles.array_collection->Preallocate(grid.counts.x*grid.counts.y*4);
    TRIANGLE_MESH& mesh=area->mesh;
    int triangle_count=0;for(int i=2;i<grid.counts.x;i++) for(int j=2;j<grid.counts.y;j++) if(levelset.phi(i,j)*sign>=0) triangle_count+=4;
    mesh.elements.Exact_Resize(triangle_count);
    int current_triangle=1;
    GRID<TV> mac_grid=grid.Get_MAC_Grid();
    for(int i=2;i<grid.counts.x;i++) for(int j=2;j<grid.counts.y;j++) if(levelset.phi(i,j)*sign>=0){
        int v0=particles.array_collection->Add_Element(),v1=particles.array_collection->Add_Element(),v2=particles.array_collection->Add_Element(),v3=particles.array_collection->Add_Element(),v4=particles.array_collection->Add_Element();
        particles.X(v0)=grid.X(i,j);
        particles.X(v1)=vertices(i-1,j-1)?geometry(vertices(i-1,j-1)):mac_grid.X(i-1,j-1);
        particles.X(v2)=vertices(i,j-1)?geometry(vertices(i,j-1)):mac_grid.X(i,j-1);
        particles.X(v3)=vertices(i,j)?geometry(vertices(i,j)):mac_grid.X(i,j);
        particles.X(v4)=vertices(i-1,j)?geometry(vertices(i-1,j)):mac_grid.X(i-1,j);
        mesh.elements(current_triangle++).Set(v2,v1,v0);
        mesh.elements(current_triangle++).Set(v3,v2,v0);
        mesh.elements(current_triangle++).Set(v4,v3,v0);
        mesh.elements(current_triangle++).Set(v1,v4,v0);}
    mesh.number_nodes=particles.array_collection->Size();
#endif
    return area;
}
//#####################################################################
template class NONMANIFOLD_DUALCONTOUR_2D<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class NONMANIFOLD_DUALCONTOUR_2D<double>;
#endif

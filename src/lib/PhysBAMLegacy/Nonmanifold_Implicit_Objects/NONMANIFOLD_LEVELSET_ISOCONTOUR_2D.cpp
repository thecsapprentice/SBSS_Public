//#####################################################################
// Copyright 2014, Raj Setaluri.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class NONMANIFOLD_LEVELSET_ISOCONTOUR_2D
//#####################################################################
#include <Nonmanifold_Implicit_Objects/NONMANIFOLD_LEVELSET_ISOCONTOUR_2D.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/SEGMENTED_CURVE_2D.h>
using namespace PhysBAM;
//#####################################################################
// Function Create_Segmented_Curve_From_Levelset
//#####################################################################
template<class T> SEGMENTED_CURVE_2D<T>* NONMANIFOLD_LEVELSET_ISOCONTOUR_2D<T>::
Create_Segmented_Curve_From_Levelset(const NONMANIFOLD_LEVELSET_2D<T>& levelset,const T isovalue)
{
    SEGMENTED_CURVE_2D<T>* curve=SEGMENTED_CURVE_2D< T >::Create();

    const T_NONMANIFOLD_MESH& nm_mesh=levelset.mesh;
    const typename NONMANIFOLD_LEVELSET_2D<T>::T_ARRAY_SCALAR& phi=levelset.phi;

    ARRAY<TV_INT> corner_index_offsets(T_NONMANIFOLD_MESH::nodes_per_cell);
    {int flat_index=1;
    for(RANGE_ITERATOR<2> iterator(RANGE<TV_INT>(TV_INT(),TV_INT()+1));iterator.Valid();iterator.Next())
        corner_index_offsets(flat_index++)=iterator.Index();}

    VECTOR<T,faces_per_cell> face_crossings_alpha;

    for(int cell=1;cell<=nm_mesh.cells.m;cell++){
        int key=0;
        RANGE<TV> cell_range(nm_mesh.aabb_hierarchy.box_hierarchy(cell));
        ARRAY<TV,TV_INT> corners;
        ARRAY<T,TV_INT> phi_local(RANGE<TV_INT>(TV_INT(),TV_INT()+1));
        cell_range.Corners(corners);
        T average=0;
        for(int i=1;i<=T_NONMANIFOLD_MESH::nodes_per_cell;i++){
            phi_local(corner_index_offsets(i))=phi(nm_mesh.cells(cell).nodes(i));
            if(phi(nm_mesh.cells(cell).nodes(i))>=isovalue) key=key|(1<<(i-1));
            average+=phi(nm_mesh.cells(cell).nodes(i));}
        average=average/T_NONMANIFOLD_MESH::nodes_per_cell;
        if(average>isovalue && key==6) key=9;
        if(average>isovalue && key==9) key=6;

        face_crossings_alpha(1)=-phi_local(TV_INT(1,0))/(phi_local(TV_INT(0,0))-phi_local(TV_INT(1,0)));
        face_crossings_alpha(2)=-phi_local(TV_INT(0,1))/(phi_local(TV_INT(0,0))-phi_local(TV_INT(0,1)));
        face_crossings_alpha(3)=-phi_local(TV_INT(1,1))/(phi_local(TV_INT(0,1))-phi_local(TV_INT(1,1)));
        face_crossings_alpha(4)=-phi_local(TV_INT(1,1))/(phi_local(TV_INT(1,0))-phi_local(TV_INT(1,1)));

        ARRAY<TV> points;
        switch(key){
        case 0: // 00, 00 --> No segment,fully outside
            break;
        case 1:{ // 00, 10
            TV E1=face_crossings_alpha(2)*corners(TV_INT(0,0))+((T)1.-face_crossings_alpha(2))*corners(TV_INT(0,1));
            TV E2=face_crossings_alpha(1)*corners(TV_INT(0,0))+((T)1.-face_crossings_alpha(1))*corners(TV_INT(1,0));
            points.Append(E1);
            points.Append(E2);}
            break;
        case 2:{ // 10, 00
            TV E1=face_crossings_alpha(2)*corners(TV_INT(0,0))+((T)1.-face_crossings_alpha(2))*corners(TV_INT(0,1));
            TV E2=face_crossings_alpha(3)*corners(TV_INT(0,1))+((T)1.-face_crossings_alpha(3))*corners(TV_INT(1,1));
            points.Append(E1);
            points.Append(E2);}
            break;
        case 3:{ // 10, 10
            TV E1=face_crossings_alpha(1)*corners(TV_INT(0,0))+((T)1.-face_crossings_alpha(1))*corners(TV_INT(1,0));
            TV E2=face_crossings_alpha(3)*corners(TV_INT(0,1))+((T)1.-face_crossings_alpha(3))*corners(TV_INT(1,1));
            points.Append(E1);
            points.Append(E2);}
            break;            
        case 4:{ // 00, 01
            TV E1=face_crossings_alpha(1)*corners(TV_INT(0,0))+((T)1.-face_crossings_alpha(1))*corners(TV_INT(1,0));
            TV E2=face_crossings_alpha(4)*corners(TV_INT(1,0))+((T)1.-face_crossings_alpha(4))*corners(TV_INT(1,1));
            points.Append(E1);
            points.Append(E2);}
            break;
            
        case 5:{ // 00, 11
            TV E1=face_crossings_alpha(2)*corners(TV_INT(0,0))+((T)1.-face_crossings_alpha(2))*corners(TV_INT(0,1));
            TV E2=face_crossings_alpha(4)*corners(TV_INT(1,0))+((T)1.-face_crossings_alpha(4))*corners(TV_INT(1,1));
            points.Append(E1);
            points.Append(E2);}
            break;
            
        case 6:{ // 10, 01
            TV E1=face_crossings_alpha(2)*corners(TV_INT(0,0))+((T)1.-face_crossings_alpha(2))*corners(TV_INT(0,1));
            TV E2=face_crossings_alpha(3)*corners(TV_INT(0,1))+((T)1.-face_crossings_alpha(3))*corners(TV_INT(1,1));
            points.Append(E1);
            points.Append(E2);}
            {TV E1=face_crossings_alpha(1)*corners(TV_INT(0,0))+((T)1.-face_crossings_alpha(1))*corners(TV_INT(1,0));
            TV E2=face_crossings_alpha(4)*corners(TV_INT(1,0))+((T)1.-face_crossings_alpha(4))*corners(TV_INT(1,1));
            points.Append(E1);
            points.Append(E2);}
            break;
        case 7:{ // 10, 11
            TV E1=face_crossings_alpha(3)*corners(TV_INT(0,1))+((T)1.-face_crossings_alpha(3))*corners(TV_INT(1,1));
            TV E2=face_crossings_alpha(4)*corners(TV_INT(1,0))+((T)1.-face_crossings_alpha(4))*corners(TV_INT(1,1));
            points.Append(E1);
            points.Append(E2);}
            break;
        case 8:{ // 01, 00
            TV E1=face_crossings_alpha(3)*corners(TV_INT(0,1))+((T)1.-face_crossings_alpha(3))*corners(TV_INT(1,1));
            TV E2=face_crossings_alpha(4)*corners(TV_INT(1,0))+((T)1.-face_crossings_alpha(4))*corners(TV_INT(1,1));
            points.Append(E1);
            points.Append(E2);}
            break;
        case 9:{ // 01, 10
            TV E1=face_crossings_alpha(3)*corners(TV_INT(0,1))+((T)1.-face_crossings_alpha(3))*corners(TV_INT(1,1));
            TV E2=face_crossings_alpha(4)*corners(TV_INT(1,0))+((T)1.-face_crossings_alpha(4))*corners(TV_INT(1,1));
            points.Append(E1);
            points.Append(E2);}
            {TV E1=face_crossings_alpha(2)*corners(TV_INT(0,0))+((T)1.-face_crossings_alpha(2))*corners(TV_INT(0,1));
            TV E2=face_crossings_alpha(1)*corners(TV_INT(0,0))+((T)1.-face_crossings_alpha(1))*corners(TV_INT(1,0));
            points.Append(E1);
            points.Append(E2);}
            break;
        case 10:{ // 11, 00
            TV E1=face_crossings_alpha(2)*corners(TV_INT(0,0))+((T)1.-face_crossings_alpha(2))*corners(TV_INT(0,1));
            TV E2=face_crossings_alpha(4)*corners(TV_INT(1,0))+((T)1.-face_crossings_alpha(4))*corners(TV_INT(1,1));
            points.Append(E1);
            points.Append(E2);}
            break;
        case 11:{ // 11, 10
            TV E1=face_crossings_alpha(1)*corners(TV_INT(0,0))+((T)1.-face_crossings_alpha(1))*corners(TV_INT(1,0));
            TV E2=face_crossings_alpha(4)*corners(TV_INT(1,0))+((T)1.-face_crossings_alpha(4))*corners(TV_INT(1,1));
            points.Append(E1);
            points.Append(E2);}
            break;
        case 12:{ // 01, 01
            TV E1=face_crossings_alpha(1)*corners(TV_INT(0,0))+((T)1.-face_crossings_alpha(1))*corners(TV_INT(1,0));
            TV E2=face_crossings_alpha(3)*corners(TV_INT(0,1))+((T)1.-face_crossings_alpha(3))*corners(TV_INT(1,1));
            points.Append(E1);
            points.Append(E2);}
            break;
        case 13:{ // 01, 11
            TV E1=face_crossings_alpha(2)*corners(TV_INT(0,0))+((T)1.-face_crossings_alpha(2))*corners(TV_INT(0,1));
            TV E2=face_crossings_alpha(3)*corners(TV_INT(0,1))+((T)1.-face_crossings_alpha(3))*corners(TV_INT(1,1));
            points.Append(E1);
            points.Append(E2);}
            break;
        case 14:{ // 11, 01
            TV E1=face_crossings_alpha(1)*corners(TV_INT(0,0))+((T)1.-face_crossings_alpha(1))*corners(TV_INT(1,0));
            TV E2=face_crossings_alpha(2)*corners(TV_INT(0,0))+((T)1.-face_crossings_alpha(2))*corners(TV_INT(0,1));
            points.Append(E1);
            points.Append(E2);}
            break;
        case 15: // 11, 11 --> No contour,fully inside
            break;}
        for(int p=1;p <= points.m;p+=2){
            int p1=curve->particles.array_collection->Add_Element();
            int p2=curve->particles.array_collection->Add_Element();
            curve->particles.X(p1)=points(p);
            curve->particles.X(p2)=points(p+1);
            curve->mesh.elements.Append(VECTOR<int,2>(p1,p2));}}
    curve->Update_Number_Nodes();
    curve->Update_Bounding_Box();
    return curve;
}
//#####################################################################
template class NONMANIFOLD_LEVELSET_ISOCONTOUR_2D<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class NONMANIFOLD_LEVELSET_ISOCONTOUR_2D<double>;
#endif

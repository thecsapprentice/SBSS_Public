//#####################################################################
// Copyright (c) 2014, Mridul Aanjaneya.
// Distributed under the FreeBSD license (see license.txt)
//#####################################################################
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_PREFERENCES.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_UNIFORM_SLICE.h>
#include <PhysBAM_Tools/Arrays_Computations/ARRAY_MIN_MAX.h>
#include <PhysBAM_Tools/Read_Write/Arrays/READ_WRITE_ARRAY.h>
#include <PhysBAM_Tools/Read_Write/Data_Structures/READ_WRITE_PAIR.h>
#include <PhysBAM_Tools/Read_Write/Grids_Uniform_Arrays/READ_WRITE_ARRAYS.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <Common_Tools/Math_Tools/RANGE_ITERATOR.h>
#include <Common_Tools/Read_Write/READ_WRITE_HYBRID_ARRAY.h>
#include <Rendering/OpenGL/OPENGL_NONMANIFOLD_LEVELSET_MESH_3D.h>
#include <iostream>
using namespace PhysBAM;
//#####################################################################
// Function Display
//#####################################################################
template<class T> void OPENGL_NONMANIFOLD_LEVELSET_MESH_3D<T>::
Display(const int in_color) const
{
    if(!draw) return;
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    Send_Transform_To_GL_Pipeline();

    glPushAttrib(GL_ENABLE_BIT);
    glDisable(GL_LIGHTING);

    glDisable(GL_CULL_FACE);
    //glColorMaterial(GL_FRONT_AND_BACK, GL_DIFFUSE);
    //glEnable(GL_COLOR_MATERIAL);

    GLint mode;
    glGetIntegerv(GL_RENDER_MODE,&mode);
    OPENGL_UNIFORM_SLICE* slice=(OPENGL_UNIFORM_SLICE*)this->slice;
    int axis=0;TV location;
    if(slice){axis=slice->axis;location=grid.domain.min_corner+grid.dX*(TV::Axis_Vector(axis)*(slice->index)-TV::Axis_Vector(axis)*(T).5);}
    
    VECTOR<int,4> node_ordering;
    node_ordering(1)=1;
    node_ordering(2)=2;
    node_ordering(3)=4;
    node_ordering(4)=3;
    VECTOR<T_INDEX,8> node_map;
    node_map(1)=T_INDEX(0,0,0);
    node_map(2)=T_INDEX(0,0,1);
    node_map(3)=T_INDEX(0,1,0);
    node_map(4)=T_INDEX(0,1,1);
    node_map(5)=T_INDEX(1,0,0);
    node_map(6)=T_INDEX(1,0,1);
    node_map(7)=T_INDEX(1,1,0);
    node_map(8)=T_INDEX(1,1,1);

    if(slice && slice->mode==OPENGL_SLICE::CELL_SLICE){
        if(mode == GL_SELECT){if(!collapse_slice){
            glPushAttrib(GL_ENABLE_BIT | GL_POINT_BIT);
            glPointSize(OPENGL_PREFERENCES::selection_point_size);

            // Draw nodes for selection
            glPushName(1);
            for(int i=1;i<=cells->Size();++i){
                TV left_corner((*cells)(i).cell_index-T_INDEX::All_Ones_Vector());
                if(location[axis]>=left_corner[axis]*grid.dX[axis] && location[axis]<=((left_corner+TV::All_Ones_Vector())*grid.dX)[axis]){
                for(int j=1;j<=NONMANIFOLD_LEVELSET_MESH<T,3>::nodes_per_cell;++j){VECTOR<T,3> X;
                    T_INDEX index=(*cells)(i).nodes(j).x;
                    if((*cells)(i).nodes(j).y>0) index=(*cells)(i).cell_index+node_map(j);
                    for(int v=1;v<=3;++v) X(v)=grid.Node(index)(v);
                    X(axis)+=(mesh_height-mesh_height*(*node_heights)((*cells)(i).nodes(j)));
                    glPushName((*cells)(i).nodes(j).y);
                    glPushName(index(1));glPushName(index(2));glPushName(index(3));
                    OPENGL_COLOR::Gray(.5).Send_To_GL_Pipeline();
                    OpenGL_Begin(GL_POINTS);
                    OpenGL_Vertex(X);
                    OpenGL_End();
                    glPopName();glPopName();glPopName();
                    glPopName();}}}
            glPopName();
            glPopAttrib();}}
        else{if(!collapse_slice){
            for(int i=1;i<=cells->Size();++i){TV left_corner((*cells)(i).cell_index-T_INDEX::All_Ones_Vector());
                if(location[axis]>=left_corner[axis]*grid.dX[axis] && location[axis]<=((left_corner+TV::All_Ones_Vector())*grid.dX)[axis]){
                    ARRAY<VECTOR<T,3> > node_locations;
                    glPushAttrib(GL_ENABLE_BIT | GL_POINT_BIT);
                    OPENGL_COLOR::Gray(.5).Send_To_GL_Pipeline();
                    OpenGL_Begin(GL_POINTS);
                    for(int j=1;j<=NONMANIFOLD_LEVELSET_MESH<T,3>::nodes_per_cell;++j){VECTOR<T,3> X;
                        T_INDEX index=(*cells)(i).nodes(j).x;
                        if((*cells)(i).nodes(j).y>0) index=(*cells)(i).cell_index+node_map(j);
                        for(int v=1;v<=3;++v) X(v)=grid.Node(index)(v);
                        X(axis)+=(mesh_height-mesh_height*(*node_heights)((*cells)(i).nodes(j)));
                        OpenGL_Vertex(X);
                        node_locations.Append(X);}
                    OpenGL_End();
                    glPopAttrib();
                    if((*cell_is_grid)(i)) OPENGL_COLOR::Gray(.5).Send_To_GL_Pipeline();
                    else OPENGL_COLOR::Blue().Send_To_GL_Pipeline();
                    OpenGL_Begin(GL_LINES);
                    for(int j=1;j<=3;++j) for(int k=1;k<=4;++k){
                        int node_A_order=node_ordering(k),node_B_order=node_ordering(k%4+1);
                        if(j==2){node_A_order=node_ordering(k)+4;node_B_order=node_ordering(k%4+1)+4;}
                        else if(j==3){node_A_order=node_ordering(k);node_B_order=node_ordering(k)+4;}
                        T nodeA_factor = (*phi)((*cells)(i).nodes(node_A_order));
                        T nodeB_factor = (*phi)((*cells)(i).nodes(node_B_order));
                        cramp->Lookup(nodeA_factor).Send_To_GL_Pipeline();
                        OpenGL_Vertex(node_locations(node_A_order));
                        cramp->Lookup(nodeB_factor).Send_To_GL_Pipeline();
                        OpenGL_Vertex(node_locations(node_B_order));}
                    OpenGL_End();}}}
            else{for(int i=1;i<=cells->Size();++i){TV left_corner((*cells)(i).cell_index-T_INDEX::All_Ones_Vector());
                if(location[axis]>=left_corner[axis]*grid.dX[axis] && location[axis]<=((left_corner+TV::All_Ones_Vector())*grid.dX)[axis]){
                    glPushAttrib(GL_ENABLE_BIT | GL_POINT_BIT);
                    OPENGL_COLOR::Gray(.5).Send_To_GL_Pipeline();
                    ARRAY<VECTOR<T,3> > node_locations;
                    OpenGL_Begin(GL_POINTS);
                    for(int t=1;t<=4;++t){VECTOR<T,3> X1,X2;int t1=t,t2=t+4;
                        if(axis==2){t1=(t<=2)?t:t+2;t2=t1+2;}
                        else if(axis==3){t1=2*t-1;t2=t1+1;}
                        T_INDEX index1=(*cells)(i).nodes(t1).x,index2=(*cells)(i).nodes(t2).x;
                        if((*cells)(i).nodes(t1).y>0) index1=(*cells)(i).cell_index+node_map(t1);
                        if((*cells)(i).nodes(t2).y>0) index2=(*cells)(i).cell_index+node_map(t2);
                        for(int v=1;v<=3;++v){X1(v)=grid.Node(index1)(v);X2(v)=grid.Node(index2)(v);}
                        X1(axis)+=(mesh_height-mesh_height*(*node_heights)((*cells)(i).nodes(t1)));
                        X2(axis)+=(mesh_height-mesh_height*(*node_heights)((*cells)(i).nodes(t2)));
                        VECTOR<T,3> X=(T).5*(X1+X2);
                        OpenGL_Vertex(X);
                        node_locations.Append(X);}
                    OpenGL_End();
                    glPopAttrib();
                    if((*cell_is_grid)(i)) OPENGL_COLOR::Gray(.5).Send_To_GL_Pipeline();
                    else OPENGL_COLOR::Blue().Send_To_GL_Pipeline();
                    OpenGL_Begin(GL_LINES);
                    for(int k=1;k<=4;++k){int node_A_order=node_ordering(k),node_B_order=node_ordering(k%4+1);
                        T nodeA_factor = (*phi)((*cells)(i).nodes(node_A_order));
                        T nodeB_factor = (*phi)((*cells)(i).nodes(node_B_order));
                        cramp->Lookup(nodeA_factor).Send_To_GL_Pipeline();
                        OpenGL_Vertex(node_locations(node_A_order));
                        cramp->Lookup(nodeB_factor).Send_To_GL_Pipeline();
                        OpenGL_Vertex(node_locations(node_B_order));}
                    OpenGL_End();}}}}}

    if(current_selection){
        if(current_selection->type == OPENGL_SELECTION::GRID_NODE_3D){
            OPENGL_SELECTION_NONMANIFOLD_LEVELSET_MESH_NODE_3D<T> *real_selection = (OPENGL_SELECTION_NONMANIFOLD_LEVELSET_MESH_NODE_3D<T>*)current_selection;
            VECTOR<T,3> X;
            for(int v=1;v<=3;++v) X(v)=grid.Node(real_selection->index)(v);
            X(axis)+=(mesh_height-mesh_height*(*node_heights)(real_selection->hindex));
            OPENGL_SELECTION::Draw_Highlighted_Vertex(X);}}

    glPopAttrib();
    glPopMatrix();
}
//#####################################################################
// Function Initialize
//#####################################################################
template<class T> void OPENGL_NONMANIFOLD_LEVELSET_MESH_3D<T>::
Initialize(const std::string& filename,const std::string& flag_filename, const std::string& value_filename, const int frame)
{
    if(cells) delete cells;
    std::string directory_filename=FILE_UTILITIES::Get_Frame_Filename(filename,frame);
    int size=0;
    std::istream* input=FILE_UTILITIES::Safe_Open_Input(directory_filename);
    Read_Binary<RW>(*input,size);
    cells=new ARRAY<T_CELL>(size);
    for(int i=1;i<=size;++i){
        Read_Binary<RW>(*input,(*cells)(i).cell_index);
        for(int j=1;j<=NONMANIFOLD_LEVELSET_MESH<T,3>::nodes_per_cell;++j)
            Read_Binary<RW>(*input,(*cells)(i).nodes(j));}
    if(cell_is_grid) delete cell_is_grid;
    cell_is_grid=new ARRAY<bool>(cells->Size());
    FILE_UTILITIES::Read_From_File(stream_type,flag_filename,*cell_is_grid);
    if(phi) delete phi;
    phi = new HYBRID_ARRAY<T,3>();
    FILE_UTILITIES::Read_From_File(stream_type,value_filename,*phi);
    

    cramp = new OPENGL_COLOR_RAMP<T>;
    T min_value = FLT_MAX;
    T max_value = -FLT_MAX;
    for( int i = 1; i <= cells->m; i++){
        int j=1;
        for( RANGE_ITERATOR<3> iter( RANGE<T_INDEX>((*cells)(i).cell_index, (*cells)(i).cell_index +1 )); iter.Valid(); iter.Next(), j++){
            T v = (*phi)( (*cells)(i).nodes(j) );
            ARRAY<HINDEX> overlapping_nodes = node_overlaps.Get_Default( iter.Index(), ARRAY<HINDEX>() );
            if( v > max_value )
                max_value = v;
            if( v < min_value )
                min_value = v;
            overlapping_nodes.Append_Unique( (*cells)(i).nodes(j) );
            node_overlaps.Set( iter.Index(), overlapping_nodes );
        }
    }
    node_heights = new HYBRID_ARRAY<T,3>(phi->Domain_Indices(), phi->Flat_Size() );
    Reset_Node_Height();

    T interval_width_inside = 0-min_value;
    T interval_width_outside = max_value-0;
    // Inside_Colors
    cramp->Add_Color(interval_width_inside*0+min_value,OPENGL_COLOR(0,0,0));
    cramp->Add_Color(interval_width_inside*0.3750+min_value,OPENGL_COLOR(1,0,0));
    cramp->Add_Color(interval_width_inside*0.7656+min_value,OPENGL_COLOR(1,1,0));
    cramp->Add_Color(interval_width_inside*1+min_value,OPENGL_COLOR(1,1,1));
    // Outside Colors
    cramp->Add_Color(interval_width_outside*0,OPENGL_COLOR(1,1,1));
    cramp->Add_Color(interval_width_inside*0.7656,OPENGL_COLOR(0,1,1));
    cramp->Add_Color(interval_width_inside*0.3750,OPENGL_COLOR(0,0,1));
    cramp->Add_Color(interval_width_outside*1,OPENGL_COLOR(0,0,0));

    //::Matlab_Hot(min(ARRAYS_COMPUTATIONS::Min(phi->flat_array),phi->nd_array.Min()),
    //                                         max(ARRAYS_COMPUTATIONS::Max(phi->flat_array),phi->nd_array.Max()));
}
/*
//#####################################################################
// Function Set_Frame
//#####################################################################
template<class T> void OPENGL_NONMANIFOLD_LEVELSET_MESH_2D<T>::
Increase_Node_Height(){

    if(current_selection){
        if(current_selection->type == OPENGL_SELECTION::GRID_NODE_3D){
            OPENGL_SELECTION_NONMANIFOLD_LEVELSET_MESH_NODE_2D<T> *real_selection = (OPENGL_SELECTION_NONMANIFOLD_LEVELSET_MESH_NODE_2D<T>*)current_selection;
            T interval = 1.0f / T( node_overlaps.Get(real_selection->index).m);
            (*node_heights)(real_selection->hindex) += interval;
        }
    }

}
//#####################################################################
// Function Set_Frame
//#####################################################################
template<class T> void OPENGL_NONMANIFOLD_LEVELSET_MESH_2D<T>::
Decrease_Node_Height(){

    if(current_selection){
        if(current_selection->type == OPENGL_SELECTION::GRID_NODE_3D){
            OPENGL_SELECTION_NONMANIFOLD_LEVELSET_MESH_NODE_2D<T> *real_selection = (OPENGL_SELECTION_NONMANIFOLD_LEVELSET_MESH_NODE_2D<T>*)current_selection;
            T interval = 1.0f / T( node_overlaps.Get(real_selection->index).m);
            (*node_heights)(real_selection->hindex) -= interval;
        }
    }
 
}
*/
//#####################################################################
// Function Set_Frame
//#####################################################################
template<class T> void OPENGL_NONMANIFOLD_LEVELSET_MESH_3D<T>::
Reset_Node_Height(){
    for( int i = 1; i <= cells->m; i++){
        int j=1;
        for( RANGE_ITERATOR<3> iter( RANGE<T_INDEX>((*cells)(i).cell_index, (*cells)(i).cell_index +1 )); iter.Valid(); iter.Next(), j++){
            T_INDEX index = iter.Index();
            HINDEX hindex = (*cells)(i).nodes(j);
            T node_height =( T(node_overlaps.Get(index).Find(hindex)) 
                             / T( node_overlaps.Get(index).m) );
            (*node_heights)(hindex) = node_height;

        }
    }
}
//#####################################################################
// Function Set_Frame
//#####################################################################
template<class T> void OPENGL_NONMANIFOLD_LEVELSET_MESH_3D<T>::
Set_Frame(int frame_input)
{
    frame=frame_input;
    return;
}
//#####################################################################
// Function Bounding_Box
//#####################################################################
template<class T> RANGE<VECTOR<float,3> > OPENGL_NONMANIFOLD_LEVELSET_MESH_3D<T>::
Bounding_Box() const
{
    return World_Space_Box(RANGE<VECTOR<float,3> >(grid.domain));
}
//#####################################################################
// Function Get_Selection
//#####################################################################
template<class T> OPENGL_SELECTION *OPENGL_NONMANIFOLD_LEVELSET_MESH_3D<T>::
Get_Selection(GLuint* buffer,int buffer_size)
{
    OPENGL_SELECTION* selection=0;
    if(buffer_size==5){
        if(buffer[0]==1) selection=new OPENGL_SELECTION_NONMANIFOLD_LEVELSET_MESH_NODE_3D<T>(this,int(buffer[1]),VECTOR<int,3>(buffer[2],buffer[3],buffer[4]));}
    return selection;
}
//#####################################################################
// Function Highlight_Selection
//#####################################################################
template<class T> void OPENGL_NONMANIFOLD_LEVELSET_MESH_3D<T>::
Highlight_Selection(OPENGL_SELECTION *selection)
{
    delete current_selection;current_selection=0;
    if (selection->type == OPENGL_SELECTION::GRID_NODE_3D){
        OPENGL_SELECTION_NONMANIFOLD_LEVELSET_MESH_NODE_3D<T>* real_selection=(OPENGL_SELECTION_NONMANIFOLD_LEVELSET_MESH_NODE_3D<T>*)selection;
        current_selection=new OPENGL_SELECTION_NONMANIFOLD_LEVELSET_MESH_NODE_3D<T>(this,real_selection->grid_mapped,real_selection->index);}
}
//#####################################################################
// Function Clear_Highlight
//#####################################################################
template<class T> void OPENGL_NONMANIFOLD_LEVELSET_MESH_3D<T>::
Clear_Highlight()
{
    delete current_selection;current_selection=0;
}
//#####################################################################
// Print_Selection_Info
//#####################################################################
template<class T> void OPENGL_NONMANIFOLD_LEVELSET_MESH_3D<T>::
Print_Selection_Info(std::ostream& stream,OPENGL_SELECTION* selection) const
{
    if(current_selection && current_selection->type==OPENGL_SELECTION::GRID_NODE_3D){
        int grid_mapped=((OPENGL_SELECTION_NONMANIFOLD_LEVELSET_MESH_NODE_3D<T>*)current_selection)->grid_mapped==0?1:0;
        int mesh_index=((OPENGL_SELECTION_NONMANIFOLD_LEVELSET_MESH_NODE_3D<T>*)current_selection)->grid_mapped;
        VECTOR<int,3> index=((OPENGL_SELECTION_NONMANIFOLD_LEVELSET_MESH_NODE_3D<T>*)current_selection)->index;
        stream<<"Grid-Mapped: "<<grid_mapped<<std::endl<<"Selected node: ["<<index<<","<<mesh_index<<"] ("<<grid.Node(index)<<")"<<std::endl<<"Phi Value: "<<(*phi)(HINDEX(index,mesh_index))<<std::endl;}
}
//#####################################################################
// Function Bounding_Box
//#####################################################################
template<class T> RANGE<VECTOR<float,3> > OPENGL_SELECTION_NONMANIFOLD_LEVELSET_MESH_NODE_3D<T>::
Bounding_Box() const
{
    PHYSBAM_ASSERT(object);
    const GRID<TV> &grid=((OPENGL_NONMANIFOLD_LEVELSET_MESH_3D<T> *)object)->grid;
    RANGE<VECTOR<T,3> > box(grid.Node(index));
    return object->World_Space_Box(box);
}
//#####################################################################
template class OPENGL_NONMANIFOLD_LEVELSET_MESH_3D<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class OPENGL_NONMANIFOLD_LEVELSET_MESH_3D<double>;
#endif

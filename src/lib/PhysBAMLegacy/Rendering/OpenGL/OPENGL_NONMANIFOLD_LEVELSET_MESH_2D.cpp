//#####################################################################
// Copyright (c) 2014, Mridul Aanjaneya.
// Distributed under the FreeBSD license (see license.txt)
//#####################################################################
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Tools/Read_Write/Arrays/READ_WRITE_ARRAY.h>
#include <PhysBAM_Tools/Read_Write/Data_Structures/READ_WRITE_PAIR.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_PREFERENCES.h>
#include <Rendering/OpenGL/OPENGL_NONMANIFOLD_LEVELSET_MESH_2D.h>
#include <PhysBAM_Tools/Read_Write/Grids_Uniform_Arrays/READ_WRITE_ARRAYS.h>
#include <Common_Tools/Read_Write/READ_WRITE_HYBRID_ARRAY.h>
#include <PhysBAM_Tools/Arrays_Computations/ARRAY_MIN_MAX.h>
#include <Common_Tools/Math_Tools/RANGE_ITERATOR.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_SHAPES.h>
#include <iostream>
using namespace PhysBAM;
//#####################################################################
// Function Display
//#####################################################################
template<class T> void OPENGL_NONMANIFOLD_LEVELSET_MESH_2D<T>::
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
    
    VECTOR<int,4> node_ordering;
    node_ordering(1)=1;
    node_ordering(2)=2;
    node_ordering(3)=4;
    node_ordering(4)=3;
    VECTOR<T_INDEX,4> node_map;
    node_map(1)=T_INDEX(0,0);
    node_map(2)=T_INDEX(0,1);
    node_map(3)=T_INDEX(1,0);
    node_map(4)=T_INDEX(1,1);

    if(mode == GL_SELECT){
        glPushAttrib(GL_ENABLE_BIT | GL_POINT_BIT);
        glPointSize(OPENGL_PREFERENCES::selection_point_size);

        // Draw nodes for selection
        glPushName(1);
        for(int i=1;i<=cells->Size();++i){
            for(int j=1;j<=NONMANIFOLD_LEVELSET_MESH<T,2>::nodes_per_cell;++j){VECTOR<T,3> X;
                T_INDEX index=(*cells)(i).nodes(j).x;
                if((*cells)(i).nodes(j).y>0)
                    index = (*cells)(i).cell_index+node_map(j);
                for(int axis=1;axis<=2;++axis) X(axis)=grid.Node(index)(axis);
                X(3)=mesh_height - mesh_height*(*node_heights)((*cells)(i).nodes(j));
                glPushName((*cells)(i).nodes(j).y);
                glPushName(index(1));glPushName(index(2));
                OPENGL_COLOR::Gray(.5).Send_To_GL_Pipeline();
                OpenGL_Begin(GL_POINTS);
                OpenGL_Vertex(X);
                OpenGL_End();
                glPopName();glPopName();
                glPopName();}}
        glPopName();
        glPopAttrib();}
    else{for(int i=1;i<=cells->Size();++i){
            ARRAY<VECTOR<T,3> > node_locations;
            glPushAttrib(GL_ENABLE_BIT | GL_POINT_BIT);
            OPENGL_COLOR::Gray(.5).Send_To_GL_Pipeline();
            OpenGL_Begin(GL_POINTS);
            for(int j=1;j<=NONMANIFOLD_LEVELSET_MESH<T,2>::nodes_per_cell;++j){VECTOR<T,3> X;
                T_INDEX index=(*cells)(i).nodes(j).x;
                if((*cells)(i).nodes(j).y>0)
                    index = (*cells)(i).cell_index+node_map(j);
                for(int axis=1;axis<=2;++axis) X(axis)=grid.Node(index)(axis);
                X(3)=mesh_height - mesh_height*(*node_heights)((*cells)(i).nodes(j));
                OpenGL_Vertex(X);
                node_locations.Append(X);}
            OpenGL_End();
            glPopAttrib();
            if((*cell_is_grid)(i)) OPENGL_COLOR::Gray(.5).Send_To_GL_Pipeline();
            else OPENGL_COLOR::Blue().Send_To_GL_Pipeline();
            if(draw_lines)
                OpenGL_Begin(GL_LINES);
            else
                OpenGL_Begin(GL_QUADS);
            for(int j=1;j<=node_locations.m;++j){
                int node_A_order = node_ordering(j);
                int node_B_order = node_ordering(j%node_locations.m+1);

                T nodeA_factor = (*phi)((*cells)(i).nodes(node_A_order));
                T nodeB_factor = (*phi)((*cells)(i).nodes(node_B_order));
                cramp->Lookup(nodeA_factor).Send_To_GL_Pipeline();
                OpenGL_Vertex(node_locations(node_A_order));
                cramp->Lookup(nodeB_factor).Send_To_GL_Pipeline();
                OpenGL_Vertex(node_locations(node_B_order));
                //OpenGL_Line<T,3>(node_locations(node_ordering(j)),node_locations(node_ordering(j%node_locations.m+1)));}
            }
            OpenGL_End();}}

    if(draw_normals){
        OpenGL_Begin(GL_LINES);
        for(int i=1;i<=cells->Size();++i){
            const TV normal=levelset->Lazy_Normal(i,(T).5*TV::All_Ones_Vector());
            const TV X=grid.Center((*cells)(i).cell_index);
            OPENGL_SHAPES::Draw_Arrow(X,X+(T)vector_size*normal);}
        OpenGL_End();}

    if(current_selection){
        if(current_selection->type == OPENGL_SELECTION::GRID_NODE_3D){
            OPENGL_SELECTION_NONMANIFOLD_LEVELSET_MESH_NODE_2D<T> *real_selection = (OPENGL_SELECTION_NONMANIFOLD_LEVELSET_MESH_NODE_2D<T>*)current_selection;
            VECTOR<T,3> X;
            for(int axis=1;axis<=2;++axis) X(axis)=grid.Node(real_selection->index)(axis);
            X(3)=mesh_height - mesh_height*(*node_heights)(real_selection->hindex);
            OPENGL_SELECTION::Draw_Highlighted_Vertex(X);}
    }

    glPopAttrib();
    glPopMatrix();
}
//#####################################################################
// Function Initialize
//#####################################################################
template<class T> void OPENGL_NONMANIFOLD_LEVELSET_MESH_2D<T>::
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
        for(int j=1;j<=NONMANIFOLD_LEVELSET_MESH<T,2>::nodes_per_cell;++j)
            Read_Binary<RW>(*input,(*cells)(i).nodes(j));}
    if(cell_is_grid) delete cell_is_grid;
    cell_is_grid=new ARRAY<bool>(cells->Size());
    FILE_UTILITIES::Read_From_File(stream_type,flag_filename,*cell_is_grid);
    if(phi) delete phi;
    phi = new HYBRID_ARRAY<T,2>();
    FILE_UTILITIES::Read_From_File(stream_type,value_filename,*phi);

    if(levelset_mesh) delete levelset_mesh;
    levelset_mesh=new NONMANIFOLD_LEVELSET_MESH<T,2>(grid.dX.x);
    levelset_mesh->cells=(*cells);
    if(levelset) delete levelset;
    levelset=new NONMANIFOLD_LEVELSET_2D<T>(*levelset_mesh,*phi);
    

    cramp = new OPENGL_COLOR_RAMP<T>;
    T min_value = FLT_MAX;
    T max_value = -FLT_MAX;
    for( int i = 1; i <= cells->m; i++){
        int j=1;
        for( RANGE_ITERATOR<2> iter( RANGE<T_INDEX>((*cells)(i).cell_index, (*cells)(i).cell_index +1 )); iter.Valid(); iter.Next(), j++){
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
    node_heights = new HYBRID_ARRAY<T,2>(phi->Domain_Indices(), phi->Flat_Size() );
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
//#####################################################################
// Function Set_Frame
//#####################################################################
template<class T> void OPENGL_NONMANIFOLD_LEVELSET_MESH_2D<T>::
Reset_Node_Height(){
    for( int i = 1; i <= cells->m; i++){
        int j=1;
        for( RANGE_ITERATOR<2> iter( RANGE<T_INDEX>((*cells)(i).cell_index, (*cells)(i).cell_index +1 )); iter.Valid(); iter.Next(), j++){
            T_INDEX index = iter.Index();
            HINDEX hindex = (*cells)(i).nodes(j);
            T node_height =( T(node_overlaps.Get(index).Find(hindex)) 
                             / T( node_overlaps.Get(index).m) );
            (*node_heights)(hindex) = node_height;

        }
    }
}
//#####################################################################
// Function Increase_Vector_Size
//#####################################################################
template<class T> void OPENGL_NONMANIFOLD_LEVELSET_MESH_2D<T>::
Increase_Vector_Size()
{
    vector_size*=(T)1.1;
}
//#####################################################################
// Function Decrease_Vector_Size
//#####################################################################
template<class T> void OPENGL_NONMANIFOLD_LEVELSET_MESH_2D<T>::
Decrease_Vector_Size()
{
    vector_size/=(T)1.1;
}
//#####################################################################
// Function Set_Frame
//#####################################################################
template<class T> void OPENGL_NONMANIFOLD_LEVELSET_MESH_2D<T>::
Set_Frame(int frame_input)
{
    frame=frame_input;
    return;
}
//#####################################################################
// Function Bounding_Box
//#####################################################################
template<class T> RANGE<VECTOR<float,3> > OPENGL_NONMANIFOLD_LEVELSET_MESH_2D<T>::
Bounding_Box() const
{
    return World_Space_Box(RANGE<VECTOR<float,2> >(grid.domain));
}
//#####################################################################
// Function Get_Selection
//#####################################################################
template<class T> OPENGL_SELECTION *OPENGL_NONMANIFOLD_LEVELSET_MESH_2D<T>::
Get_Selection(GLuint* buffer,int buffer_size)
{
    OPENGL_SELECTION* selection=0;
    if(buffer_size==4){
        if(buffer[0]==1) selection=new OPENGL_SELECTION_NONMANIFOLD_LEVELSET_MESH_NODE_2D<T>(this,int(buffer[1]),VECTOR<int,2>(buffer[2],buffer[3]));
    }
    return selection;
}
//#####################################################################
// Function Highlight_Selection
//#####################################################################
template<class T> void OPENGL_NONMANIFOLD_LEVELSET_MESH_2D<T>::
Highlight_Selection(OPENGL_SELECTION *selection)
{
    delete current_selection;current_selection=0;
    if (selection->type == OPENGL_SELECTION::GRID_NODE_3D){
        OPENGL_SELECTION_NONMANIFOLD_LEVELSET_MESH_NODE_2D<T>* real_selection=(OPENGL_SELECTION_NONMANIFOLD_LEVELSET_MESH_NODE_2D<T>*)selection;
        current_selection=new OPENGL_SELECTION_NONMANIFOLD_LEVELSET_MESH_NODE_2D<T>(this,real_selection->grid_mapped,real_selection->index);}
}
//#####################################################################
// Function Clear_Highlight
//#####################################################################
template<class T> void OPENGL_NONMANIFOLD_LEVELSET_MESH_2D<T>::
Clear_Highlight()
{
    delete current_selection;current_selection=0;
}
//#####################################################################
// Print_Selection_Info
//#####################################################################
template<class T> void OPENGL_NONMANIFOLD_LEVELSET_MESH_2D<T>::
Print_Selection_Info(std::ostream& stream,OPENGL_SELECTION* selection) const
{
    if(current_selection && current_selection->type==OPENGL_SELECTION::GRID_NODE_3D){
        int grid_mapped=((OPENGL_SELECTION_NONMANIFOLD_LEVELSET_MESH_NODE_2D<T>*)current_selection)->grid_mapped==0?1:0;
        int mesh_index=((OPENGL_SELECTION_NONMANIFOLD_LEVELSET_MESH_NODE_2D<T>*)current_selection)->grid_mapped;
        VECTOR<int,2> index=((OPENGL_SELECTION_NONMANIFOLD_LEVELSET_MESH_NODE_2D<T>*)current_selection)->index;
        stream<<"Grid-Mapped: "<<grid_mapped<<std::endl<<"Selected node: ["<<index<<","<<mesh_index<<"] ("<<grid.Node(index)<<")"<<std::endl<<"Phi Value: "<<(*phi)(HINDEX(index,mesh_index))<<std::endl;}
}
//#####################################################################
// Function Bounding_Box
//#####################################################################
template<class T> RANGE<VECTOR<float,3> > OPENGL_SELECTION_NONMANIFOLD_LEVELSET_MESH_NODE_2D<T>::
Bounding_Box() const
{
    PHYSBAM_ASSERT(object);
    const GRID<TV> &grid=((OPENGL_NONMANIFOLD_LEVELSET_MESH_2D<T> *)object)->grid;
    RANGE<VECTOR<T,2> > box(grid.Node(index));
    return object->World_Space_Box(box);
}
//#####################################################################
template class OPENGL_NONMANIFOLD_LEVELSET_MESH_2D<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class OPENGL_NONMANIFOLD_LEVELSET_MESH_2D<double>;
#endif

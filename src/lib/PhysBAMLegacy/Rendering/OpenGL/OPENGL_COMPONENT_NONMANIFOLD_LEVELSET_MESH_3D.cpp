//#####################################################################
// Copyright (c) 2014, Mridul Aanjaneya.
// Distributed under the FreeBSD license (see license.txt)
//#####################################################################
#include <Rendering/OpenGL/OPENGL_COMPONENT_NONMANIFOLD_LEVELSET_MESH_3D.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T,class T2,class RW> OPENGL_COMPONENT_NONMANIFOLD_LEVELSET_MESH_3D<T,T2,RW>::
OPENGL_COMPONENT_NONMANIFOLD_LEVELSET_MESH_3D(GRID<TV>& grid,const std::string& filename_input,const std::string& flag_filename_input, const std::string& value_filename_input)
    : OPENGL_COMPONENT("Non-manifold Levelset Mesh 3D"),opengl_nonmanifold_levelset_mesh(grid),filename(filename_input),flag_filename(flag_filename_input),value_filename(value_filename_input),frame_loaded(-1),valid(false)
{
    is_animation = FILE_UTILITIES::Is_Animated(filename);
}
//#####################################################################
// Function Valid_Frame
//#####################################################################
template<class T,class T2,class RW> bool OPENGL_COMPONENT_NONMANIFOLD_LEVELSET_MESH_3D<T,T2,RW>::
Valid_Frame(int frame_input) const
{
    return FILE_UTILITIES::Frame_File_Exists(filename,frame_input);
}
//#####################################################################
// Function Set_Frame
//#####################################################################
template<class T,class T2,class RW> void OPENGL_COMPONENT_NONMANIFOLD_LEVELSET_MESH_3D<T,T2,RW>::
Set_Frame(int frame_input)
{
    OPENGL_COMPONENT::Set_Frame(frame_input);
    Reinitialize();
}
//#####################################################################
// Function Set_Draw
//#####################################################################
template<class T,class T2,class RW> void OPENGL_COMPONENT_NONMANIFOLD_LEVELSET_MESH_3D<T,T2,RW>::
Set_Draw(bool draw_input)
{
    OPENGL_COMPONENT::Set_Draw(draw_input);
    Reinitialize();
}
//#####################################################################
// Function Display
//#####################################################################
template<class T,class T2,class RW> void OPENGL_COMPONENT_NONMANIFOLD_LEVELSET_MESH_3D<T,T2,RW>::
Display(const int in_color) const
{
    if (valid && draw) opengl_nonmanifold_levelset_mesh.Display(in_color);
}
//#####################################################################
// Function Bounding_Box
//#####################################################################
template<class T,class T2,class RW> RANGE<VECTOR<float,3> > OPENGL_COMPONENT_NONMANIFOLD_LEVELSET_MESH_3D<T,T2,RW>::
Bounding_Box() const
{
    if (valid && draw) return opengl_nonmanifold_levelset_mesh.Bounding_Box();
    else return RANGE<VECTOR<float,3> >::Centered_Box();
}
//#####################################################################
// Function Print_Selection_Info
//#####################################################################
template<class T,class T2,class RW> void OPENGL_COMPONENT_NONMANIFOLD_LEVELSET_MESH_3D<T,T2,RW>::
Print_Selection_Info(std::ostream& output_stream,OPENGL_SELECTION* current_selection) const
{
    if(Is_Up_To_Date(frame)){
        output_stream<<component_name<<":"<<std::endl;
        opengl_nonmanifold_levelset_mesh.Print_Selection_Info(output_stream,current_selection);}
}
//#####################################################################
// Function Reinitialize
//#####################################################################
template<class T,class T2,class RW> void OPENGL_COMPONENT_NONMANIFOLD_LEVELSET_MESH_3D<T,T2,RW>::
Reinitialize()
{
    if(draw){
        if((is_animation && frame_loaded != frame) || (!is_animation && frame_loaded<0)){
            valid=false;
            std::string current_filename=FILE_UTILITIES::Get_Frame_Filename(filename,frame);
            std::string current_flag_filename=FILE_UTILITIES::Get_Frame_Filename(flag_filename,frame);
            std::string current_value_filename=FILE_UTILITIES::Get_Frame_Filename(value_filename,frame);

            if(FILE_UTILITIES::File_Exists(current_filename)){
                if(opengl_nonmanifold_levelset_mesh.frame!=frame){
                    opengl_nonmanifold_levelset_mesh.Set_Frame(frame);
                    opengl_nonmanifold_levelset_mesh.Initialize(current_filename,current_flag_filename,current_value_filename,frame);}}
            else return;

            frame_loaded=frame;
            valid=true;}}
}
//#####################################################################
template class OPENGL_COMPONENT_NONMANIFOLD_LEVELSET_MESH_3D<float,int,float>;
template class OPENGL_COMPONENT_NONMANIFOLD_LEVELSET_MESH_3D<float,bool,float>;
template class OPENGL_COMPONENT_NONMANIFOLD_LEVELSET_MESH_3D<float,float,float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class OPENGL_COMPONENT_NONMANIFOLD_LEVELSET_MESH_3D<double,int,double>;
template class OPENGL_COMPONENT_NONMANIFOLD_LEVELSET_MESH_3D<double,bool,double>;
template class OPENGL_COMPONENT_NONMANIFOLD_LEVELSET_MESH_3D<double,double,double>;
#endif

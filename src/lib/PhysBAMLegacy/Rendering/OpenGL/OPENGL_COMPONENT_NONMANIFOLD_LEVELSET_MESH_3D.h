//#####################################################################
// Copyright (c) 2014, Mridul Aanjaneya.
// Distributed under the FreeBSD license (see license.txt)
//#####################################################################
// Class OPENGL_COMPONENT_NONMANIFOLD_LEVELSET_MESH_3D
//##################################################################### 
#ifndef __OPENGL_COMPONENT_NONMANIFOLD_LEVELSET_MESH_3D__
#define __OPENGL_COMPONENT_NONMANIFOLD_LEVELSET_MESH_3D__

#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT.h>
#include <Rendering/OpenGL/OPENGL_NONMANIFOLD_LEVELSET_MESH_3D.h>

namespace PhysBAM{
template<class T,class T2=T,class RW=T>
class OPENGL_COMPONENT_NONMANIFOLD_LEVELSET_MESH_3D : public OPENGL_COMPONENT
{
    typedef VECTOR<T,3> TV;
public:
    OPENGL_COMPONENT_NONMANIFOLD_LEVELSET_MESH_3D(GRID<TV>& grid,const std::string& filename_input,const std::string& flag_filename_input, const std::string& value_filename_input);
    ~OPENGL_COMPONENT_NONMANIFOLD_LEVELSET_MESH_3D() {}

    OPENGL_NONMANIFOLD_LEVELSET_MESH_3D<T> opengl_nonmanifold_levelset_mesh;
private:
    std::string filename,flag_filename,value_filename;
    int frame_loaded;
    bool valid;

    virtual OPENGL_SELECTION *Get_Selection(GLuint *buffer, int buffer_size)
    {return opengl_nonmanifold_levelset_mesh.Get_Selection(buffer,buffer_size);}

    virtual void Highlight_Selection(OPENGL_SELECTION *selection) PHYSBAM_OVERRIDE
    {opengl_nonmanifold_levelset_mesh.Highlight_Selection(selection);}

    virtual void Clear_Highlight() PHYSBAM_OVERRIDE
    {opengl_nonmanifold_levelset_mesh.Clear_Highlight();}

    virtual void Set_Slice(OPENGL_SLICE *slice_input) PHYSBAM_OVERRIDE {slice=slice_input;opengl_nonmanifold_levelset_mesh.Set_Slice(slice_input);}

//#####################################################################
public:
    bool Valid_Frame(int frame_input) const PHYSBAM_OVERRIDE;
    bool Is_Up_To_Date(int frame) const PHYSBAM_OVERRIDE {return valid && frame_loaded==frame;}
    void Set_Frame(int frame_input) PHYSBAM_OVERRIDE;
    void Set_Draw(bool draw_input=true) PHYSBAM_OVERRIDE;
    void Display(const int in_color=1) const PHYSBAM_OVERRIDE;
    bool Use_Bounding_Box() const PHYSBAM_OVERRIDE {return draw&&valid;}
    virtual RANGE<VECTOR<float,3> > Bounding_Box() const PHYSBAM_OVERRIDE;
    void Print_Selection_Info(std::ostream& output_stream,OPENGL_SELECTION* current_selection) const PHYSBAM_OVERRIDE;

    void Increase_Mesh_Height(){opengl_nonmanifold_levelset_mesh.mesh_height+=opengl_nonmanifold_levelset_mesh.grid.dX.Max();};
    void Decrease_Mesh_Height(){opengl_nonmanifold_levelset_mesh.mesh_height-=opengl_nonmanifold_levelset_mesh.grid.dX.Max();};
    void Toggle_Collapse_Slice(){opengl_nonmanifold_levelset_mesh.Toggle_Collapse_Slice();}
    //void Increase_Node_Height(){opengl_nonmanifold_levelset_mesh.Increase_Node_Height();};
    //void Decrease_Node_Height(){opengl_nonmanifold_levelset_mesh.Decrease_Node_Height();};
    //void Reset_Node_Height(){opengl_nonmanifold_levelset_mesh.Reset_Node_Height();};

    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_NONMANIFOLD_LEVELSET_MESH_3D,Increase_Mesh_Height,"Raise Mesh Nodes");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_NONMANIFOLD_LEVELSET_MESH_3D,Decrease_Mesh_Height,"Lower Mesh Nodes");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_NONMANIFOLD_LEVELSET_MESH_3D,Toggle_Collapse_Slice,"Collapse Slice");
    //DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_NONMANIFOLD_LEVELSET_MESH_3D,Increase_Node_Height,"Raise Selected Node");
    //DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_NONMANIFOLD_LEVELSET_MESH_3D,Decrease_Node_Height,"Lower Selected Node");
    //DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_NONMANIFOLD_LEVELSET_MESH_3D,Reset_Node_Height,"Reset Height of all Nodes");


private:
    void Reinitialize();
//#####################################################################
};
}
#endif

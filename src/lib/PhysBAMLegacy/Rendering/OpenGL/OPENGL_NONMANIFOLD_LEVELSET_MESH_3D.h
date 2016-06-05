//#####################################################################
// Copyright (c) 2014, Mridul Aanjaneya.
// Distributed under the FreeBSD license (see license.txt)
//#####################################################################
// Class OPENGL_NONMANIFOLD_LEVELSET_MESH_3D
//##################################################################### 
#ifndef __OPENGL_NONMANIFOLD_LEVELSET_MESH_3D__
#define __OPENGL_NONMANIFOLD_LEVELSET_MESH_3D__

#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_CALLBACK.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_COLOR_RAMP.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_OBJECT.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_SELECTION.h>
#include <PhysBAM_Tools/Data_Structures/HASHTABLE.h>
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <Common_Tools/Arrays/HYBRID_ARRAY.h>
#include <Common_Tools/Arrays/HYBRID_ARRAY_ITERATOR.h>
#include <Nonmanifold_Implicit_Objects/NONMANIFOLD_LEVELSET_MESH.h>

namespace PhysBAM{

template<class T>
class OPENGL_NONMANIFOLD_LEVELSET_MESH_3D : public OPENGL_OBJECT
{
  public:
    typedef T RW;
    typedef VECTOR<T,3> TV;typedef VECTOR<int,3> T_INDEX;
    typedef PAIR< T_INDEX, int > HINDEX;
    typedef typename NONMANIFOLD_LEVELSET_MESH<T,3>::T_CELL T_CELL;

    GRID<TV>& grid;
    ARRAY<T_CELL>* cells;
    ARRAY<bool>* cell_is_grid;
    HYBRID_ARRAY<T,3>* phi;
    HYBRID_ARRAY<T,3>* node_heights;
    OPENGL_COLOR_RAMP<T>* cramp;
    HASHTABLE<T_INDEX, ARRAY< HINDEX > > node_overlaps;

    bool draw,collapse_slice;
    T mesh_height;
    int frame;
  private:
    STREAM_TYPE stream_type;
    OPENGL_SELECTION *current_selection;

  public:
    OPENGL_NONMANIFOLD_LEVELSET_MESH_3D(GRID<TV>& grid_input,const int frame_input=-1)
        :grid(grid_input),cells(0),cell_is_grid(0),phi(0),cramp(0),draw(true),collapse_slice(false),mesh_height(0),frame(frame_input),stream_type((RW())),current_selection(0)
    {
        mesh_height = grid.dX.Max();
    }

    ~OPENGL_NONMANIFOLD_LEVELSET_MESH_3D()
    {
        if(cells) delete cells;
        if(cell_is_grid) delete cell_is_grid;
        if(phi) delete phi;
        if(cramp) delete cramp;
    }

    void Toggle_Collapse_Slice() {collapse_slice=!collapse_slice;}

//##################################################################### 
    void Display(const int in_color=1) const PHYSBAM_OVERRIDE;
    virtual void Set_Frame(int frame_input);
    virtual RANGE<VECTOR<float,3> > Bounding_Box() const PHYSBAM_OVERRIDE;
    virtual OPENGL_SELECTION *Get_Selection(GLuint *buffer, int buffer_size);
    void Highlight_Selection(OPENGL_SELECTION *selection) PHYSBAM_OVERRIDE;
    void Clear_Highlight() PHYSBAM_OVERRIDE;
    void Initialize(const std::string& filename,const std::string& flag_filename, const std::string& value_filename,  const int frame);
    void Print_Selection_Info(std::ostream& stream,OPENGL_SELECTION* selection) const PHYSBAM_OVERRIDE;
    //void Increase_Node_Height();
    //void Decrease_Node_Height();
    void Reset_Node_Height();

//##################################################################### 
};

template<class T>
class OPENGL_SELECTION_NONMANIFOLD_LEVELSET_MESH_NODE_3D : public OPENGL_SELECTION
{
  private:
    typedef VECTOR<T,3> TV;typedef VECTOR<int,3> T_INDEX;
    typedef PAIR<T_INDEX,int> HINDEX;

  public:
    int grid_mapped;
    T_INDEX index;
    HINDEX hindex;

    OPENGL_SELECTION_NONMANIFOLD_LEVELSET_MESH_NODE_3D(OPENGL_OBJECT *object,const int grid_mapped_input=0,const T_INDEX& index_input=T_INDEX())
     : OPENGL_SELECTION(OPENGL_SELECTION::GRID_NODE_3D,object),grid_mapped(grid_mapped_input),index(index_input),hindex(index_input,grid_mapped_input) {}

    RANGE<VECTOR<float,3> > Bounding_Box() const PHYSBAM_OVERRIDE;
};
//##################################################################### 
}
#endif

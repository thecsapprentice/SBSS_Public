#ifndef __REGION_GENERATOR_H__
#define __REGION_GENERATOR_H__

#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Data_Structures/UNION_FIND.h>
#include <PhysBAM_Tools/Data_Structures/HASHTABLE.h>
#include <PhysBAM_Tools/Data_Structures/HASHTABLE_ITERATOR.h>


namespace PhysBAM
{

    enum AXIS { A_X, A_mX, A_Y, A_mY, A_Z, A_mZ };


template<int d>
struct T_FACETYPE
{
};

template<>
struct T_FACETYPE<2>
{
    static const int vertices_per_face=2;
    int indices[2];
};


template<>
struct T_FACETYPE<3>
{
    static const int vertices_per_face=4;
    int indices[4];
};

template<class T, int d>
struct T_CELLTYPE
{
};

template<class T>
struct T_CELLTYPE<T, 2>
{
    static const int vertices_per_cell=4;
    int indices[4];
    void Get_X_Face(T_FACETYPE<2>& face) const { face.indices[0] = indices[3]; face.indices[1] = indices[2]; }
    void Get_mX_Face(T_FACETYPE<2>& face) const { face.indices[0] = indices[1]; face.indices[1] = indices[0]; }
    void Get_Y_Face(T_FACETYPE<2>& face) const { face.indices[0] = indices[1]; face.indices[1] = indices[3]; }
    void Get_mY_Face(T_FACETYPE<2>& face) const { face.indices[0] = indices[0]; face.indices[1] = indices[2]; }

    void GetAxisFace(T_FACETYPE<2>& face, AXIS axis) const { 
        switch(axis){
        case A_X:
            Get_X_Face(face);
            break;
        case A_mX:
            Get_mX_Face(face);
            break;
        case A_Y:
            Get_Y_Face(face);
            break;
        case A_mY:
            Get_mY_Face(face);
            break;
        default:
            PHYSBAM_FATAL_ERROR("Axis not acceptable for this dimension.");
        }
    }
 };



template<class T>
struct T_CELLTYPE<T, 3>
{
    static const int vertices_per_cell=8;
    int indices[8];

    void Get_X_Face(T_FACETYPE<3>& face) const  { face.indices[0] = indices[4]; face.indices[1] = indices[5];
                                           face.indices[2] = indices[6]; face.indices[3] = indices[7];}
    void Get_mX_Face(T_FACETYPE<3>& face) const { face.indices[0] = indices[0]; face.indices[1] = indices[1];
                                           face.indices[2] = indices[2]; face.indices[3] = indices[3];}

    void Get_Y_Face(T_FACETYPE<3>& face) const  { face.indices[0] = indices[2]; face.indices[1] = indices[3];
                                           face.indices[2] = indices[6]; face.indices[3] = indices[7];}
    void Get_mY_Face(T_FACETYPE<3>& face) const { face.indices[0] = indices[0]; face.indices[1] = indices[1];
                                           face.indices[2] = indices[4]; face.indices[3] = indices[5];}

    void Get_Z_Face(T_FACETYPE<3>& face) const  { face.indices[0] = indices[1]; face.indices[1] = indices[3];
                                           face.indices[2] = indices[5]; face.indices[3] = indices[7];}
    void Get_mZ_Face(T_FACETYPE<3>& face) const { face.indices[0] = indices[0]; face.indices[1] = indices[2];
                                           face.indices[2] = indices[4]; face.indices[3] = indices[6];}

    void GetAxisFace(T_FACETYPE<3>& face, AXIS axis) const { 
        switch(axis){
        case A_X:
            Get_X_Face(face);
            break;
        case A_mX:
            Get_mX_Face(face);
            break;
        case A_Y:
            Get_Y_Face(face);
            break;
        case A_mY:
            Get_mY_Face(face);
            break;
        case A_Z:
            Get_Z_Face(face);
            break;
        case A_mZ:
            Get_mZ_Face(face);
            break;
        default:
            PHYSBAM_FATAL_ERROR("Axis not acceptable for this dimension.");
        }
    }

 };


std::ostream& operator<<(std::ostream& out, const T_CELLTYPE<float, 2>& cell);
std::ostream& operator<<(std::ostream& out, const T_CELLTYPE<double, 2>& cell);
std::ostream& operator<<(std::ostream& out, const T_CELLTYPE<float,3>& cell);
std::ostream& operator<<(std::ostream& out, const T_CELLTYPE<double,3>& cell);

template< class T, int d >
class REGIONS
{
 public:
    typedef VECTOR<T,d> TV;
    typedef T_CELLTYPE<T,d> T_CELL;
    typedef ARRAY<T_CELL> T_REGION;
    typedef ARRAY<T_REGION> T_REGION_ARRAY; 
    typedef VECTOR<int,d> T_INDEX;
    typedef ARRAY<T_INDEX> T_GRID_REGION;
    typedef ARRAY<T_GRID_REGION> T_GRID_REGION_ARRAY;

    T_REGION_ARRAY regions;
    T_GRID_REGION_ARRAY grid_regions;
    int max_index;
    ARRAY<TV> vertices;
};

template<class T, int d>
class REGION_GENERATOR
{
 protected:
    static const int vertices_per_face=(d==2)?2:4;
    UNION_FIND<int>  region_sets;
    UNION_FIND<int>  equivalence_sets;
    REGIONS<T,d> regions;

 public:
    typedef VECTOR<int, d> T_INDEX;
    typedef VECTOR<T, d> TV;
    typedef GRID<TV> T_GRID;
   
    typedef T_CELLTYPE<T,d> T_CELL;
    typedef T_FACETYPE<d> T_FACE;

    typedef ARRAY<T_CELL> T_REGION;
    typedef ARRAY<T_CELL> T_CELL_ARRAY;


    REGION_GENERATOR(const T_GRID& grid_input);
    
    virtual void Generate();
    virtual const REGIONS<T,d>* GetRegionData() const;
    int DuplicatesAtCoarseIndex( const T_INDEX& cell_index ) const; 
    bool IsIncident( const T_INDEX& cell_index, const int subcell, const T_INDEX& test_index );
    bool IsMeshMappable( const T_INDEX& cell_index) const;
    VECTOR<int,T_CELLTYPE<T, d>::vertices_per_cell> GetCellVertices( const T_INDEX& cell_index, int duplicate) const;
    virtual VECTOR<bool,T_CELLTYPE<T, d>::vertices_per_cell> CornersInsideOutside( const T_INDEX& cell_index, const int subcell ) const = 0;

    //protected:
 public :

    const T_GRID grid;

    ARRAY<T_INDEX> root_grid_cell_mapping;
    HASHTABLE< T_INDEX, int > root_index_to_root_id; 

    ARRAY<T_INDEX> root_grid_node_mapping;
    ARRAY<ARRAY<int>,T_INDEX> node_root_grid_mapping; // For each grid node, list of mesh vertices

    T_CELL_ARRAY root_cells;
    
    ARRAY<int> root_sub_mapping; // Maps sub_cell indexes to root_cell indexes, i.e. 
                                 //    root_sub_mapping( sub_id ) = root_id;

    ARRAY< ARRAY<int> > rootnode_to_subnode; // Maps root node indexes to subcell node indexes, i.e. 
                                             //    rootnode_to_subnode( root_node_id ) = list(subnodes);

    ARRAY< ARRAY< int > > root_cell_to_subcell; // For each root cell, return list of subcells

    ARRAY< ARRAY< int > > subcell_incident_list; // For each subcell, provide a list of incident, touching, 
                                        // cells
    
    ARRAY< ARRAY< int > > node_to_subcell; // For each mesh vertex, list cells it belongs to.
   

    T_CELL_ARRAY sub_cells;

    

    virtual void ResolveCuts() = 0;
    virtual bool IsMaterialContinous(const int root_cellA, const int sub_cellA, const AXIS faceA,
                                     const int root_cellB, const int sub_cellB, const AXIS faceB) const = 0;

    void Generate_EquivalenceClass_Range(int cell_start, int cell_end, UNION_FIND<int> & ec, UNION_FIND<int>& rc );

 private:

    void CreateRootCells();
    void BuildIncidentLists();

};


}


#endif

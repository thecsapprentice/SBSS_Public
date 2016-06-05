//#####################################################################
// Copyright 2014, Raj Setaluri, Nathan Mitchell.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class NONMANIFOLD_DUALCONTOUR_2D
//#####################################################################
#ifndef __NONMANIFOLD_DUALCONTOUR_2D__
#define __NONMANIFOLD_DUALCONTOUR_2D__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <Nonmanifold_Implicit_Objects/NONMANIFOLD_LEVELSET_2D.h>
namespace PhysBAM{

template<class T> class SEGMENTED_CURVE_2D;
template<class T> class TRIANGULATED_AREA;

template<class T>
class NONMANIFOLD_DUALCONTOUR_2D:public NONCOPYABLE
{
    typedef VECTOR<T,2> TV;typedef VECTOR<int,2> TV_INT;
    typedef typename NONMANIFOLD_LEVELSET_2D<T>::T_ARRAY_SCALAR::template REBIND<int>::TYPE T_ARRAY_INT;
    typedef typename NONMANIFOLD_LEVELSET_2D<T>::T_ARRAY_SCALAR::INDEX T_INDEX;

private:
    NONMANIFOLD_LEVELSET_2D<T>& levelset;
    T contour_value;
    HASHTABLE<PAIR<TV_INT,int>,int> node_to_int;
    ARRAY<PAIR<TV_INT,int> > int_to_node;
    ARRAY<VECTOR<int,2> > topology;
    //HASHTABLE<PAIR<int,int>,int> vertices;
    T_ARRAY_INT vertices;
    
    ARRAY<TV> geometry;
    ARRAY<TV> normals;
    //GRID<TV> grid;
    bool is_distance_field;
public:

    NONMANIFOLD_DUALCONTOUR_2D(NONMANIFOLD_LEVELSET_2D<T>& levelset_input,const T contour_value_input=0,const bool is_distance_field_input=true)
        :levelset(levelset_input),contour_value(contour_value_input),is_distance_field(is_distance_field_input)
    {}

    void Dualcontour()
    {Generate_Topology();Generate_Vertices();Ensure_Vertices_In_Correct_Cells();}
    
    static SEGMENTED_CURVE_2D<T>* Create_Segmented_Curve_From_Levelset(NONMANIFOLD_LEVELSET_2D<T>& levelset,const T contour_value=0,const bool is_distance_field_input=true)
    {NONMANIFOLD_DUALCONTOUR_2D<T> dualcontour(levelset,contour_value,is_distance_field_input);dualcontour.Dualcontour();return dualcontour.Get_Segmented_Curve();} 

    static TRIANGULATED_AREA<T>* Create_Triangulated_Area_From_Levelset(NONMANIFOLD_LEVELSET_2D<T>& levelset,const int sign=-1,const T contour_value=0,const bool is_distance_field_input=true)
    {NONMANIFOLD_DUALCONTOUR_2D<T> dualcontour(levelset,contour_value,is_distance_field_input);dualcontour.Dualcontour();return dualcontour.Get_Triangulated_Area(sign);} 

//#####################################################################
    void Generate_Topology();
    void Generate_Vertices();
    void Ensure_Vertices_In_Correct_Cells();
    SEGMENTED_CURVE_2D<T>* Get_Segmented_Curve();
    TRIANGULATED_AREA<T>* Get_Triangulated_Area(const int sign=-1);
//#####################################################################
};
}
#endif

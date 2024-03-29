//#####################################################################
// Copyright 2005, Eran Guendelman.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class NONMANIFOLD_FAST_MARCHING_METHOD_UNIFORM  
//#####################################################################
#ifndef __NONMANIFOLD_FAST_MARCHING_METHOD_UNIFORM__
#define __NONMANIFOLD_FAST_MARCHING_METHOD_UNIFORM__

#include <PhysBAM_Tools/Arrays/ARRAYS_FORWARD.h>
#include <PhysBAM_Tools/Parallel_Computation/THREAD_QUEUE.h>
#include <Nonmanifold_Implicit_Objects/NONMANIFOLD_LEVELSET_POLICY_UNIFORM.h>
#include <Nonmanifold_Implicit_Objects/NONMANIFOLD_LEVELSET_UNIFORM.h>
namespace PhysBAM{


    template<class T, int d>
    class NONMANIFOLD_FAST_MARCHING_METHOD_UNIFORM
    {

        typedef typename NONMANIFOLD_LEVELSET_POLICY<T,d>::NONMANIFOLD_LEVELSET T_NONMANIFOLD_LEVELSET;
        typedef typename NONMANIFOLD_LEVELSET_MESH<T,d>::AXIS_ALIGNED_NEIGHBORS AA_NEIGHBORS;
        typedef HYBRID_ARRAY<T,d> T_HYBRID_ARRAY;
        typedef PAIR< VECTOR<int,d>, int> HINDEX;
        typedef VECTOR<int,d> T_INDEX;

    public:
        const T_NONMANIFOLD_LEVELSET& nonmanifold_levelset;
        const T stopping_distance;
        const NONMANIFOLD_LEVELSET_MESH<T,d>& cell_mesh;
        const RANGE<T_INDEX> domain_range;
        const int flat_range;

        NONMANIFOLD_FAST_MARCHING_METHOD_UNIFORM(const T_NONMANIFOLD_LEVELSET& nonmanifold_levelset);
        ~NONMANIFOLD_FAST_MARCHING_METHOD_UNIFORM();

        void Fast_Marching_Method(HYBRID_ARRAY<T,d>& phi, HYBRID_ARRAY<bool,d>& done);
    private:
        void Update_Or_Add_Neighbor(HYBRID_ARRAY<T,d>& phi, HYBRID_ARRAY<bool,d>& done,
                                    HYBRID_ARRAY<int,d>& close_k,ARRAY<HINDEX>& heap, int& heap_length,
                                    const HINDEX& neighbor);
        void Initialize_Interface(HYBRID_ARRAY<T,d>& phi, HYBRID_ARRAY<bool,d>& done,
                                  HYBRID_ARRAY<int,d>& close_k,ARRAY<HINDEX>& heap, int& heap_length);
        void Update_Close_Point(HYBRID_ARRAY<T,d>& phi,HYBRID_ARRAY<bool,d>& done, const HINDEX& hindex);
        void Add_To_Initial(HYBRID_ARRAY<bool,d>& done,HYBRID_ARRAY<int,d>& close_k,const HINDEX& hindex);

    };





#if 0
template<class T_GRID>
class NONMANIFOLD_FAST_MARCHING_METHOD_UNIFORM:public NONCOPYABLE
{
    typedef typename T_GRID::VECTOR_T TV;typedef typename TV::SCALAR T;
    typedef typename T_GRID::VECTOR_INT TV_INT;typedef typename NONMANIFOLD_LEVELSET_POLICY<T_GRID>::NONMANIFOLD_LEVELSET T_NONMANIFOLD_LEVELSET;
    typedef typename T_GRID::CELL_ITERATOR CELL_ITERATOR;typedef typename T_GRID::FACE_ITERATOR FACE_ITERATOR;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_SCALAR T_ARRAYS_SCALAR;typedef typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_SCALAR::template REBIND<bool>::TYPE T_ARRAYS_BOOL;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_SCALAR::template REBIND<int>::TYPE T_ARRAYS_INT;
public:
    const T_NONMANIFOLD_LEVELSET& nonmanifold_levelset;
private:
    T_GRID cell_grid;
    TV_INT dimension_start,dimension_end;
    int ghost_cells;
public:
    THREAD_QUEUE* thread_queue;

    NONMANIFOLD_FAST_MARCHING_METHOD_UNIFORM(const T_NONMANIFOLD_LEVELSET& nonmanifold_levelset,const int ghost_cells,THREAD_QUEUE* thread_queue_input=0);
    ~NONMANIFOLD_FAST_MARCHING_METHOD_UNIFORM();

    bool Neighbor_Visible(const int neighbor_number,const TV_INT& current_index) // neighbor_number between 1 and 3 -- right, top, back
    {return !nonmanifold_levelset.collision_body_list->cell_neighbors_visible.Valid_Index(current_index) || nonmanifold_levelset.collision_body_list->cell_neighbors_visible(current_index)(neighbor_number);}

//#####################################################################
    //void Fast_Marching_Method_Threaded(RANGE<TV_INT>& domain,T_ARRAYS_SCALAR& phi_ghost,const T stopping_distance,const ARRAY<TV_INT>* seed_indices,const bool add_seed_indices_for_ghost_cells=false);
    //void Fast_Marching_Method(T_ARRAYS_SCALAR& phi_ghost,const T stopping_distance=0,const ARRAY<TV_INT>* seed_indices=0,const bool add_seed_indices_for_ghost_cells=false);
    void Fast_Marching_Method(T_ARRAYS_SCALAR& phi_ghost,T_ARRAYS_BOOL& seed_indices,const T stopping_distance=0,const bool add_seed_indices_for_ghost_cells=false);
    //void Initialize_Interface_Threaded(RANGE<TV_INT>& domain,T_ARRAYS_SCALAR& phi_ghost,T_ARRAYS_SCALAR& phi_new,T_ARRAYS_BOOL& done);
private:
    void Update_Or_Add_Neighbor(T_ARRAYS_SCALAR& phi_ghost,T_ARRAYS_BOOL& done,T_ARRAYS_INT& close_k,ARRAY<TV_INT>& heap,int& heap_length,const TV_INT& neighbor);
    //void Initialize_Interface(RANGE<TV_INT>& domain,T_ARRAYS_SCALAR& phi_ghost,T_ARRAYS_BOOL& done,T_ARRAYS_INT& close_k,ARRAY<TV_INT>& heap,int& heap_length,const ARRAY<TV_INT>* seed_indices,const bool add_seed_indices_for_ghost_cells=false);
    //void Initialize_Interface(T_ARRAYS_SCALAR& phi_ghost,T_ARRAYS_BOOL& done,T_ARRAYS_INT& close_k,ARRAY<TV_INT>& heap,int& heap_length,const ARRAY<TV_INT>* seed_indices=0, const bool add_seed_indices_for_ghost_cells=false);
    void Initialize_Interface(T_ARRAYS_SCALAR& phi_ghost,T_ARRAYS_BOOL& done,T_ARRAYS_INT& close_k,ARRAY<TV_INT>& heap,int& heap_length,const bool add_seed_indices_for_ghost_cells=false);
    void Update_Close_Point(T_ARRAYS_SCALAR& phi_ghost,const T_ARRAYS_BOOL& done,const TV_INT& index);
    void Add_To_Initial(T_ARRAYS_BOOL& done,T_ARRAYS_INT& close_k,const TV_INT& index);
//#####################################################################
};
#endif

}
#endif

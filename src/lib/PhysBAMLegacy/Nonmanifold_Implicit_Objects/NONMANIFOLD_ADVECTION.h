//#####################################################################
// Copyright 2014, Nathan Mitchell
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class NONMANIFOLD_ADVECTION  
//#####################################################################
#ifndef __NONMANIFOLD_ADVECTION__
#define __NONMANIFOLD_ADVECTION__

#include <PhysBAM_Tools/Arrays/ARRAYS_FORWARD.h>
#include <PhysBAM_Tools/Parallel_Computation/THREAD_QUEUE.h>
#include <Nonmanifold_Implicit_Objects/NONMANIFOLD_LEVELSET_POLICY_UNIFORM.h>
#include <Nonmanifold_Implicit_Objects/NONMANIFOLD_LEVELSET_UNIFORM.h>
namespace PhysBAM{


    template<class T, int d>
    class NONMANIFOLD_ADVECTION
    {

        typedef typename NONMANIFOLD_LEVELSET_POLICY<T,d>::NONMANIFOLD_LEVELSET T_NONMANIFOLD_LEVELSET;
        typedef typename NONMANIFOLD_LEVELSET_MESH<T,d>::AXIS_ALIGNED_NEIGHBORS AA_NEIGHBORS;
        typedef HYBRID_ARRAY<T,d> T_HYBRID_ARRAY;
        typedef PAIR< VECTOR<int,d>, int> HINDEX;
        typedef VECTOR<int,d> T_INDEX;

    public:
        T_NONMANIFOLD_LEVELSET& nonmanifold_levelset;
        const T stopping_distance;
        const NONMANIFOLD_LEVELSET_MESH<T,d>& cell_mesh;
        const RANGE<T_INDEX> domain_range;
        const int flat_range;

        NONMANIFOLD_ADVECTION(T_NONMANIFOLD_LEVELSET& nonmanifold_levelset);
        ~NONMANIFOLD_ADVECTION();

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


#endif

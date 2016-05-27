//#####################################################################
// Copyright 2011-2013, Nathan Mitchell, Taylor Patterson, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class CG_POLICY
//#####################################################################
#ifndef __CG_POLICY__
#define __CG_POLICY__

namespace PhysBAM{

template<class T_ELASTICITY> struct CG_POLICY;

template<class T, int d> class NONLINEAR_ELASTICITY;
template<class T, int d, bool enable_constraints, bool enable_muscles> class SKINNING_NONLINEAR_ELASTICITY;
template<class T_ELASTICITY> class HYBRID_NONLINEAR_ELASTICITY;

template<class T,int d> class NONLINEAR_ELASTICITY_STATE;
template<class T_ELASTICITY> class HYBRID_NONLINEAR_ELASTICITY_STATE;

template<class ARRAY, int d, int stride> class HYBRID_BINDER;
template<class ARRAY, int d, int stride> class BINDER;
template<class T, class ID, unsigned int align> class ALIGNED_ARRAY;

template<class T,int d,bool enable_constraints,bool enable_muscles>
struct CG_POLICY<HYBRID_NONLINEAR_ELASTICITY<SKINNING_NONLINEAR_ELASTICITY<T,d,enable_constraints,enable_muscles> > >
{
    typedef HYBRID_NONLINEAR_ELASTICITY<SKINNING_NONLINEAR_ELASTICITY<T,d,enable_constraints,enable_muscles> > ELASTICITY;
    typedef HYBRID_NONLINEAR_ELASTICITY_STATE<SKINNING_NONLINEAR_ELASTICITY<T,d,enable_constraints,enable_muscles> > STATE;
    typedef HYBRID_BINDER<ALIGNED_ARRAY<T, int, 64>, d, 16> STATE_BINDER;
    static const bool SUPPORTS_MESH = true;
};

template<class T,int d>
struct CG_POLICY<HYBRID_NONLINEAR_ELASTICITY<NONLINEAR_ELASTICITY<T,d> > >
{
    typedef HYBRID_NONLINEAR_ELASTICITY<NONLINEAR_ELASTICITY<T,d> > ELASTICITY;
    typedef HYBRID_NONLINEAR_ELASTICITY_STATE<NONLINEAR_ELASTICITY<T,d> > STATE;
    typedef HYBRID_BINDER<ALIGNED_ARRAY<T, int, 64>, d, 16> STATE_BINDER;
    static const bool SUPPORTS_MESH = true;
};

template<class T,int d,bool enable_constraints,bool enable_muscles>
struct CG_POLICY<SKINNING_NONLINEAR_ELASTICITY<T,d,enable_constraints,enable_muscles> >
{
    typedef NONLINEAR_ELASTICITY<T,d> ELASTICITY;
    typedef NONLINEAR_ELASTICITY_STATE<T,d> STATE;
    typedef BINDER<ALIGNED_ARRAY<T, int, 64>, d, 16> STATE_BINDER;
    static const bool SUPPORTS_MESH = false;
};

template<class T,int d>
struct CG_POLICY<NONLINEAR_ELASTICITY<T,d> >
{
    typedef NONLINEAR_ELASTICITY<T,d> ELASTICITY;
    typedef NONLINEAR_ELASTICITY_STATE<T,d> STATE;
    typedef BINDER<ALIGNED_ARRAY<T, int, 64>, d, 16> STATE_BINDER;
    static const bool SUPPORTS_MESH = false;
};


}

#endif

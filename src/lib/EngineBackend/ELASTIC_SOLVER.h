#ifndef __ELASTIC_SOLVER_H__
#define __ELASTIC_SOLVER_H__

#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <Common/STENCIL.h>
#include <Common/STENCIL_ITERATOR.h>

namespace PhysBAM{

    template<class BINDER, class STATE, class DISCRETIZATION>
    class ELASTIC_SOLVER{
    public:
        typedef typename DISCRETIZATION::SCALAR SCALAR;
        typedef typename DISCRETIZATION::SCALAR T;
        enum {d=DISCRETIZATION::dim};

        typedef VECTOR<T,d> TV;
        typedef VECTOR<int,d> T_INDEX;
        
        typedef STENCIL<T,d> T_STENCIL;
        typedef STENCIL_ITERATOR<T,d> T_STENCIL_ITERATOR;
        typedef STENCIL_ITERATOR<const T,d> T_CONST_STENCIL_ITERATOR;

    private:
        
        DISCRETIZATION* engine;

        // Internal State
        BINDER *u_bind, *f_bind;
        STATE *u_internal, *f_internal;

        // Solver Scratch Space
        BINDER *x_b, *b_b, *q_b, *s_b, *r_b, *k_b, *z_b;        
        STATE *x, *b, *q, *s, *r, *k, *z;

        // Debug Variables
        BINDER *b_debug_b,*b2_debug_b;
        STATE *b_debug, *b2_debug;

        // Preconditioner
        BINDER *D_b;
        BINDER *D_primal_b;
        BINDER *D_dual_b;
        BINDER *omega_b;
        STATE *D;
        STATE *D_primal;
        STATE *D_dual;
        STATE *omega;



        bool solver_initialized;
        bool owns_engine;

        void UpdateOmega();
        void Extract_D(STATE& d_state, int subdomain, const T_INDEX s_origin, int radius, bool build_complement);
        
    public:

        ELASTIC_SOLVER( DISCRETIZATION* engine, bool owns_engine=false );
        ~ELASTIC_SOLVER();

        void Initialize(int partitions);

        float Exact_Solve( int krylov_iterations, int newton_iterations, T krylov_tolerance, T newton_tolerance, bool no_cut_cells);       

        const STATE& U_const() const;
        const STATE& F_const() const;
        STATE& U();
        STATE& F();       
    };
}

#endif

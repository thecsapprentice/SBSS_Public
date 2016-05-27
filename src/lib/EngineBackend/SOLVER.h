#ifndef __ELASTIC_SOLVER_H__
#define __ELASTIC_SOLVER_H__

namespace PhysBAM{

    template<class BINDER, class STATE, class DISCRETIZATION>
    class ELASTIC_SOLVER{
    private:
        
        DISCRETIZATION* engine;

        BINDER *u_bind, *f_bind;
        BINDER *b_debug_b,*b2_debug_b;
        BINDER *x_b, *b_b, *q_b, *s_b, *r_b, *k_b, *z_b;
        
        STATE *u_internal, *f_internal;
        STATE *b_debug, *b2_debug;
        STATE *x, *b, *q, *s, *r, *k, *z;

        // For preconditioner
        STATE_BINDER *D_b;
        STATE *D;
        STATE_BINDER *D_primal_b;
        STATE *D_primal;
        STATE_BINDER *D_dual_b;
        STATE *D_dual;
        STATE_BINDER *omega_b;
        STATE *omega;

        bool solver_initialized;

        void UpdateOmega();
        void Extract_D(int subdomain, const T_INDEX s_origin, int radius, bool build_complement);
        
    public:

        ELASTIC_SOLVER( DISCRETIZATION* engine, bool owns_engine=false );
        ~ELASTIC_SOLVER();

        void Initialize();

        float Exact_Solve( int krylov_iterations, int newton_iterations, T krylov_tolerance, T newton_tolerance, bool no_cut_cells);       

        const STATE& U_const() const;
        const STATE& F_const() const;
        STATE& U();
        STATE& F();       
    };
}

#endif

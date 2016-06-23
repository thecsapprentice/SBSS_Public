#include <PhysBAM_Tools/Krylov_Solvers/SYMMQMR.h>
#include <PhysBAM_Tools/Krylov_Solvers/CONJUGATE_GRADIENT.h>
#include <PhysBAM_Tools/Data_Structures/HASHTABLE.h>

#include <Common/ALIGNED_ARRAY.h>
#include <Common_Tools/Math_Tools/RANGE_ITERATOR.h>

#include "ELASTIC_SOLVER.h"

#include "CG_SYSTEM.h"
#include "CG_VECTOR.h"

#include "NONLINEAR_ELASTICITY.h"
#include "SKINNING_NONLINEAR_ELASTICITY.h"
#include "HYBRID_NONLINEAR_ELASTICITY.h"
#include "Write_Output.h"


using namespace PhysBAM;

template<class BINDER, class STATE, class DISCRETIZATION>
ELASTIC_SOLVER<BINDER,STATE,DISCRETIZATION>::
ELASTIC_SOLVER( DISCRETIZATION* engine_in, bool owns_engine_in ) :
    engine(engine_in), owns_engine( owns_engine_in ), solver_initialized(false)
{
    PHYSBAM_ASSERT( engine != NULL );

    u_bind = new BINDER();
    f_bind = new BINDER();
    //b_debug_b = new BINDER();
    //b2_debug_b = new BINDER();
    x_b = new BINDER();
    b_b = new BINDER();
    q_b = new BINDER();
    s_b = new BINDER();
    //r_b = new BINDER();
    //k_b = new BINDER();
    //z_b = new BINDER();
    D_b = new BINDER();
    //D_primal_b = new BINDER();
    //D_dual_b = new BINDER();
    //omega_b = new BINDER();
    
    u_internal = new STATE();
    f_internal = new STATE();
    //b_debug = new STATE();
    //b2_debug = new STATE();
    x = new STATE();
    b = new STATE();
    q = new STATE();
    s = new STATE();
    //r = new STATE();
    //k = new STATE();
    //z = new STATE();
    D = new STATE();
    //D_primal = new STATE();
    //D_dual = new STATE();
    //omega = new STATE();
}


template<class BINDER, class STATE, class DISCRETIZATION>
ELASTIC_SOLVER<BINDER,STATE,DISCRETIZATION>::~ELASTIC_SOLVER(){

    delete u_bind; 
    delete f_bind;
    //delete b_debug_b;
    //delete b2_debug_b;
    delete x_b;
    delete b_b;
    delete q_b;
    delete s_b;
    //delete r_b;
    //delete k_b;
    //delete z_b;
    delete D_b;
    //delete D_primal_b;
    //delete D_dual_b;
    //delete omega_b;
    delete u_internal;
    delete f_internal;
    //delete b_debug;
    //delete b2_debug;
    delete x;
    delete b;
    delete q;
    delete s;
    //delete r;
    //delete k;
    //delete z;
    delete D;
    //delete D_primal;
    //delete D_dual;
    //delete omega;

    if( owns_engine ) 
        if( engine != NULL )
            delete engine;
}


template<class BINDER, class STATE, class DISCRETIZATION> void
ELASTIC_SOLVER<BINDER,STATE,DISCRETIZATION>::UpdateOmega()
{
    const float CONSTRAINT_OMEGA = 0.2;
    
    VECTOR< const ARRAY<CONSTRAINT_SEGMENT<T,d> >*, 3 > constraint_v(&engine->dynamic_point_constraints,
                                                                     &engine->static_point_constraints,
                                                                     &engine->collision_constraints);
    LOG::SCOPE scope( "Adding Constraint Contributions..." );
    for(int ctype = 1; ctype <=3 ; ctype ++){            
        const ARRAY<CONSTRAINT_SEGMENT<T,d> >& constraints = *constraint_v(ctype);
        for(int c=1;c<=constraints.m;c++){
            
            if( constraints(c).isDualPoint() ){
                for(int e=1; e<=2; e++){
                    const CONSTRAINT_NODE<T,d>& endpoint = constraints(c).endpoints[e-1];
                    switch(endpoint.type){
                    case CONSTRAINT_NODE<T,d>::GRID_FIXED:
                        {
                            const T_INDEX& grid_index = endpoint.grid_index();
                            for( RANGE_ITERATOR<d> iter( RANGE<T_INDEX>(grid_index, grid_index+1) ); iter.Valid(); iter.Next() )
                                for( int v = 1; v<=d; v++)
                                    omega->x(v)(iter.Index()) = CONSTRAINT_OMEGA;                        
                        }                                        
                        break;
                    case CONSTRAINT_NODE<T,d>::MESH_FIXED:
                        {
                            const T_INDEX& grid_index = engine->cell_indices_mesh( endpoint.mesh_index() );
                            int flat_index = 1;
                            for( RANGE_ITERATOR<d> iter( RANGE<T_INDEX>(grid_index, grid_index+1) ); iter.Valid(); iter.Next(), flat_index++ ){
                                int mesh_node = engine->cells_mesh(endpoint.mesh_index())(flat_index);
                                if(mesh_node == 0)
                                    for( int v = 1; v<=d; v++)
                                        omega->x(v)(iter.Index()) = CONSTRAINT_OMEGA;    
                                else
                                    for( int v = 1; v<=d; v++)
                                        omega->x_mesh(v)(mesh_node) = CONSTRAINT_OMEGA;                                    
                            }
                        }                                        
                        break;
                    default:
                        PHYSBAM_FATAL_ERROR();
                    }
                }
            }
        }
    }  
}


template<class BINDER, class STATE, class DISCRETIZATION> void
ELASTIC_SOLVER<BINDER,STATE,DISCRETIZATION>::Extract_D(STATE& d_state, int subdomain, const T_INDEX s_origin, int radius, bool build_complement)
{
    CG_SYSTEM<DISCRETIZATION> krylov_system(*engine);
    CG_VECTOR<CG_POLICY<DISCRETIZATION> > krylov_b_debug(*engine,*b_debug);
    CG_VECTOR<CG_POLICY<DISCRETIZATION> > krylov_b2(*engine,*b2_debug);
    CG_VECTOR<CG_POLICY<DISCRETIZATION> > krylov_d(*engine,d_state);

    LOG::SCOPE scope("Extracting Diagonal");
    

    T norm = 0;

    HASHTABLE<int, T_INDEX> node_to_index;
    for(int cell=1; cell <= engine->Number_Of_Mesh_Cells(); cell++){
        const T_INDEX& index = engine->cell_indices_mesh(cell);
        int vertex=1;
        for( RANGE_ITERATOR<d> node_iterator(RANGE<T_INDEX>(index, index+1));node_iterator.Valid();node_iterator.Next(),vertex++){
            const int node=engine->cells_mesh(cell)(vertex);
            if( node )
                PHYSBAM_ASSERT(node_to_index.Get_Or_Insert(node,node_iterator.Index())==node_iterator.Index());
        }
    }

    // Check for projection faults
    {
        LOG::SCOPE scope( "Projection..." );
#if 0        
        for( int w = 1; w <= d; w++){
            for( RANGE_ITERATOR<d> n_iter(engine->unpadded_node_domain);n_iter.Valid();n_iter.Next()){
                if( engine->node_is_dirichlet(n_iter.Index()) || !engine->node_is_active(n_iter.Index()) ){
                    d_state.x(w)(n_iter.Index()) = 0;
                }
            }
            for(int node=1; node <= engine->Number_Of_Mesh_Nodes(); node++){
                if( engine->node_is_dirichlet_mesh(node) || !engine->node_is_active_mesh(node) ) {
                    d_state.x_mesh(w)(node) = 0;
                }
            }
        }
#else
        krylov_system.Project( krylov_d );
#endif
    }

#if 0
    // Check for projection faults
    {
        LOG::SCOPE scope( "Checking for Projection Faults A ..." );
        for( int w = 1; w <= d; w++){
            for( RANGE_ITERATOR<d> n_iter(engine->unpadded_node_domain);n_iter.Valid();n_iter.Next()){
                if( engine->node_is_dirichlet(n_iter.Index()) )
                    LOG::cout << n_iter.Index() << "  " << d_state.x(w)(n_iter.Index()) << std::endl;
                PHYSBAM_ASSERT( !engine->node_is_dirichlet(n_iter.Index()) || (engine->node_is_dirichlet(n_iter.Index()) && d_state.x(w)(n_iter.Index()) == 0));
                PHYSBAM_ASSERT( engine->node_is_active(n_iter.Index()) || (!engine->node_is_active(n_iter.Index()) && d_state.x(w)(n_iter.Index()) == 0));
            }
            for(int node=1; node <= engine->Number_Of_Mesh_Nodes(); node++){
                PHYSBAM_ASSERT( !engine->node_is_dirichlet_mesh(node) || (engine->node_is_dirichlet_mesh(node) &&  d_state.x_mesh(w)(node) == 0));
                PHYSBAM_ASSERT(engine->node_is_active_mesh(node)  || (!engine->node_is_active_mesh(node) && d_state.x_mesh(w)(node) == 0));
            }
        }
    }
#endif

    // Compute Diagonal Effects from dynamic constraints // There are no static constraints at the moment

    VECTOR< const ARRAY<CONSTRAINT_SEGMENT<T,d> >*, 3 > constraint_v(&engine->dynamic_point_constraints,
                                                                     &engine->static_point_constraints,
                                                                     &engine->collision_constraints);
#if 1
    {
        LOG::SCOPE scope( "Adding Constraint Contributions..." );
        for(int ctype = 1; ctype <=3 ; ctype ++){            
            const ARRAY<CONSTRAINT_SEGMENT<T,d> >& constraints = *constraint_v(ctype);
            for(int c=1;c<=constraints.m;c++){
            
                T K;
                if( constraints(c).is_reference )
                    K = engine->collision_spring_constants(constraints(c).spring_coefficient_ptr);
                else
                    K = constraints(c).spring_coefficient;

                for(int e=1; e<=2; e++){
                    const CONSTRAINT_NODE<T,d>& endpoint = constraints(c).endpoints[e-1];
                    switch(endpoint.type){                       
                    case CONSTRAINT_NODE<T,d>::GRID_FIXED:
                        {
                            const T_INDEX& grid_index = endpoint.grid_index();
                            const T_STENCIL interpolation_stencil=engine->Multilinear_Interpolation_Stencil_Grid(T_INDEX(),endpoint.multilinear_coordinates());
                            for(T_CONST_STENCIL_ITERATOR stencil_iterator(interpolation_stencil);stencil_iterator.Valid();stencil_iterator.Next())
                                if( engine->node_is_active(stencil_iterator.Key()+grid_index) )
                                    for(int v=1;v<=d;v++)
                                        d_state.x(v)(stencil_iterator.Key()+grid_index)+= -K*stencil_iterator.Data()*stencil_iterator.Data();
                        }                                        
                        break;
                    
                    case CONSTRAINT_NODE<T,d>::MESH_FIXED:
                        {
                            const T_INDEX& grid_index = engine->cell_indices_mesh( endpoint.mesh_index() );
                            const T_STENCIL interpolation_stencil=engine->Multilinear_Interpolation_Stencil_Mesh(endpoint.multilinear_coordinates());
                            int flat_index=1;
                            for(T_CONST_STENCIL_ITERATOR stencil_iterator(interpolation_stencil);stencil_iterator.Valid();stencil_iterator.Next(),flat_index++){
                                int mesh_node = engine->cells_mesh(endpoint.mesh_index())(flat_index);
                                if( mesh_node == 0){
                                    if( engine->node_is_active(stencil_iterator.Key()+grid_index) )
                                        for(int v=1;v<=d;v++)  d_state.x(v)(stencil_iterator.Key()+grid_index)+=-K*stencil_iterator.Data()*stencil_iterator.Data();
                                }
                                else{
                                    if( engine->node_is_active_mesh(mesh_node) ) 
                                        for(int v=1;v<=d;v++)  d_state.x_mesh(v)(mesh_node)+=-K*stencil_iterator.Data()*stencil_iterator.Data();
                                }
                            }
                        }
                        break;
                    
                    case CONSTRAINT_NODE<T,d>::KINEMATIC:
                        // Do nothing. Free space locations do not receive forces.
                        break;
                    case CONSTRAINT_NODE<T,d>::KINEMATIC_PTR:
                        // Do nothing. Free space locations do not receive forces.
                        break;
                    }
                }
            }
        }
    }
#endif

    // Invert and Negate to create a positive D^(-1)
    {
        LOG::SCOPE scope( "Invert and Negate..." );
        RANGE<T_INDEX> domain = engine->unpadded_node_domain;
        for( int w=1; w <= d; w++){
            for( RANGE_ITERATOR<d> n_iter(domain);n_iter.Valid();n_iter.Next()){
                if(engine->node_is_active(n_iter.Index())){
                    d_state.x(w)(n_iter.Index()) = d_state.x(w)(n_iter.Index()) == 0 ? 0 : -1.0/d_state.x(w)(n_iter.Index());
                    //PHYSBAM_ASSERT( d_state.x(w)(n_iter.Index()) == d_state.x(w)(n_iter.Index()) );
                }
            }
            for(int node=1; node <= engine->Number_Of_Mesh_Nodes(); node++){
                if(engine->node_is_active_mesh(node)){
                    const T_INDEX& nindex = node_to_index.Get( node ); 
                    d_state.x_mesh(w)(node) = d_state.x_mesh(w)(node) == 0 ? 0 : -1.0/d_state.x_mesh(w)(node);
                    //PHYSBAM_ASSERT( d_state.x_mesh(w)(node) == d_state.x_mesh(w)(node) );
                }
            }
        }
    }
    
#if 0    
    // Check for projection faults
    {
        LOG::SCOPE scope( "Checking for Projection Faults B ..." );
        for( int w = 1; w <= d; w++){
            for( RANGE_ITERATOR<d> n_iter(engine->unpadded_node_domain);n_iter.Valid();n_iter.Next()){
                if( engine->node_is_dirichlet(n_iter.Index()) || !engine->node_is_active(n_iter.Index()) )
                    PHYSBAM_ASSERT( d_state.x(w)(n_iter.Index()) == 0);
            }
            for(int node=1; node <= engine->Number_Of_Mesh_Nodes(); node++){
                if( engine->node_is_dirichlet_mesh(node) || !engine->node_is_active_mesh(node) ) 
                    PHYSBAM_ASSERT( d_state.x_mesh(w)(node) == 0);
            }
        }
    }
#endif

#if 0
    // Check for sign faults
    {
        LOG::SCOPE scope( "Checking for Sign Faults..." );
        for( int w = 1; w <= d; w++){
            for( RANGE_ITERATOR<d> n_iter(engine->unpadded_node_domain);n_iter.Valid();n_iter.Next()){
                PHYSBAM_ASSERT( d_state.x(w)(n_iter.Index()) >= 0);
            }
            for(int node=1; node <= engine->Number_Of_Mesh_Nodes(); node++){
                PHYSBAM_ASSERT( d_state.x_mesh(w)(node) >= 0);
            }
        }
    }
#endif


    {
        LOG::SCOPE scope( "Modifing for Stencils..." );
        RANGE<T_INDEX> domain = engine->unpadded_node_domain;
        RANGE<T_INDEX> offset_range( s_origin, s_origin );
        RANGE<T_INDEX> shift_domain = domain + offset_range;
        for( int w = 1; w <= d; w++){
            for( RANGE_ITERATOR<d> stencil_iterator(domain); stencil_iterator.Valid();stencil_iterator.Next()){
                const T_INDEX nindex=stencil_iterator.Index() - s_origin;
                const T_INDEX& rindex=stencil_iterator.Index();

                

                if(!build_complement && ((nindex.x % radius == 0) || (nindex.y % radius == 0) ||
                                         (nindex.z % radius == 0))){
                    bool exclude = true;
                    if(
                       (nindex.x % 2 == 0) &&
                       (nindex.y % 2 == 0) &&
                       (nindex.x != nindex.y) )
                        exclude = false;

                    if(
                       (nindex.x % 2 == 0) &&
                       (nindex.z % 2 == 0) &&
                       (nindex.x != nindex.z) )
                        exclude = false;

                    if(
                       (nindex.z % 2 == 0) &&
                       (nindex.y % 2 == 0) &&
                       (nindex.z != nindex.y) )
                        exclude = false;

                    if( exclude )
                        d_state.x(w)(rindex) = 0;
                }
                else if(build_complement && !((nindex.x % radius == 0) || (nindex.y % radius == 0) ||
                                              (nindex.z % radius == 0))){
                    d_state.x(w)(rindex) = 0;                        
                }
            }

            for(int node=1; node <= engine->Number_Of_Mesh_Nodes(); node++){
                if(engine->node_is_active_mesh(node)){
                    const T_INDEX nindex = node_to_index.Get( node ) - s_origin; 
                    const T_INDEX& rindex = node_to_index.Get( node ); 
                    if(!build_complement && ((nindex.x % radius == 0) || (nindex.y % radius == 0) || (nindex.z % radius == 0))){
                        bool exclude = true;
                        if(
                           (nindex.x % 2 == 0) &&
                           (nindex.y % 2 == 0) &&
                           (nindex.x != nindex.y) )
                            exclude = false;
                        
                        if(
                           (nindex.x % 2 == 0) &&
                           (nindex.z % 2 == 0) &&
                           (nindex.x != nindex.z) )
                            exclude = false;
                        
                        if(
                           (nindex.z % 2 == 0) &&
                           (nindex.y % 2 == 0) &&
                           (nindex.z != nindex.y) )
                            exclude = false;
                        
                        if( exclude )
                            d_state.x_mesh(w)(node) = 0;
                    }
                    else if(build_complement && !((nindex.x % radius == 0) || (nindex.y % radius == 0) || (nindex.z % radius == 0))){
                        d_state.x_mesh(w)(node) = 0;
                    }
                }           
            }    
        }
        
    }

    krylov_b_debug.Clear();
    krylov_b2.Clear();
}

template<class BINDER, class STATE, class DISCRETIZATION> void
ELASTIC_SOLVER<BINDER,STATE,DISCRETIZATION>::Initialize(int partitions)
{
    if( solver_initialized ){
        // PHYSBAM_FATAL_ERROR( "Can not initialize solver more than once." );
        return;
    }

    // Internal State
    LOG::cout << "Initialize Internal State..." << std::endl;
    engine->Initialize_State( *u_bind, *u_internal);
    engine->Initialize_State( *f_bind, *f_internal);
    
    // Debugging
    LOG::cout << "Initialize Debug temporaries..." << std::endl;
    //engine->Initialize_State(*b_debug_b,*b_debug);
    //engine->Initialize_State(*b2_debug_b,*b2_debug);

    // Preconditioner
    LOG::cout << "Initialize Preconditioner temporaries..." << std::endl;
    //engine->Initialize_State(*D_b,*D);
    //engine->Initialize_State(*D_primal_b,*D_primal);
    //engine->Initialize_State(*D_dual_b,*D_dual);
    //engine->Initialize_State(*omega_b,*omega);
    
    LOG::cout << "Initialize Additional data..." << std::endl;
    engine->Initialize_Domain_Structures();
    engine->Initialize_Blocks(partitions);
    engine->Initialize_Undeformed_Configuration(*u_internal);

    // Krylov-related temporaries
    LOG::cout << "Initialize Krylov-related temporaries..." << std::endl;
    engine->Initialize_State(*x_b,*x);
    engine->Initialize_State(*b_b,*b);
    engine->Initialize_State(*q_b,*q);
    engine->Initialize_State(*s_b,*s);
    //engine->Initialize_State(*r_b,*r);
    //engine->Initialize_State(*k_b,*k);
    //engine->Initialize_State(*z_b,*z);
    
    engine->SetAllowBoundaryCells(false);
    engine->SetFirstOrder(true);      
    solver_initialized = true;
}


template<class BINDER, class STATE, class DISCRETIZATION> float
ELASTIC_SOLVER<BINDER,STATE,DISCRETIZATION>::Exact_Solve( int krylov_iterations, int newton_iterations, T krylov_tolerance, T newton_tolerance, bool no_cut_cells)
{
    LOG::SCOPE scope("ELASTIC_SOLVER::Exact_Solve()");


    if( !solver_initialized )
        PHYSBAM_FATAL_ERROR( "Can not perform exact solve until solver is initialized.");

    static int frame = 0;
    STREAM_TYPE stream_type((float)0);
    CG_SYSTEM<DISCRETIZATION> krylov_system(*engine);
    CG_VECTOR<CG_POLICY<DISCRETIZATION> > krylov_u(*engine,*u_internal);
    CG_VECTOR<CG_POLICY<DISCRETIZATION> > krylov_b(*engine,*b);
    CG_VECTOR<CG_POLICY<DISCRETIZATION> > krylov_q(*engine,*q);
    CG_VECTOR<CG_POLICY<DISCRETIZATION> > krylov_s(*engine,*s);
    CG_VECTOR<CG_POLICY<DISCRETIZATION> > krylov_x(*engine,*x);
    //CG_VECTOR<CG_POLICY<DISCRETIZATION> > krylov_r(*engine,*r);
    //CG_VECTOR<CG_POLICY<DISCRETIZATION> > krylov_k(*engine,*k);
    //CG_VECTOR<CG_POLICY<DISCRETIZATION> > krylov_z(*engine,*z);


    //CG_VECTOR<CG_POLICY<DISCRETIZATION> > DiagonalMatrix(*engine,*D);
    //CG_VECTOR<CG_POLICY<DISCRETIZATION> > DiagonalPrimalMatrix(*engine,*D_primal);
    //CG_VECTOR<CG_POLICY<DISCRETIZATION> > DiagonalDualMatrix(*engine,*D_dual);
    //CG_VECTOR<CG_POLICY<DISCRETIZATION> > Omega(*engine,*omega);


    //krylov_system.SetDiagonalMatrix(&DiagonalMatrix);
    //krylov_system.SetDiagonalPrimalMatrix(&DiagonalPrimalMatrix);
    //krylov_system.SetDiagonalDualMatrix(&DiagonalDualMatrix);
    //krylov_system.SetOmega(&Omega);

    T norm;
    for(int newton_iteration=1;newton_iteration<=newton_iterations;newton_iteration++){
        LOG::SCOPE scope("Newton iteration");


        {
            LOG::SCOPE scope("Computing off-balance forces");
            krylov_b.Clear();

            //DiagonalMatrix.Clear();
            //DiagonalPrimalMatrix.Clear();
            //DiagonalDualMatrix.Clear();

            //Omega.Set(.5);
            //UpdateOmega();

            engine->Update_Position_Based_State(*u_internal, *D);  
        
            //DiagonalPrimalMatrix += DiagonalMatrix;
            //DiagonalDualMatrix += DiagonalMatrix;
        
            //Extract_D( *D, -1,         T_INDEX(), 10000000, false); 
            //Extract_D( *D_primal, -1,     T_INDEX(0,0,0), 4,  false);
            //Extract_D( *D_dual, -1,        T_INDEX(2,2,2), 4,  false);

            //DiagonalMatrix *= Omega;
            //DiagonalPrimalMatrix *= Omega;
            //DiagonalDualMatrix *= Omega;

            engine->Add_Force(*u_internal, *b);
        
            norm = krylov_system.Convergence_Norm(krylov_b);
            LOG::cout << norm << std::endl;
            if( norm < newton_tolerance)
                break;
	    }

        {
            LOG::SCOPE scope_new("Zeroing out displacements");
            krylov_x.Clear();
        }

        {
#define USE_CG
#ifndef USE_CG
            LOG::SCOPE scope_new("Symmetric QMR");
            SYMMQMR<T> symmqmr;
            symmqmr.restart_iterations=150;
            symmqmr.print_diagnostics=true;
            symmqmr.print_residuals=true;
            symmqmr.nullspace_tolerance=0;
            symmqmr.Solve(krylov_system,
                          krylov_x,krylov_b,krylov_q,krylov_s,krylov_r,krylov_k,krylov_z,
                          krylov_tolerance,0,krylov_iterations);
#else
            LOG::SCOPE scope_new("Conjugate Gradient");
            CONJUGATE_GRADIENT<T> cg;
            cg.restart_iterations=100;
            cg.print_diagnostics=true;
            cg.print_residuals=true;
            //symmqmr.nullspace_tolerance=0;
            cg.Solve(krylov_system,
                     krylov_x,krylov_b,krylov_q,krylov_s,krylov_b,krylov_b,krylov_b,
                     krylov_tolerance,0,krylov_iterations);            
#endif
        }

        {
            LOG::SCOPE scope_new("Updating displacements");
            krylov_u += krylov_x;
        }       

        {
            LOG::SCOPE scope( "Computing Final Off-Balance Forces");
            norm = krylov_system.Convergence_Norm(krylov_b);
        }
 
    }
    //Write_Output_Minimal<DISCRETIZATION>(stream_type, *engine, *u_internal,"output",frame++, 0);
    return norm;
}

template<class BINDER, class STATE, class DISCRETIZATION> const STATE&
ELASTIC_SOLVER<BINDER,STATE,DISCRETIZATION>::U_const() const 
{
    return *u_internal;
}


template<class BINDER, class STATE, class DISCRETIZATION> const STATE&
ELASTIC_SOLVER<BINDER,STATE,DISCRETIZATION>::F_const() const
{
    return *f_internal;
}

template<class BINDER, class STATE, class DISCRETIZATION> STATE&
ELASTIC_SOLVER<BINDER,STATE,DISCRETIZATION>::U()
{
    return *u_internal;
}

template<class BINDER, class STATE, class DISCRETIZATION> STATE&
ELASTIC_SOLVER<BINDER,STATE,DISCRETIZATION>::F()
{
    return *f_internal;
}


//template class ELASTIC_SOLVER<BINDER<ALIGNED_ARRAY<float, int, 64>, 3, 16>, NONLINEAR_ELASTICITY_STATE<float, 3>, NONLINEAR_ELASTICITY<float, 3> >;

//template class ELASTIC_SOLVER<BINDER<ALIGNED_ARRAY<float, int, 64>, 3, 16>, NONLINEAR_ELASTICITY_STATE<float, 3>, SKINNING_NONLINEAR_ELASTICITY<float, 3, true,false> >;

template class ELASTIC_SOLVER<HYBRID_BINDER<ALIGNED_ARRAY<float, int, 64>, 3, 16>, HYBRID_NONLINEAR_ELASTICITY_STATE< SKINNING_NONLINEAR_ELASTICITY<float,3,true,false> >, HYBRID_NONLINEAR_ELASTICITY<SKINNING_NONLINEAR_ELASTICITY<float,3,true,false> > >;

template class ELASTIC_SOLVER<HYBRID_BINDER<ALIGNED_ARRAY<float, int, 64>, 3, 16>, HYBRID_NONLINEAR_ELASTICITY_STATE< SKINNING_NONLINEAR_ELASTICITY<float,3,true,true> >, HYBRID_NONLINEAR_ELASTICITY<SKINNING_NONLINEAR_ELASTICITY<float,3,true,true> > >;

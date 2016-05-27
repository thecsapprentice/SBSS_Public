//#####################################################################
// Copyright 2011-2013, Nathan Mitchell, Taylor Patterson, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class NONLINEAR_ELASTICITY
//#####################################################################
#ifndef __NONLINEAR_ELASTICITY__
#define __NONLINEAR_ELASTICITY__

#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Matrices/MATRIX_2X2.h>
#include <PhysBAM_Tools/Matrices/MATRIX_3X3.h>
#include <PhysBAM_Tools/Matrices/DIAGONAL_MATRIX_2X2.h>
#include <PhysBAM_Tools/Matrices/DIAGONAL_MATRIX_3X3.h>
#include <PhysBAM_Tools/Matrices/MATRIX_MXN.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND_VIEW.h>

#include <Common/POWER.h>
#include <Common/STENCIL.h>

#include <EngineInterface/CELL_TYPE.h>

#include "MATERIAL_MODEL.h"
#include "COROTATED.h"
#include "NEOHOOKEAN.h"
#include "ST_VENANT_KIRCHHOFF.h"
#include "MOONEY_RIVLIN.h"
#include "BIPHASIC.h"
#include "CG_POLICY.h"
#include "SPECIALIZED_KERNELS_DATA.h"

namespace PhysBAM{

#define USE_SPECIALIZED_KERNELS
#define USE_THREADED_KERNELS
//#define LOG_DETAILED_PERFORMANCE

//#define USE_CG

template<class T, int d>
    VECTOR< ARRAY_VIEW<T, VECTOR<int,d> >, d> operator+=(VECTOR< ARRAY_VIEW<T, VECTOR<int,d> >, d> base,
                                                         const VECTOR< ARRAY<T, VECTOR<int,d> >, d>& other){
    for( int i=1; i<=d; i++)
        base(i) += other(i);
    return base;}


template<class T_ARRAY, int DIM, int stride=16> 
    struct BINDER
    {
        static const int STRIDE=stride;
        VECTOR<T_ARRAY,DIM> x;
        T_ARRAY p;
    };

template<class T,int d>
struct NONLINEAR_ELASTICITY_STATE
{
    typedef VECTOR<T,d> TV;
    typedef VECTOR<int,d> T_INDEX; 
    typedef T SCALAR;
    static const int DIM=d;

    typedef ARRAY<SCALAR,T_INDEX> T_SCALAR_VARIABLE_RAW;
    typedef VECTOR<T_SCALAR_VARIABLE_RAW,DIM> T_VECTOR_VARIABLE_RAW;    

    typedef ARRAY_VIEW<SCALAR,T_INDEX> T_SCALAR_VARIABLE;
    typedef VECTOR<T_SCALAR_VARIABLE,DIM> T_VECTOR_VARIABLE;   

    T_VECTOR_VARIABLE x;
    T_SCALAR_VARIABLE p;

    NONLINEAR_ELASTICITY_STATE& operator+=( const NONLINEAR_ELASTICITY_STATE& other ){
        x += other.x;
        p += other.p;
    }

    NONLINEAR_ELASTICITY_STATE& operator*=( const NONLINEAR_ELASTICITY_STATE& other ){
        x *= other.x;
        p *= other.p;
    }

    template<class T_ARRAY, int stride>
        void Initialize(BINDER<T_ARRAY, DIM, stride>& binder, const RANGE<T_INDEX> domain, bool resize_binder=true )
    {
        //LOG::SCOPE scope("NONLINEAR_ELASTICITY_STATE::Initialize()");
        //LOG::cout << "Initializing with data stride: " << stride << std::endl;
        int size = (domain.Edge_Lengths()+1).Product() + stride*2;       

        if(resize_binder){
            for(int i=1;i<=d;i++){
                binder.x(i).Resize(size,true,false);
                binder.x(i).Fill(SCALAR(0));
            }
            //binder.p.Resize(size,true,false);
            //binder.p.Fill(SCALAR(0));
            binder.p.Resize(0,true,false);
            binder.p.Fill(SCALAR(0));
        }

        for(int i=1;i<=d;i++){
            T_SCALAR_VARIABLE t(domain, binder.x(i).base_pointer);
            //x(i) = t;
            T_SCALAR_VARIABLE::Exchange_Arrays(t,x(i));
        }
        T_SCALAR_VARIABLE t(domain, binder.p.base_pointer);
        //p = t;
        T_SCALAR_VARIABLE::Exchange_Arrays(t,p);
    }
};

#define __MATERIAL__ COROTATED

template<class T,int d>
class NONLINEAR_ELASTICITY
    :public MATERIAL_MODEL<__MATERIAL__<T,d> >
{
public:
    typedef __MATERIAL__<T,d> T_CONSTITUTIVE_MODEL;
    typedef MATERIAL_MODEL<T_CONSTITUTIVE_MODEL> T_MATERIAL;
    typedef typename CG_POLICY<NONLINEAR_ELASTICITY<T,d> >::STATE T_STATE;
    enum {DIM=d};
    typedef T SCALAR;
    
    static const bool supports_constraints = false;
    static const bool supports_muscles = false;

protected:
    typedef VECTOR<T,d> TV;
    typedef VECTOR<int,d> T_INDEX;
    typedef STENCIL<T,d> T_STENCIL;
protected:
    typedef GRID<TV> T_GRID;  

private:
    typedef ARRAY<CELL_TYPE,T_INDEX> T_CELL_TYPE_FIELD;
    typedef ARRAY<int,T_INDEX> T_FLAG;

protected:
    typedef ARRAY<T,T_INDEX> T_SCALAR_VARIABLE_RAW;
    typedef VECTOR<T_SCALAR_VARIABLE_RAW,d> T_VECTOR_VARIABLE_RAW;

    typedef ARRAY_VIEW<T,T_INDEX> T_SCALAR_VARIABLE_VIEW;
    typedef VECTOR<T_SCALAR_VARIABLE_VIEW,d> T_VECTOR_VARIABLE_VIEW;
    typedef ARRAY_VIEW<const T,T_INDEX> T_SCALAR_VARIABLE_VIEW_CONST;
    typedef VECTOR<T_SCALAR_VARIABLE_VIEW_CONST,d> T_VECTOR_VARIABLE_VIEW_CONST;
private:

    typedef ARRAY<MATRIX<T,d>,T_INDEX> T_MATRIX_ARRAY;
    typedef ARRAY<DIAGONAL_MATRIX<T,d>,T_INDEX> T_DIAGONAL_MATRIX_ARRAY;
    typedef ARRAY<ROTATED_STRESS_DERIVATIVE<T,d>,T_INDEX> T_STRESS_DERIVATIVE_ARRAY;

protected:
    enum WORKAROUND {
        vertices_per_cell=POWER<2,d>::value,
        cells_per_block=POWER<2,d>::value,
        lower_cell_padding=1,
        upper_cell_padding=3,
        lower_node_padding=1,
        upper_node_padding=2
    };

public:
    const T h;
    const T_INDEX n;
    const T_GRID grid;

    const RANGE<T_INDEX> unpadded_cell_domain;
    const RANGE<T_INDEX> unpadded_node_domain;
    const RANGE<T_INDEX> padded_cell_domain;
    const RANGE<T_INDEX> padded_node_domain;

public:
    T constant_mu;
    T constant_mu_10;
    T constant_mu_01;
    T constant_kappa;
    T constant_alpha;
    T cutoff_value;
    T stabilization_factor;

    T constant_mu_stab;
    T constant_kappa_max;
    T constant_alpha_squared_over_kappa;
    T constant_lambda;

    int constant_partitions;

public:
    MATRIX_MXN<T> G_One;
protected:
    VECTOR<T_STENCIL,d> cell_centered_derivative_operator;
    MATRIX_MXN<T> K_stab;

protected:
    bool allow_boundary_cells;
    bool first_order;

public:    
    // The following are initialized by external entities prior to Initialize_Domain()
    T_CELL_TYPE_FIELD cell_type;
    ARRAY<ARRAY<TV>,T_INDEX> quadrature_point_coordinates;
    ARRAY<ARRAY<T>,T_INDEX> quadrature_point_weights;
    // TODO: Include support for spatially varying material parameters

    T_FLAG node_is_active;
    T_FLAG node_is_dirichlet;


//#ifndef USE_SPECIALIZED_KERNELS
    T_MATRIX_ARRAY U;
    T_MATRIX_ARRAY V;
    T_DIAGONAL_MATRIX_ARRAY Sigma;
    T_DIAGONAL_MATRIX_ARRAY Q_hat;
    T_STRESS_DERIVATIVE_ARRAY dP_dF;
//#endif

    // Blocking data structures

    ARRAY<T_INDEX> cell_block_base_indices;
    ARRAY<int> cell_block_partition_offsets;
    ARRAY<int> cell_block_base_offsets;

    ARRAY<T_INDEX> node_block_base_indices;
    ARRAY<int> node_block_partition_offsets;
    ARRAY<int> node_block_base_offsets;

    // Flattened data structures

#ifdef USE_SPECIALIZED_KERNELS
    SPECIALIZED_KERNEL_DATA<T,d> specialized_data;
#endif

public:
    ARRAY<T> muscle_fiber_max_stresses;
    ARRAY<T> muscle_activations;
    ARRAY<ARRAY<int>, T_INDEX> cell_muscles;
    ARRAY<ARRAY<TV>, T_INDEX> cell_fibers;
    ARRAY<ARRAY<TV>, T_INDEX> cell_F_fibers;
    ARRAY<ARRAY<T>, T_INDEX> cell_densities;
    ARRAY<ARRAY<T>, T_INDEX> cell_c1;
    ARRAY<ARRAY<T>, T_INDEX> cell_c2;

//#####################################################################    
public:
    NONLINEAR_ELASTICITY(const T_INDEX& n_input,const T h_input,const TV origin);
private:
    // Initialization routines dependent on init-time constants
    void Build_Cell_Centered_Derivative_Operator();
    void Initialize_Cell_Centered_Quadrature();
    void Initialize_Stabilization_Kernel();

 protected:
    void CompactData_Specialized(T_VECTOR_VARIABLE_VIEW_CONST x, T_SCALAR_VARIABLE_VIEW_CONST p) const;
    void UnCompactData_Specialized(T_VECTOR_VARIABLE_VIEW x,T_SCALAR_VARIABLE_VIEW p) const;
    void CompactData_Specialized(T_VECTOR_VARIABLE_VIEW_CONST x) const;
    void UnCompactData_Specialized(T_VECTOR_VARIABLE_VIEW x) const;

public:
    static T_VECTOR_VARIABLE_VIEW_CONST View_Convert(const T_VECTOR_VARIABLE_VIEW& in);

    static int Get_Lower_Cell_Padding() {return lower_cell_padding;};
    static int Get_Upper_Cell_Padding() {return upper_cell_padding;};
    static int Get_Lower_Node_Padding() {return lower_node_padding;};
    static int Get_Upper_Node_Padding() {return upper_node_padding;};

    void Initialize_Parameters(const T mu,const T kappa,const T alpha,const T cutoff_value,const T stabilization_factor);
    void Initialize_Parameters(const T mu_10,const T mu_01,const T kappa,const T alpha,const T cutoff_value,const T stabilization_factor);
    void Initialize_Domain_Structures();
    void Initialize_Blocks_Specialized(const int number_of_partitions);
    void Initialize_Blocks(const int number_of_partitions);

    template<class T_ARRAY, int stride>
        void Initialize_State(BINDER<T_ARRAY, d, stride>& binder, T_STATE& state, bool resize_binder=true) const
    {
        PHYSBAM_ASSERT(padded_node_domain == padded_cell_domain);
        state.Initialize(binder, padded_node_domain, resize_binder);
    }


    static bool Supports_Constraints() {return supports_constraints;}
    static bool Supports_Muscles() {return supports_muscles;}

    // Set state functions
    void SetAllowBoundaryCells( bool allow_boundary_cells_input ) { allow_boundary_cells = allow_boundary_cells_input; } ;
    void SetFirstOrder( bool first_order_input ) { first_order = first_order_input; };

    // Update state functions
    void Update_Position_Based_State(const T_STATE& state_u, T_STATE& state_d);
    void Update_Position_Based_State_Elasticity(T_VECTOR_VARIABLE_VIEW_CONST u,T_SCALAR_VARIABLE_VIEW_CONST p, T_VECTOR_VARIABLE_VIEW diag);
    template<bool enable_constraints,bool enable_muscles>
    void Update_Position_Based_State_Specialized(T_VECTOR_VARIABLE_VIEW_CONST u,T_SCALAR_VARIABLE_VIEW_CONST p, T_VECTOR_VARIABLE_VIEW diag);

    // Force routines
    void Add_Force(const T_STATE& state_u, T_STATE& state_f);
    void Add_Force_Second_Order_Elasticity(T_VECTOR_VARIABLE_VIEW_CONST u,T_SCALAR_VARIABLE_VIEW_CONST p, T_VECTOR_VARIABLE_VIEW f,T_SCALAR_VARIABLE_VIEW q) const;
    void Add_Force_First_Order_Elasticity(T_VECTOR_VARIABLE_VIEW_CONST u,T_SCALAR_VARIABLE_VIEW_CONST p, T_VECTOR_VARIABLE_VIEW f,T_SCALAR_VARIABLE_VIEW q) const;
    template<bool enable_constraints,bool enable_muscles>
        void Add_Force_First_Order_Elasticity_Specialized(T_VECTOR_VARIABLE_VIEW_CONST u,T_SCALAR_VARIABLE_VIEW_CONST p, T_VECTOR_VARIABLE_VIEW f,T_SCALAR_VARIABLE_VIEW q) ;
    void Add_Force_Stab(T_VECTOR_VARIABLE_VIEW_CONST dx,T_VECTOR_VARIABLE_VIEW df,const bool first_order) const;
    void Add_Force_Differential(const T_STATE& state_u, T_STATE& state_f) const;
    void Add_Force_Differential(const T_STATE& state_u, T_STATE& state_f, const int subdomain) const;
    void Add_Force_Differential_Elasticity(T_VECTOR_VARIABLE_VIEW_CONST du, T_SCALAR_VARIABLE_VIEW_CONST dp,T_VECTOR_VARIABLE_VIEW df,T_SCALAR_VARIABLE_VIEW dq) const;
    template<bool enable_constraints,bool enable_muscles>
        void Add_Force_Differential_Elasticity_Specialized(T_VECTOR_VARIABLE_VIEW_CONST du, T_SCALAR_VARIABLE_VIEW_CONST dp,T_VECTOR_VARIABLE_VIEW df,T_SCALAR_VARIABLE_VIEW dq) const;
    template<bool enable_constraints,bool enable_muscles>
    void Add_Force_Differential_Elasticity_Specialized(T_VECTOR_VARIABLE_VIEW_CONST du, T_SCALAR_VARIABLE_VIEW_CONST dp,T_VECTOR_VARIABLE_VIEW df,T_SCALAR_VARIABLE_VIEW dq, const int subdomain) const;

//#####################################################################
};

}
#endif

#include "SKINNING_NONLINEAR_ELASTICITY.h"
#include <PhysBAM_Tools/Math_Tools/cube.h>
#include <PhysBAM_Tools/Matrices/SYMMETRIC_MATRIX_2X2.h>
#include <PhysBAM_Tools/Matrices/SYMMETRIC_MATRIX_3X3.h>
//#include <PhysBAM_Tools/Read_Write/Grids_Uniform_Arrays/READ_WRITE_ARRAYS.h>
//#include <PhysBAM_Tools/Read_Write/Matrices_And_Vectors/READ_WRITE_MATRIX.h>
//#include <PhysBAM_Tools/Read_Write/Arrays/READ_WRITE_ARRAY.h>
//#include <PhysBAM_Tools/Read_Write/Vectors/READ_WRITE_VECTOR_3D.h>
#include <Common/RANGE_ITERATOR.h>
#include <Common/STENCIL_ITERATOR.h>

#include <iostream>
#include <fstream>

using namespace PhysBAM;


//#####################################################################
// Function Update_Position_Based_State_Muscles
//#####################################################################

#ifndef USE_SPECIALIZED_KERNELS
#include <SIMD_Optimized_Kernels/Kernel_Wrappers/Muscle_Tension/Muscle_Tension.h>
#include <SIMD_Optimized_Kernels/Kernel_Wrappers/Central_Gradient/Central_Gradient.h>
#endif

template<class T,int d,bool enable_constraints,bool enable_muscles> void SKINNING_NONLINEAR_ELASTICITY<T,d,enable_constraints,enable_muscles>::
Update_Position_Based_State_Muscles(T_VECTOR_VARIABLE_VIEW_CONST u,T_SCALAR_VARIABLE_VIEW_CONST p)
{
    if(!enable_muscles || d==2) return;
    LOG::SCOPE scope("SKINNING_NONLINEAR_ELASTICITY::Update_Position_Based_State_Muscles()");

    // Muscles

#ifndef USE_SPECIALIZED_KERNELS

    for(RANGE_ITERATOR<d> cell_iterator(unpadded_cell_domain);cell_iterator.Valid();cell_iterator.Next()){
        const T_INDEX& cell_index=cell_iterator.Index();

        switch(this->cell_type(cell_index))
        {
            case BOUNDARY_CELL_TYPE:
                PHYSBAM_ASSERT(allow_boundary_cells);
            case INTERIOR_CELL_TYPE:
                {

                    MATRIX_MXN<T> Du;Du.Resize(d,8);
                    {int vertex=1;
                        for(RANGE_ITERATOR<d> iterator(RANGE<T_INDEX>(cell_index,cell_index+1));iterator.Valid();iterator.Next(),vertex++){
                            const T_INDEX& index=iterator.Index();
                            for(int v=1;v<=d;v++)
                                Du(v,vertex)=u(v)(index);}}
                    MATRIX<T,d> Fe;
                    Central_Gradient(Du,Fe,1.f/h);
                    Fe+=1;

                    for(int m=1;m<=cell_muscles(cell_index).m;m++){
                        const int muscle=this->cell_muscles(cell_index)(m);                
                        const TV& fiber=this->cell_fibers(cell_index)(m);
                        const T density=this->cell_densities(cell_index)(m);
                        const T activation=BASE::muscle_activations(muscle);
                        const T fiber_max_stress=BASE::muscle_fiber_max_stresses(muscle);
                        TV& F_fiber = this->cell_F_fibers(cell_index)(m);
                        F_fiber=Fe*fiber;
                        T stretch_squared=F_fiber.Magnitude_Squared(),stretch=sqrt(stretch_squared);
                        stretch = min((T)3.0f, stretch);
                        T tension, tension_derivative;
#ifndef USE_SPECIALIZED_KERNELS
                        PHYSBAM_NOT_IMPLEMENTED("TODO: Tension functions need to be fixed for non-specialized kernels.");
                        //Tension(tension, stretch,activation,density,fiber_max_stress);
                        //Tension_Derivative(tension_derivative, stretch,activation,density,fiber_max_stress);
#else
                        K_Tension(tension, stretch,activation,density,fiber_max_stress);
                        K_Tension_Derivative(tension_derivative, stretch,activation,density,fiber_max_stress);
#endif
                        tension_derivative = max( T(0), tension_derivative );

                        T& c1=cell_c1(cell_index)(m);
                        T& c2=cell_c2(cell_index)(m);
                        c1=tension/stretch;
                        c2=(tension_derivative-c1)/stretch_squared;

                        if( density == 0 ) { F_fiber(1) = 0; F_fiber(2) = 0; F_fiber(3) = 0; c1=0; c2=0; }
                    }
                }
                break;
            case EXTERIOR_CELL_TYPE:
            case DIRICHLET_CELL_TYPE:
                break;
            default:
                PHYSBAM_FATAL_ERROR("Improper cell type in SKINNING_NONLINEAR_ELASTICITY::Update_Position_Based_State");            
        }
    }

#endif
}
//#####################################################################
// Function Add_Force_First_Order_Muscles
//#####################################################################
template<class T,int d,bool enable_constraints,bool enable_muscles> void SKINNING_NONLINEAR_ELASTICITY<T,d,enable_constraints,enable_muscles>::
Add_Force_First_Order_Muscles(T_VECTOR_VARIABLE_VIEW_CONST u,T_SCALAR_VARIABLE_VIEW_CONST p, T_VECTOR_VARIABLE_VIEW f,T_SCALAR_VARIABLE_VIEW q) const
{
    if(!enable_muscles || d==2) return;
    LOG::SCOPE scope("SKINNING_NONLINEAR_ELASTICITY::Add_Force_First_Order_Muscles()");

    // Muscles
#ifndef USE_SPECIALIZED_KERNELS
    const T cell_volume=pow<d>(h);
    for(RANGE_ITERATOR<d> iterator(unpadded_cell_domain);iterator.Valid();iterator.Next()){
        const T_INDEX& cell_index=iterator.Index();
        if(this->cell_type(cell_index)==INTERIOR_CELL_TYPE || (/*first_order &&*/ this->cell_type(cell_index)==BOUNDARY_CELL_TYPE)){

            for(int m=1;m<=this->cell_muscles(cell_index).m;m++){
                const T c1=this->cell_c1(cell_index)(m);
                const TV& fiber=this->cell_fibers(cell_index)(m);
                const TV& F_fiber = this->cell_F_fibers(cell_index)(m);
                MATRIX<T,d> Pe=MATRIX<T,d>::Outer_Product(F_fiber,fiber)*c1;
                //MATRIX<T,d> Pe=M*c1;

                for(int v=1;v<=d;v++)
                    for(int w=1;w<=d;w++)
                        Accumulation(f(v),cell_centered_derivative_operator(w),cell_index,-Pe(v,w)*cell_volume);
            }

        }
        else if(this->cell_type(cell_index)==BOUNDARY_CELL_TYPE)
            PHYSBAM_FATAL_ERROR("Cut cells not yet supported with muscles");
    }
 #endif

}
//#####################################################################
// Function Add_Force_Second_Order_Muscles
//#####################################################################
template<class T,int d,bool enable_constraints,bool enable_muscles> void SKINNING_NONLINEAR_ELASTICITY<T,d,enable_constraints,enable_muscles>::
Add_Force_Second_Order_Muscles(T_VECTOR_VARIABLE_VIEW_CONST u,T_SCALAR_VARIABLE_VIEW_CONST p, T_VECTOR_VARIABLE_VIEW f,T_SCALAR_VARIABLE_VIEW q) const
{
    if(!enable_muscles || d==2) return;
    LOG::SCOPE scope("SKINNING_NONLINEAR_ELASTICITY::Add_Force_Second_Order_Muscles()");

    const T cell_volume=(d==2)?sqr(h):cube(h);

    // Muscles

    for(RANGE_ITERATOR<d> iterator(unpadded_cell_domain);iterator.Valid();iterator.Next()){
        const T_INDEX& cell_index=iterator.Index();
        if(this->cell_type(cell_index)==INTERIOR_CELL_TYPE || (/*first_order &&*/ this->cell_type(cell_index)==BOUNDARY_CELL_TYPE)){

            for(int m=1;m<=this->cell_muscles(cell_index).m;m++){
                const T c1=this->cell_c1(cell_index)(m);
                const TV& fiber=this->cell_fibers(cell_index)(m);
                const TV& F_fiber = this->cell_F_fibers(cell_index)(m);
                //const MATRIX<T,d>& M=cell_M(cell_index)(m);
                MATRIX<T,d> Pe=MATRIX<T,d>::Outer_Product(F_fiber,fiber)*c1;
                //MATRIX<T,d> Pe=M*c1;

                for(int v=1;v<=d;v++)
                    for(int w=1;w<=d;w++)
                        Accumulation(f(v),cell_centered_derivative_operator(w),cell_index,-Pe(v,w)*cell_volume);
            }

        }
        else if(this->cell_type(cell_index)==BOUNDARY_CELL_TYPE)
            PHYSBAM_FATAL_ERROR("Cut cells not yet supported with muscles");
    }
 

}
//#####################################################################
// Function Add_Force_Differential_Muscles
//#####################################################################
template<class T,int d,bool enable_constraints,bool enable_muscles> void SKINNING_NONLINEAR_ELASTICITY<T,d,enable_constraints,enable_muscles>::
Add_Force_Differential_Muscles(T_VECTOR_VARIABLE_VIEW_CONST du,T_SCALAR_VARIABLE_VIEW_CONST dp,T_VECTOR_VARIABLE_VIEW df,T_SCALAR_VARIABLE_VIEW dq) const
{
    if(!enable_muscles || d==2) return;
    // LOG::SCOPE scope("SKINNING_NONLINEAR_ELASTICITY::Add_Force_Differential_Muscles()");
    // Muscles
#ifndef USE_SPECIALIZED_KERNELS
    const T cell_volume=pow<d>(h);
    for(RANGE_ITERATOR<d> iterator(unpadded_cell_domain);iterator.Valid();iterator.Next()){
        const T_INDEX& cell_index=iterator.Index();
        if(cell_type(cell_index)==INTERIOR_CELL_TYPE || cell_type(cell_index)==BOUNDARY_CELL_TYPE){
            
            MATRIX_MXN<T> Du;Du.Resize(d,BASE::vertices_per_cell);
            {int vertex=1;
                for(RANGE_ITERATOR<d> iterator(RANGE<T_INDEX>(cell_index,cell_index+1));iterator.Valid();iterator.Next(),vertex++){
                    const T_INDEX& index=iterator.Index();
                    for(int v=1;v<=d;v++)
                        Du(v,vertex)=du(v)(index);}}
            MATRIX<T,d> dF;
            Central_Gradient(Du,dF,1.f/h);

            //for(int v=1;v<=d;v++)
            //    for(int w=1;w<=d;w++)
            //        dF(v,w)=cell_centered_derivative_operator(w).Contraction(du(v),cell_index);

            for(int m=1;m<=cell_muscles(cell_index).m;m++){
                const TV& fiber=cell_fibers(cell_index)(m);
                const TV& F_fiber = cell_F_fibers(cell_index)(m);
                const T c1=cell_c1(cell_index)(m);
                const T c2=cell_c2(cell_index)(m);
                //const MATRIX<T,d>& M=cell_M(cell_index)(m);
                //MATRIX<T,d> dP_fiber=c2*MATRIX<T,d>::Inner_Product(M,dF)*M+c1*MATRIX<T,d>::Outer_Product(dF*fiber,fiber);
                MATRIX<T,d> M = MATRIX<T,d>::Outer_Product(F_fiber,fiber);
                MATRIX<T,d> dP_fiber=c2*MATRIX<T,d>::Inner_Product(M,dF)*M+c1*MATRIX<T,d>::Outer_Product(dF*fiber,fiber);
                
                for(int v=1;v<=d;v++)
                    for(int w=1;w<=d;w++)
                        Accumulation(df(v),cell_centered_derivative_operator(w),cell_index,-dP_fiber(v,w)*cell_volume);
            }

        }
        else if(cell_type(cell_index)==BOUNDARY_CELL_TYPE)
            PHYSBAM_FATAL_ERROR("Cut cells not yet supported with muscles");        
    }

#endif
}


template class SKINNING_NONLINEAR_ELASTICITY<float,3,true,true>;
template class SKINNING_NONLINEAR_ELASTICITY<float,3,true,false>;
template class SKINNING_NONLINEAR_ELASTICITY<float,2,true,true>;
template class SKINNING_NONLINEAR_ELASTICITY<float,2,true,false>;
#ifndef USE_SPECIALIZED_KERNELS
template class SKINNING_NONLINEAR_ELASTICITY<double,3,true,true>;
template class SKINNING_NONLINEAR_ELASTICITY<double,3,true,false>;
template class SKINNING_NONLINEAR_ELASTICITY<double,2,true,true>;
template class SKINNING_NONLINEAR_ELASTICITY<double,2,true,false>;
#endif

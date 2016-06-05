//#####################################################################
//  Copyright (c) 2011-2013 Nathan Mitchell, Eftychios Sifakis.
//  This file is covered by the FreeBSD license. Please refer to the 
//  license.txt file for more information.
//#####################################################################


#include <PhysBAM_Tools/Matrices/MATRIX_3X3.h>

using namespace PhysBAM;

template<class T>
void Tension_Reference(T& tension, const T stretch, const T activation, const T density, const T fiber_max_stress)
{
    T fiber_p1=(T).05;
    T fiber_p2=(T)6.6;
    T fiber_cutoff=(T)1.4;
    T cutoff_scaled=fiber_p2*(fiber_cutoff-1);
    T fiber_p3=fiber_p1*fiber_p2*(exp(cutoff_scaled)-1);
    T fiber_p4=fiber_p1*(exp(cutoff_scaled)*(1-fiber_p2*fiber_cutoff)+fiber_p2-1);

    T strain=stretch-1,strain_abs=abs(strain),active_tension=0,passive_tension=0,scale=(T)25/(T)6;
    if(stretch>fiber_cutoff)passive_tension=fiber_p3*stretch+fiber_p4;else if(stretch>1)passive_tension=fiber_p1*(exp(fiber_p2*strain)-fiber_p2*strain-1);
    if(strain_abs<.4)active_tension=activation*density*(1-scale*sqr(strain));else if(strain_abs<.6)active_tension=2*scale*activation*density*sqr(strain_abs-(T).6);
    tension=fiber_max_stress*(active_tension+passive_tension);
}

template<class T>
void Tension_Derivative_Reference(T& tension_derivative, const T stretch, const T activation, const T density, const T fiber_max_stress)
{
    T fiber_p1=(T).05;
    T fiber_p2=(T)6.6;
    T fiber_cutoff=(T)1.4;
    T cutoff_scaled=fiber_p2*(fiber_cutoff-1);
    T fiber_p3=fiber_p1*fiber_p2*(exp(cutoff_scaled)-1);

    T strain=stretch-1,strain_abs=abs(strain),active_tension_derivative=0,passive_tension_derivative=0,scale=(T)25/(T)6;
    if(stretch>fiber_cutoff)passive_tension_derivative=fiber_p3;else if(stretch>1)passive_tension_derivative=fiber_p1*fiber_p2*(exp(fiber_p2*strain)-1);
    if(strain_abs<.4)active_tension_derivative=-2*scale*activation*density*strain;else if(strain_abs<.6)active_tension_derivative=4*scale*activation*density*(strain-sign(strain)*(T).6);
    tension_derivative=fiber_max_stress*(active_tension_derivative+passive_tension_derivative);
}

template<class T>
bool Tension_Compare(const T tension, const T tension_reference)
{
    std::cout<<"tension = "<<tension<<std::endl;
    std::cout<<"tension_reference = "<<tension_reference<<std::endl;
    std::cout<<"Difference="<<fabs(tension-tension_reference)<<std::endl;;
    return fabs(tension-tension_reference)<(T)1e-6;
}

template<class T>
bool Tension_Derivative_Compare(const T tension_derivative, const T tension_derivative_reference)
{
    std::cout<<"tension_derivative = "<<tension_derivative<<std::endl;
    std::cout<<"tension_derivative_reference = "<<tension_derivative_reference<<std::endl;
    std::cout<<"Difference="<<fabs(tension_derivative-tension_derivative_reference)<<std::endl;
    return fabs(tension_derivative-tension_derivative_reference)<(T)1e-6;
}

template void Tension_Reference(float& tension, const float stretch, const float activation, const float density, const float fiber_max_stress);
template void Tension_Derivative_Reference(float& tension_derivative, const float stretch, const float activation, const float density, const float fiber_max_stress);
template bool Tension_Compare(const float tension, const float tension_reference);
template bool Tension_Derivative_Compare(const float tension_derivative, const float tension_derivative_reference);

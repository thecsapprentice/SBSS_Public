//#####################################################################
// Copyright 2002-2008, Ronald Fedkiw, Geoffrey Irving, Igor Neverov, Craig Schroeder, Tamar Shinar, Eftychios Sifakis, Huamin Wang, Rachel Weinstein.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SYMMCGNE
//#####################################################################
#include <PhysBAM_Tools/Krylov_Solvers/SYMMCGNE.h>
#include <PhysBAM_Tools/Krylov_Solvers/KRYLOV_VECTOR_BASE.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <cfloat>
#include <limits>
using namespace PhysBAM;
//#####################################################################
// Destructor
//#####################################################################
template<class T> SYMMCGNE<T>::
~SYMMCGNE()
{}
//#####################################################################
// Function Solve
//#####################################################################
template<class T> bool SYMMCGNE<T>::
Solve(const KRYLOV_SYSTEM_BASE<T>& system,KRYLOV_VECTOR_BASE<T>& x,const KRYLOV_VECTOR_BASE<T>& b,KRYLOV_VECTOR_BASE<T>& q,KRYLOV_VECTOR_BASE<T>& s,KRYLOV_VECTOR_BASE<T>& r,
    KRYLOV_VECTOR_BASE<T>& k,KRYLOV_VECTOR_BASE<T>& z,const T tolerance,const int min_iterations,const int max_iterations)
{
    // NOTE: you should never try to make copies of VECTOR_T's inside here as they could be indirect.
    static const T small_number=std::numeric_limits<T>::epsilon();
    system.Set_Boundary_Conditions(x);
    T rho_old=(T)FLT_MAX;T convergence_norm=0;
    int iterations;for(iterations=0;;iterations++){
        bool restart=!iterations || (restart_iterations && iterations%restart_iterations==0);
        if(restart){
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
            if(print_residuals) LOG::cout<<"restarting cg"<<std::endl;
#endif
            r=b;system.Multiply(x,q);r-=q;system.Project(r);}
        // stopping conditions
        system.Project_Nullspace(r);
        convergence_norm=system.Convergence_Norm(r);
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
        if(print_residuals) LOG::cout<<convergence_norm<<std::endl;
#endif
        if(convergence_norm<=tolerance && (iterations>=min_iterations || convergence_norm<small_number)){
            if(print_diagnostics) LOG::Stat("cg iterations",iterations);if(iterations_used) *iterations_used=iterations;return true;}
        if(iterations==max_iterations) break;
        // actual iteration
        T rho=(T)system.Inner_Product(r,r);
        if(restart){
            system.Multiply(r,s);
            system.Project(s);}
        else{
            system.Multiply(r,q);
            system.Project(q);
            s.Copy(rho/rho_old,s,q);}
        T s_dot_s=(T)system.Inner_Product(s,s);
        T alpha=s_dot_s?rho/s_dot_s:(T)FLT_MAX;
        x.Copy(alpha,s,x);
        system.Multiply(s,q);
        system.Project(q);
        r.Copy(-alpha,q,r);
        rho_old=rho;}

#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
    if(print_diagnostics) LOG::Stat("cg iterations",iterations);if(iterations_used) *iterations_used=iterations;
    if(print_diagnostics) LOG::cout<<"cg not converged after "<<max_iterations<<" iterations, error = "<<convergence_norm<<std::endl;
#endif
    return false;
}
//#####################################################################
template class SYMMCGNE<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class SYMMCGNE<double>;
#endif

//#####################################################################
// Copyright 2011, Taylor Patterson, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MATERIAL_MODEL
//#####################################################################
#ifndef __MATERIAL_MODEL__
#define __MATERIAL_MODEL__

#include "ROTATED_STRESS_DERIVATIVE.h"

namespace PhysBAM{

#define DEBUG_NONMIXED

template<class T_MATERIAL>
class MATERIAL_MODEL:public T_MATERIAL
{
private:
    typedef typename T_MATERIAL::SCALAR T;
    enum {d=T_MATERIAL::dim};
    typedef T_MATERIAL BASE;

public:
    static bool Is_Mooney_Rivlin()
    {return T_MATERIAL::is_mooney_rivlin;}

    // Non Mooney-Rivlin
    static DIAGONAL_MATRIX<T,d> P_hat(const DIAGONAL_MATRIX<T,d>& Sigma,const DIAGONAL_MATRIX<T,d>& Q_hat,const T p,const T mu,const T alpha,const T kappa)
    {
#ifndef DEBUG_NONMIXED
       return BASE::Deviatoric_P_hat(Sigma,mu)+alpha*p*Q_hat;
#else
       return BASE::Deviatoric_P_hat(Sigma,mu)+(kappa*M(Sigma))*Q_hat;
#endif
    }

    // Mooney-Rivlin
    static DIAGONAL_MATRIX<T,d> P_hat(const DIAGONAL_MATRIX<T,d>& Sigma,const DIAGONAL_MATRIX<T,d>& Q_hat,const T p,const T mu_10,const T mu_01,const T alpha,const T kappa)
    {
#ifndef DEBUG_NONMIXED
        return BASE::Deviatoric_P_hat(Sigma,mu_10,mu_01)+alpha*p*Q_hat;
#else
        PHYSBAM_NOT_IMPLEMENTED();
#endif
    }

    static T q(const DIAGONAL_MATRIX<T,d>& Sigma,const T p,const T alpha,const T alpha_squared_over_kappa)
    {
#ifndef DEBUG_NONMIXED
        return -alpha*BASE::M(Sigma)+alpha_squared_over_kappa*p;
#else
        return -alpha_squared_over_kappa*p;
#endif
    }

    static MATRIX<T,d> dP_hat(const ROTATED_STRESS_DERIVATIVE<T,d>& rotated_dPdF,const DIAGONAL_MATRIX<T,d>& Q_hat,const MATRIX<T,d>& dF_hat,const T dp,const T alpha)
    {
#ifndef DEBUG_NONMIXED
        return rotated_dPdF.dP_hat(dF_hat)+alpha*dp*Q_hat;
#else
        return rotated_dPdF.dP_hat(dF_hat);
#endif
    }

    static T dq(const DIAGONAL_MATRIX<T,d>& Q_hat,const MATRIX<T,d>& dF_hat,const T dp,const T alpha,const T alpha_squared_over_kappa)
    {
#ifndef DEBUG_NONMIXED
        return -alpha*Q_hat.Times_Transpose(dF_hat).Trace()+alpha_squared_over_kappa*dp;
#else
        return -alpha_squared_over_kappa*dp;
#endif
    }

    // Non Mooney-Rivlin
    static ROTATED_STRESS_DERIVATIVE<T,d> Isotropic_Stress_Derivative(const DIAGONAL_MATRIX<T,d>& Sigma,const T p,const T mu,const T kappa,const T alpha,const bool apply_definiteness_fix)
    {
        PHYSBAM_ASSERT(!Is_Mooney_Rivlin());
#ifndef DEBUG_NONMIXED
        ROTATED_STRESS_DERIVATIVE<T,d> rotated_dPdF=BASE::Augmented_Rotated_Stress_Derivative(Sigma,p,mu,alpha);
#else
        ROTATED_STRESS_DERIVATIVE<T,d> rotated_dPdF=BASE::Rotated_Stress_Derivative(Sigma,mu,kappa);
        rotated_dPdF.Make_Positive_Definite();
        return rotated_dPdF;
#endif
        ROTATED_STRESS_DERIVATIVE<T,d> rotated_dPdF_indefinite=BASE::Rotated_Stress_Derivative(Sigma,mu,kappa);
        ROTATED_STRESS_DERIVATIVE<T,d> rotated_dPdF_incremental_fix=rotated_dPdF_indefinite;
        rotated_dPdF_incremental_fix.Make_Positive_Definite();
        rotated_dPdF_incremental_fix-=rotated_dPdF_indefinite;
        if(apply_definiteness_fix) rotated_dPdF+=rotated_dPdF_incremental_fix;

        return rotated_dPdF;
    }

    // Mooney-Rivlin
    static ROTATED_STRESS_DERIVATIVE<T,d> Isotropic_Stress_Derivative(const DIAGONAL_MATRIX<T,d>& Sigma,const T p,const T mu_10,const T mu_01,const T kappa,const T alpha,const bool apply_definiteness_fix)
    {
        PHYSBAM_ASSERT(Is_Mooney_Rivlin());

        ROTATED_STRESS_DERIVATIVE<T,d> rotated_dPdF=BASE::Augmented_Rotated_Stress_Derivative(Sigma,p,mu_10,mu_01,alpha);
        ROTATED_STRESS_DERIVATIVE<T,d> rotated_dPdF_indefinite=BASE::Rotated_Stress_Derivative(Sigma,mu_10,mu_01,kappa);
        ROTATED_STRESS_DERIVATIVE<T,d> rotated_dPdF_incremental_fix=rotated_dPdF_indefinite;
        rotated_dPdF_incremental_fix.Make_Positive_Definite();
        rotated_dPdF_incremental_fix-=rotated_dPdF_indefinite;
        if(apply_definiteness_fix) rotated_dPdF+=rotated_dPdF_incremental_fix;

        return rotated_dPdF;
    }

};

}

#endif

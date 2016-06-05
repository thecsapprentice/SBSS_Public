//#####################################################################
// Copyright 2011, Taylor Patterson, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class NEOHOOKEAN
//#####################################################################
#ifndef __NEOHOOKEAN__
#define __NEOHOOKEAN__

#include <PhysBAM_Tools/Matrices/DIAGONAL_MATRIX_3X3.h>

#include "ROTATED_STRESS_DERIVATIVE.h"

namespace PhysBAM{

template<class T,int d>
class NEOHOOKEAN
{
protected:
    typedef T T_SCALAR;
    enum {dim=d};
    enum {is_mooney_rivlin=0};

public:

    // Consider making protected
    static T M(const DIAGONAL_MATRIX<T,d>& Sigma)
    {return log(Sigma.Determinant());}

    static DIAGONAL_MATRIX<T,d> Q_hat(const DIAGONAL_MATRIX<T,d>& Sigma)
    {return Sigma.Inverse();}

public:

    static DIAGONAL_MATRIX<T,d> Deviatoric_P_hat(const DIAGONAL_MATRIX<T,d>& Sigma,const T mu)
    {return mu*(Sigma-Sigma.Inverse());}

    static ROTATED_STRESS_DERIVATIVE<T,2> Rotated_Stress_Derivative(const DIAGONAL_MATRIX<T,2>& Sigma,const T mu,const T kappa)
    {
        T J=Sigma.Determinant();

        T alpha11=mu;
        T alpha12=mu;
        T alpha22=mu;

        T beta_const=(mu-kappa*log(J));
        T beta11=(1/sqr(Sigma(1,1)))*beta_const;
        T beta12=(1/(Sigma(1,1)*Sigma(2,2)))*beta_const;
        T beta22=(1/sqr(Sigma(2,2)))*beta_const;

        T gamma11=kappa/sqr(Sigma(1,1));
        T gamma12=kappa/(Sigma(1,1)*Sigma(2,2));
        T gamma22=kappa/sqr(Sigma(2,2));

        ROTATED_STRESS_DERIVATIVE<T,2> rotated_dPdF;

        rotated_dPdF.a1111=alpha11+beta11+gamma11;
        rotated_dPdF.a1122=gamma12;
        rotated_dPdF.a2222=alpha22+beta22+gamma22;

        rotated_dPdF.a1212=alpha12;
        rotated_dPdF.a1221=beta12;

        return rotated_dPdF;
    }

    static ROTATED_STRESS_DERIVATIVE<T,3> Rotated_Stress_Derivative(const DIAGONAL_MATRIX<T,3>& Sigma,const T mu,const T kappa)
    {
        T J=Sigma.Determinant();

        T alpha11=mu;
        T alpha12=mu;
        T alpha13=mu;
        T alpha22=mu;
        T alpha23=mu;
        T alpha33=mu;

        T beta_const=(mu-kappa*log(J));
        T beta11=beta_const/sqr(Sigma(1,1));
        T beta12=beta_const/(Sigma(1,1)*Sigma(2,2));
        T beta13=beta_const/(Sigma(1,1)*Sigma(3,3));
        T beta22=beta_const/sqr(Sigma(2,2));
        T beta23=beta_const/(Sigma(2,2)*Sigma(3,3));
        T beta33=beta_const/sqr(Sigma(3,3));

        T gamma11=kappa/sqr(Sigma(1,1));
        T gamma12=kappa/(Sigma(1,1)*Sigma(2,2));
        T gamma13=kappa/(Sigma(1,1)*Sigma(3,3));
        T gamma22=kappa/sqr(Sigma(2,2));
        T gamma23=kappa/(Sigma(2,2)*Sigma(3,3));
        T gamma33=kappa/sqr(Sigma(3,3));
    
        ROTATED_STRESS_DERIVATIVE<T,3> rotated_dPdF;

        rotated_dPdF.a1111=alpha11+beta11+gamma11;
        rotated_dPdF.a1122=gamma12;
        rotated_dPdF.a1133=gamma13;
        rotated_dPdF.a2222=alpha22+beta22+gamma22;
        rotated_dPdF.a2233=gamma23;
        rotated_dPdF.a3333=alpha33+beta33+gamma33;

        rotated_dPdF.a1212=alpha12;
        rotated_dPdF.a1221=beta12;
        rotated_dPdF.a1313=alpha13;
        rotated_dPdF.a1331=beta13;
        rotated_dPdF.a2323=alpha23;
        rotated_dPdF.a2332=beta23;

        return rotated_dPdF;
    }

    static ROTATED_STRESS_DERIVATIVE<T,2> Augmented_Rotated_Stress_Derivative(const DIAGONAL_MATRIX<T,2>& Sigma,const T p,const T mu,const T alpha)
    {
        T alpha11=mu;
        T alpha12=mu;
        T alpha22=mu;

        T beta11=(mu-alpha*p)/sqr(Sigma(1,1));
        T beta12=(mu-alpha*p)/(Sigma(1,1)*Sigma(2,2));
        T beta22=(mu-alpha*p)/sqr(Sigma(2,2));
    
        ROTATED_STRESS_DERIVATIVE<T,2> rotated_dPdF;

        rotated_dPdF.a1111=alpha11+beta11;
        rotated_dPdF.a1122=T();
        rotated_dPdF.a2222=alpha22+beta22;

        rotated_dPdF.a1212=alpha12;
        rotated_dPdF.a1221=beta12;

        return rotated_dPdF;
    }

    static ROTATED_STRESS_DERIVATIVE<T,3> Augmented_Rotated_Stress_Derivative(const DIAGONAL_MATRIX<T,3>& Sigma,const T p,const T mu,const T alpha)
    {
        T alpha11=mu;
        T alpha12=mu;
        T alpha13=mu;
        T alpha22=mu;
        T alpha23=mu;
        T alpha33=mu;

        T beta11=(mu-alpha*p)/sqr(Sigma(1,1));
        T beta12=(mu-alpha*p)/(Sigma(1,1)*Sigma(2,2));
        T beta13=(mu-alpha*p)/(Sigma(1,1)*Sigma(3,3));
        T beta22=(mu-alpha*p)/sqr(Sigma(2,2));
        T beta23=(mu-alpha*p)/(Sigma(2,2)*Sigma(3,3));
        T beta33=(mu-alpha*p)/sqr(Sigma(3,3));
    
        ROTATED_STRESS_DERIVATIVE<T,3> rotated_dPdF;

        rotated_dPdF.a1111=alpha11+beta11;
        rotated_dPdF.a1122=T();
        rotated_dPdF.a1133=T();
        rotated_dPdF.a2222=alpha22+beta22;
        rotated_dPdF.a2233=T();
        rotated_dPdF.a3333=alpha33+beta33;

        rotated_dPdF.a1212=alpha12;
        rotated_dPdF.a1221=beta12;
        rotated_dPdF.a1313=alpha13;
        rotated_dPdF.a1331=beta13;
        rotated_dPdF.a2323=alpha23;
        rotated_dPdF.a2332=beta23;

        return rotated_dPdF;
    }

//#################################################
    static DIAGONAL_MATRIX<T,d> Deviatoric_P_hat(const DIAGONAL_MATRIX<T,d>& Sigma,const T mu_10,const T mu_01)
    {PHYSBAM_FATAL_ERROR();}
    static ROTATED_STRESS_DERIVATIVE<T,d> Rotated_Stress_Derivative(const DIAGONAL_MATRIX<T,d>& Sigma,const T mu_10,const T mu_01,const T kappa) {PHYSBAM_FATAL_ERROR();}
    static ROTATED_STRESS_DERIVATIVE<T,d> Augmented_Rotated_Stress_Derivative(const DIAGONAL_MATRIX<T,d>& Sigma,const T p,const T mu_10,const T mu_01,const T alpha) {PHYSBAM_FATAL_ERROR();}
//#################################################

};

}

#endif

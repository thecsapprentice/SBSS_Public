//#####################################################################
// Copyright 2011, Taylor Patterson, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ST_VENANT_KIRCHHOFF
//#####################################################################
#ifndef __ST_VENANT_KIRCHHOFF__
#define __ST_VENANT_KIRCHHOFF__

#include "MATERIAL_MODEL.h"

namespace PhysBAM{

template<class T,int d>
class ST_VENANT_KIRCHHOFF
{
public:

    typedef T SCALAR;
    enum {dim=d};

protected:

    enum {is_mooney_rivlin=0};

public:

    // Consider making protected
    static T M(const DIAGONAL_MATRIX<T,d>& Sigma)
    {return (sqr(Sigma)-1.).Trace()*0.5;}

    static DIAGONAL_MATRIX<T,d> Q_hat(const DIAGONAL_MATRIX<T,d>& Sigma)
    {return Sigma;}

protected:

    static DIAGONAL_MATRIX<T,d> Deviatoric_P_hat(const DIAGONAL_MATRIX<T,d>& Sigma,const T mu)
    {
        DIAGONAL_MATRIX<T,d> strain=(sqr(Sigma)-1.)*.5;
        DIAGONAL_MATRIX<T,d> P_hat=Sigma*(strain*2.*mu);
        return P_hat;
    }

    static ROTATED_STRESS_DERIVATIVE<T,2> Rotated_Stress_Derivative(const DIAGONAL_MATRIX<T,2>& Sigma,const T mu,const T kappa)
    {
        T alpha_const=(T)0.5*kappa*Sigma.Transpose_Times(Sigma).Trace()-kappa-mu;
        T alpha11=mu*(sqr(Sigma(1,1))+sqr(Sigma(1,1)))+alpha_const;
        T alpha12=mu*(sqr(Sigma(1,1))+sqr(Sigma(2,2)))+alpha_const;
        T alpha22=mu*(sqr(Sigma(2,2))+sqr(Sigma(2,2)))+alpha_const;

        T beta11=mu*sqr(Sigma(1,1));
        T beta12=mu*(Sigma(1,1)*Sigma(2,2));
        T beta22=mu*sqr(Sigma(2,2));

        T gamma11=kappa*sqr(Sigma(1,1));
        T gamma12=kappa*(Sigma(1,1)*Sigma(2,2));
        T gamma22=kappa*sqr(Sigma(2,2));
    
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
        T alpha_const=(T)0.5*kappa*Sigma.Transpose_Times(Sigma).Trace()-(T)1.5*kappa-mu;
        T alpha11=mu*(sqr(Sigma(1,1))+sqr(Sigma(1,1)))+alpha_const;
        T alpha12=mu*(sqr(Sigma(1,1))+sqr(Sigma(2,2)))+alpha_const;
        T alpha13=mu*(sqr(Sigma(1,1))+sqr(Sigma(3,3)))+alpha_const;
        T alpha22=mu*(sqr(Sigma(2,2))+sqr(Sigma(2,2)))+alpha_const;
        T alpha23=mu*(sqr(Sigma(2,2))+sqr(Sigma(3,3)))+alpha_const;
        T alpha33=mu*(sqr(Sigma(3,3))+sqr(Sigma(3,3)))+alpha_const;

        T beta11=mu*sqr(Sigma(1,1));
        T beta12=mu*(Sigma(1,1)*Sigma(2,2));
        T beta13=mu*(Sigma(1,1)*Sigma(3,3));
        T beta22=mu*sqr(Sigma(2,2));
        T beta23=mu*(Sigma(2,2)*Sigma(3,3));
        T beta33=mu*sqr(Sigma(3,3));

        T gamma11=kappa*sqr(Sigma(1,1));
        T gamma12=kappa*(Sigma(1,1)*Sigma(2,2));
        T gamma13=kappa*(Sigma(1,1)*Sigma(3,3));
        T gamma22=kappa*sqr(Sigma(2,2));
        T gamma23=kappa*(Sigma(2,2)*Sigma(3,3));
        T gamma33=kappa*sqr(Sigma(3,3));
    
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
        T alpha_const=alpha*p-mu;
        T alpha11=mu*(sqr(Sigma(1,1))+sqr(Sigma(1,1)))+alpha_const;
        T alpha12=mu*(sqr(Sigma(1,1))+sqr(Sigma(2,2)))+alpha_const;
        T alpha22=mu*(sqr(Sigma(2,2))+sqr(Sigma(2,2)))+alpha_const;

        T beta11=mu*sqr(Sigma(1,1));
        T beta12=mu*(Sigma(1,1)*Sigma(2,2));
        T beta22=mu*sqr(Sigma(2,2));
    
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
        T alpha_const=alpha*p-mu;
        T alpha11=mu*(sqr(Sigma(1,1))+sqr(Sigma(1,1)))+alpha_const;
        T alpha12=mu*(sqr(Sigma(1,1))+sqr(Sigma(2,2)))+alpha_const;
        T alpha13=mu*(sqr(Sigma(1,1))+sqr(Sigma(3,3)))+alpha_const;
        T alpha22=mu*(sqr(Sigma(2,2))+sqr(Sigma(2,2)))+alpha_const;
        T alpha23=mu*(sqr(Sigma(2,2))+sqr(Sigma(3,3)))+alpha_const;
        T alpha33=mu*(sqr(Sigma(3,3))+sqr(Sigma(3,3)))+alpha_const;

        T beta11=mu*sqr(Sigma(1,1));
        T beta12=mu*(Sigma(1,1)*Sigma(2,2));
        T beta13=mu*(Sigma(1,1)*Sigma(3,3));
        T beta22=mu*sqr(Sigma(2,2));
        T beta23=mu*(Sigma(2,2)*Sigma(3,3));
        T beta33=mu*sqr(Sigma(3,3));
    
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

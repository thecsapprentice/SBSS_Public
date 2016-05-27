//#####################################################################
// Copyright 2011, Taylor Patterson, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MOONEY_RIVLIN_MATERIAL
//#####################################################################
#ifndef __MOONEY_RIVLIN_MATERIAL__
#define __MOONEY_RIVLIN_MATERIAL__

#include <PhysBAM_Tools/Matrices/MATRIX_2X3.h>
#include <PhysBAM_Tools/Matrices/MATRIX_3X2.h>

#include "ROTATED_STRESS_DERIVATIVE.h"

namespace PhysBAM{

template<class T,int d>
class MOONEY_RIVLIN_MATERIAL
{
protected:
    T alpha;
    T mu_10;
    T mu_01;
    T kappa;

    T mu;
    T lambda;

public:
    MOONEY_RIVLIN_MATERIAL() {}

    static bool Is_Mooney_Rivlin() {return true;}

    void Initialize_Parameters(const T alpha_input,const T mu_input,const T lambda_input)
    {PHYSBAM_FATAL_ERROR("Improper parameters");}

    void Initialize_Parameters(const T alpha_input,const T mu_10_input,const T mu_01_input,const T kappa_input)
    {
        alpha=alpha_input;mu_10=mu_10_input;mu_01=mu_01_input;kappa=kappa_input;
        mu=2*(mu_10+mu_01);
        lambda=(T)one_third*(8*mu_01-4*mu_10)+kappa;
    }

    DIAGONAL_MATRIX<T,d> P_hat(const DIAGONAL_MATRIX<T,d>& Sigma,const T p) const
    {
        DIAGONAL_MATRIX<T,d> C=Sigma*Sigma,F_cube=C*Sigma,F_inverse=Sigma.Inverse();
        T I_C=C.Trace(),II_C=(C*C).Trace(),J=Sigma.Determinant(),Jcc=pow(J,-(T)two_thirds);
        DIAGONAL_MATRIX<T,d> P_hat=(2*Jcc*(mu_10+Jcc*mu_01*I_C))*Sigma-(2*Jcc*Jcc*mu_01)*F_cube+((-(T)two_thirds*Jcc*(mu_10*I_C+Jcc*mu_01*(I_C*I_C-II_C))))*F_inverse;

        P_hat+=alpha*p*Sigma.Inverse();
        return P_hat;
    }

    ROTATED_STRESS_DERIVATIVE<T,d> Isotropic_Stress_Derivative(const DIAGONAL_MATRIX<T,d>& Sigma,const T p,const bool apply_definiteness_fix) const
    {
        ROTATED_STRESS_DERIVATIVE<T,d> rotated_dPdF=Augmented_Rotated_Stress_Derivative(Sigma,p,mu_10,mu_01,alpha);

#if 1
        ROTATED_STRESS_DERIVATIVE<T,d> rotated_dPdF_indefinite=Rotated_Stress_Derivative(Sigma,mu_10,mu_01,std::min((T)5*(mu_10+mu_01),kappa));
        ROTATED_STRESS_DERIVATIVE<T,d> rotated_dPdF_incremental_fix=rotated_dPdF_indefinite;
        rotated_dPdF_incremental_fix.Make_Positive_Definite();
        rotated_dPdF_incremental_fix-=rotated_dPdF_indefinite;
        if(apply_definiteness_fix) rotated_dPdF+=rotated_dPdF_incremental_fix;
#else
        if(apply_definiteness_fix) rotated_dPdF.Make_Positive_Definite();
#endif

        return rotated_dPdF;
    }

    MATRIX<T,d> dP_hat(const ROTATED_STRESS_DERIVATIVE<T,d>& rotated_dPdF,const DIAGONAL_MATRIX<T,d>& Sigma,const MATRIX<T,d>& dF_hat,const T dp) const
    {
        MATRIX<T,d> dP_hat=rotated_dPdF.dP_hat(dF_hat);
        dP_hat+=MATRIX<T,d>(alpha*dp*Sigma.Inverse());
        return dP_hat;
    }

    T q(const DIAGONAL_MATRIX<T,d>& Sigma,const T p) const
    {
        T J=Sigma.Determinant();
        return -alpha*log(J)+sqr(alpha)*p/kappa;
    }

    T dq(const DIAGONAL_MATRIX<T,d>& Sigma,const MATRIX<T,d>& dF_hat,const T dp) const
    {
        return -alpha*Sigma.Inverse().Times_Transpose(dF_hat).Trace()+sqr(alpha)*dp/kappa;
    }

    T M(const DIAGONAL_MATRIX<T,d>& Sigma) const
    {
        T J=Sigma.Determinant();
        return log(J);
    }

    DIAGONAL_MATRIX<T,d> Q_hat(const DIAGONAL_MATRIX<T,d>& Sigma) const
    {
        return Sigma.Inverse();
    }

private:


    static ROTATED_STRESS_DERIVATIVE<T,2> Rotated_Stress_Derivative(const DIAGONAL_MATRIX<T,2>& Sigma,const T mu_10,const T mu_01,const T kappa)
    {
        ROTATED_STRESS_DERIVATIVE<T,2> rotated_dPdF;

        DIAGONAL_MATRIX<T,2> C=Sigma*Sigma,F_cube=C*Sigma,F_inverse=Sigma.Inverse();
        T I_C=C.Trace(),II_C=(C*C).Trace(),J=Sigma.Determinant(),Jcc=pow(J,-(T)two_thirds);
        SYMMETRIC_MATRIX<T,2> alpha;
        alpha.x11=2*Jcc*(mu_10+Jcc*mu_01*(I_C-C.x11-C.x11));
        alpha.x21=2*Jcc*(mu_10+Jcc*mu_01*(I_C-C.x22-C.x11));
        alpha.x22=2*Jcc*(mu_10+Jcc*mu_01*(I_C-C.x22-C.x22));
        SYMMETRIC_MATRIX<T,2> beta;
        beta.x11=(Jcc*(T)two_thirds*(mu_10*I_C+Jcc*mu_01*(I_C*I_C-II_C))-kappa*log(J))*F_inverse.x11*F_inverse.x11-2*Jcc*Jcc*mu_01*Sigma.x11*Sigma.x11;
        beta.x21=(Jcc*(T)two_thirds*(mu_10*I_C+Jcc*mu_01*(I_C*I_C-II_C))-kappa*log(J))*F_inverse.x22*F_inverse.x11-2*Jcc*Jcc*mu_01*Sigma.x22*Sigma.x11;
        beta.x22=(Jcc*(T)two_thirds*(mu_10*I_C+Jcc*mu_01*(I_C*I_C-II_C))-kappa*log(J))*F_inverse.x22*F_inverse.x22-2*Jcc*Jcc*mu_01*Sigma.x22*Sigma.x22;
        SYMMETRIC_MATRIX<T,3> eta;
        eta.x11=4*Jcc*Jcc*mu_01;
        eta.x31=-Jcc*(T)one_third*(4*mu_10+8*Jcc*mu_01*I_C);
        eta.x32=8*(T)one_third*Jcc*Jcc*mu_01;
        eta.x33=Jcc*(T)one_ninth*(4*mu_10*I_C+Jcc*8*mu_01*(I_C*I_C-II_C))+kappa;
        MATRIX<T,2,3> F_base(Sigma.x11,Sigma.x22,F_cube.x11,F_cube.x22,F_inverse.x11,F_inverse.x22);
        SYMMETRIC_MATRIX<T,2> gamma=SYMMETRIC_MATRIX<T,2>::Conjugate(F_base,eta);
        rotated_dPdF.a1111=alpha.x11+beta.x11+gamma.x11;
        rotated_dPdF.a1122=gamma.x21;
        rotated_dPdF.a2222=alpha.x22+beta.x22+gamma.x22;
        rotated_dPdF.a1212=alpha.x21;
        rotated_dPdF.a1221=beta.x21;

        return rotated_dPdF;
    }

    static ROTATED_STRESS_DERIVATIVE<T,3> Rotated_Stress_Derivative(const DIAGONAL_MATRIX<T,3>& Sigma,const T mu_10,const T mu_01,const T kappa)
    {
        ROTATED_STRESS_DERIVATIVE<T,3> rotated_dPdF;

        DIAGONAL_MATRIX<T,3> C=Sigma*Sigma,F_cube=C*Sigma,F_inverse=Sigma.Inverse();
        T I_C=C.Trace(),II_C=(C*C).Trace(),J=Sigma.Determinant(),Jcc=pow(J,-(T)two_thirds);
        SYMMETRIC_MATRIX<T,3> alpha;
        alpha.x11=2*Jcc*(mu_10+Jcc*mu_01*(I_C-C.x11-C.x11));alpha.x21=2*Jcc*(mu_10+Jcc*mu_01*(I_C-C.x22-C.x11));alpha.x31=2*Jcc*(mu_10+Jcc*mu_01*(I_C-C.x33-C.x11));
        alpha.x22=2*Jcc*(mu_10+Jcc*mu_01*(I_C-C.x22-C.x22));alpha.x32=2*Jcc*(mu_10+Jcc*mu_01*(I_C-C.x33-C.x22));alpha.x33=2*Jcc*(mu_10+Jcc*mu_01*(I_C-C.x33-C.x33));
        SYMMETRIC_MATRIX<T,3> beta;
        beta.x11=(Jcc*(T)two_thirds*(mu_10*I_C+Jcc*mu_01*(I_C*I_C-II_C))-kappa*log(J))*F_inverse.x11*F_inverse.x11-2*Jcc*Jcc*mu_01*Sigma.x11*Sigma.x11;
        beta.x21=(Jcc*(T)two_thirds*(mu_10*I_C+Jcc*mu_01*(I_C*I_C-II_C))-kappa*log(J))*F_inverse.x22*F_inverse.x11-2*Jcc*Jcc*mu_01*Sigma.x22*Sigma.x11;
        beta.x31=(Jcc*(T)two_thirds*(mu_10*I_C+Jcc*mu_01*(I_C*I_C-II_C))-kappa*log(J))*F_inverse.x33*F_inverse.x11-2*Jcc*Jcc*mu_01*Sigma.x33*Sigma.x11;
        beta.x22=(Jcc*(T)two_thirds*(mu_10*I_C+Jcc*mu_01*(I_C*I_C-II_C))-kappa*log(J))*F_inverse.x22*F_inverse.x22-2*Jcc*Jcc*mu_01*Sigma.x22*Sigma.x22;
        beta.x32=(Jcc*(T)two_thirds*(mu_10*I_C+Jcc*mu_01*(I_C*I_C-II_C))-kappa*log(J))*F_inverse.x33*F_inverse.x22-2*Jcc*Jcc*mu_01*Sigma.x33*Sigma.x22;
        beta.x33=(Jcc*(T)two_thirds*(mu_10*I_C+Jcc*mu_01*(I_C*I_C-II_C))-kappa*log(J))*F_inverse.x33*F_inverse.x33-2*Jcc*Jcc*mu_01*Sigma.x33*Sigma.x33;
        SYMMETRIC_MATRIX<T,3> eta;
        eta.x11=4*Jcc*Jcc*mu_01;
        eta.x31=-Jcc*(T)one_third*(4*mu_10+8*Jcc*mu_01*I_C);
        eta.x32=8*(T)one_third*Jcc*Jcc*mu_01;
        eta.x33=Jcc*(T)one_ninth*(4*mu_10*I_C+Jcc*8*mu_01*(I_C*I_C-II_C))+kappa;
        MATRIX<T,3> F_base(Sigma.x11,Sigma.x22,Sigma.x33,F_cube.x11,F_cube.x22,F_cube.x33,F_inverse.x11,F_inverse.x22,F_inverse.x33);
        SYMMETRIC_MATRIX<T,3> gamma=SYMMETRIC_MATRIX<T,3>::Conjugate(F_base,eta);
        rotated_dPdF.a1111=alpha.x11+beta.x11+gamma.x11;rotated_dPdF.a2222=alpha.x22+beta.x22+gamma.x22;rotated_dPdF.a3333=alpha.x33+beta.x33+gamma.x33;
        rotated_dPdF.a1122=gamma.x21;rotated_dPdF.a1133=gamma.x31;rotated_dPdF.a2233=gamma.x32;
        rotated_dPdF.a1212=alpha.x21;rotated_dPdF.a1313=alpha.x31;rotated_dPdF.a2323=alpha.x32;
        rotated_dPdF.a1221=beta.x21;rotated_dPdF.a1331=beta.x31;rotated_dPdF.a2332=beta.x31;

        return rotated_dPdF;
    }

    static ROTATED_STRESS_DERIVATIVE<T,2> Augmented_Rotated_Stress_Derivative(const DIAGONAL_MATRIX<T,2>& Sigma,const T p,const T mu_10,const T mu_01,const T alpha)
    {
        ROTATED_STRESS_DERIVATIVE<T,2> rotated_dPdF;
        rotated_dPdF=Rotated_Stress_Derivative(Sigma,mu_10,mu_01,T());

        T beta11=-alpha*p/sqr(Sigma(1,1));
        T beta12=-alpha*p/(Sigma(1,1)*Sigma(2,2));
        T beta22=-alpha*p/sqr(Sigma(2,2));    

        rotated_dPdF.a1111+=beta11;
        rotated_dPdF.a2222+=beta22;
        rotated_dPdF.a1221+=beta12;

        return rotated_dPdF;
    }

    static ROTATED_STRESS_DERIVATIVE<T,3> Augmented_Rotated_Stress_Derivative(const DIAGONAL_MATRIX<T,3>& Sigma,const T p,const T mu_10,const T mu_01,const T alpha)
    {
        ROTATED_STRESS_DERIVATIVE<T,3> rotated_dPdF;
        rotated_dPdF=Rotated_Stress_Derivative(Sigma,mu_10,mu_01,T());

        T beta11=-alpha*p/sqr(Sigma(1,1));
        T beta12=-alpha*p/(Sigma(1,1)*Sigma(2,2));
        T beta13=-alpha*p/(Sigma(1,1)*Sigma(3,3));
        T beta22=-alpha*p/sqr(Sigma(2,2));
        T beta23=-alpha*p/(Sigma(2,2)*Sigma(3,3));
        T beta33=-alpha*p/sqr(Sigma(3,3));

        rotated_dPdF.a1111+=beta11;
        rotated_dPdF.a2222+=beta22;
        rotated_dPdF.a3333+=beta33;

        rotated_dPdF.a1221+=beta12;
        rotated_dPdF.a1331+=beta13;
        rotated_dPdF.a2332+=beta23;

        return rotated_dPdF;
    }

};

}

#endif

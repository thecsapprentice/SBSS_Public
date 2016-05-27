//#####################################################################
// Copyright 2011, Taylor Patterson, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MOONEY_RIVLIN
//#####################################################################
#ifndef __MOONEY_RIVLIN__
#define __MOONEY_RIVLIN__

#include <PhysBAM_Tools/Math_Tools/cube.h>

#include "MATERIAL_MODEL.h"

namespace PhysBAM{

template<class T,int d>
class MOONEY_RIVLIN
{
public:

    typedef T SCALAR;
    enum {dim=d};

protected:

    enum {is_mooney_rivlin=1};

public:

    // Consider making protected
    static T M(const DIAGONAL_MATRIX<T,d>& Sigma)
    {return log(Sigma.Determinant());}

    static DIAGONAL_MATRIX<T,d> Q_hat(const DIAGONAL_MATRIX<T,d>& Sigma)
    {return Sigma.Inverse();}

protected:

    static T Deviatoric_Energy(const DIAGONAL_MATRIX<T,d>& Sigma,const T mu_10,const T mu_01)
    {
        DIAGONAL_MATRIX<T,d> C=Sigma*Sigma;
        T I_C=C.Trace(),II_C=(C*C).Trace(),III_C=C.Determinant();
	T I_C_hat=pow(III_C,-1./(T)d)*I_C,II_C_hat=pow(III_C,-2./(T)d)*II_C;

        return mu_10*(I_C_hat-d)+.5*mu_01*(sqr(I_C_hat)-II_C_hat-sqr(d)+d);
    }

    static DIAGONAL_MATRIX<T,d> Deviatoric_P_hat(const DIAGONAL_MATRIX<T,d>& Sigma,const T mu_10,const T mu_01)
    {
        DIAGONAL_MATRIX<T,d> C=Sigma*Sigma,Sigma_cube=C*Sigma,Sigma_inverse=Sigma.Inverse();
        T I_C=C.Trace(),II_C=(C*C).Trace(),J=Sigma.Determinant(),III_C=sqr(J);

    	T Psi_I=(mu_10*pow(III_C,-1./(T)d))+(mu_01*I_C*pow(III_C,-2./(T)d));
    	T Psi_II=-.5*mu_01*pow(III_C,-2./(T)d);
    	T Psi_III=-(mu_10/(T)d)*I_C*pow(III_C,-(T)(d+1)/(T)d)-(mu_01/(T)d)*(sqr(I_C)-II_C)*pow(III_C,-(T)(d+2)/(T)d);

        DIAGONAL_MATRIX<T,d> P_hat=2*Psi_I*Sigma+4*Psi_II*Sigma_cube+2*III_C*Psi_III*Sigma_inverse;

        return P_hat;
    }

    static ROTATED_STRESS_DERIVATIVE<T,2> Rotated_Stress_Derivative(const DIAGONAL_MATRIX<T,2>& Sigma,const T mu_10,const T mu_01,const T kappa)
    {
        ROTATED_STRESS_DERIVATIVE<T,2> rotated_dPdF;
        
        DIAGONAL_MATRIX<T,d> C=Sigma*Sigma,Sigma_cube=C*Sigma,Sigma_inverse=Sigma.Inverse();
        T I_C=C.Trace(),II_C=(C*C).Trace(),III_C=C.Determinant();
        
	T III_C_inverse=1./III_C,III_C_root_minus_d=pow(III_C_inverse,1./(T)d);
        
    	T Psi_I=mu_10*III_C_root_minus_d+mu_01*I_C*sqr(III_C_root_minus_d);
    	T Psi_II=-.5*mu_01*sqr(III_C_root_minus_d);
	T Psi_III=-(mu_10/(T)d)*I_C*III_C_inverse*III_C_root_minus_d-(mu_01/(T)d)*(sqr(I_C)-II_C)*III_C_inverse*sqr(III_C_root_minus_d)+.25*kappa*log(III_C)/III_C;
        
        SYMMETRIC_MATRIX<T,3> hessian;
        hessian.x11=mu_01*sqr(III_C_root_minus_d);
        hessian.x21=0;
        hessian.x31=-(mu_10/(T)d)*III_C_inverse*III_C_root_minus_d-(2*mu_01/(T)d)*I_C*III_C_inverse*sqr(III_C_root_minus_d);
        hessian.x22=0;
        hessian.x32=(mu_01/(T)d)*III_C_inverse*sqr(III_C_root_minus_d);
        hessian.x33=(mu_10*(T)(d+1)/sqr((T)d))*I_C*sqr(III_C_inverse)*III_C_root_minus_d+(mu_01*(T)(d+2)/sqr((T)d))*(sqr(I_C)-II_C)*sqr(III_C_inverse)*sqr(III_C_root_minus_d)+.25*kappa*(1.-log(III_C))/sqr(III_C);

        SYMMETRIC_MATRIX<T,d> Alpha;
        for(int i=1;i<=d;i++)
            for(int j=1;j<=i;j++)
                Alpha(i,j)=2*Psi_I+4*Psi_II*(C(i,i)+C(j,j));
        
        SYMMETRIC_MATRIX<T,d> Beta;
        for(int i=1;i<=d;i++)
            for(int j=1;j<=i;j++)
                Beta(i,j)=4*Psi_II*Sigma(i,i)*Sigma(j,j)-(2*III_C*Psi_III)*Sigma_inverse(i,i)*Sigma_inverse(j,j);

        SYMMETRIC_MATRIX<T,d> Gamma;
        Gamma(1,1)=4.*(hessian.x11*Sigma(1,1)*Sigma(1,1)+hessian.x31*III_C*(Sigma(1,1)*Sigma_inverse(1,1)+Sigma(1,1)*Sigma_inverse(1,1))
		       +2.*hessian.x32*III_C*(Sigma_cube(1,1)*Sigma_inverse(1,1)+Sigma_cube(1,1)*Sigma_inverse(1,1))
		       +(hessian.x33*III_C+Psi_III)*III_C*Sigma_inverse(1,1)*Sigma_inverse(1,1));
        Gamma(2,1)=4.*(hessian.x11*Sigma(2,2)*Sigma(1,1)+hessian.x31*III_C*(Sigma(1,1)*Sigma_inverse(2,2)+Sigma(2,2)*Sigma_inverse(1,1))
		       +2.*hessian.x32*III_C*(Sigma_cube(1,1)*Sigma_inverse(2,2)+Sigma_cube(2,2)*Sigma_inverse(1,1))
		       +(hessian.x33*III_C+Psi_III)*III_C*Sigma_inverse(2,2)*Sigma_inverse(1,1));
        Gamma(2,2)=4.*(hessian.x11*Sigma(2,2)*Sigma(2,2)+hessian.x31*III_C*(Sigma(2,2)*Sigma_inverse(2,2)+Sigma(2,2)*Sigma_inverse(2,2))
		       +2.*hessian.x32*III_C*(Sigma_cube(2,2)*Sigma_inverse(2,2)+Sigma_cube(2,2)*Sigma_inverse(2,2))
		       +(hessian.x33*III_C+Psi_III)*III_C*Sigma_inverse(2,2)*Sigma_inverse(2,2));

        rotated_dPdF.a1111=Alpha.x11+Beta.x11+Gamma.x11;
    	rotated_dPdF.a2222=Alpha.x22+Beta.x22+Gamma.x22;
        rotated_dPdF.a1122=Gamma.x21;
        rotated_dPdF.a1212=Alpha.x21;
        rotated_dPdF.a1221=Beta.x21;

	return rotated_dPdF;
    }

    static ROTATED_STRESS_DERIVATIVE<T,3> Rotated_Stress_Derivative(const DIAGONAL_MATRIX<T,3>& Sigma,const T mu_10,const T mu_01,const T kappa)
    {
        ROTATED_STRESS_DERIVATIVE<T,3> rotated_dPdF;
        
        DIAGONAL_MATRIX<T,d> C=Sigma*Sigma,Sigma_cube=C*Sigma,Sigma_inverse=Sigma.Inverse();
        T I_C=C.Trace(),II_C=(C*C).Trace(),III_C=C.Determinant();
        
	T III_C_inverse=1./III_C,III_C_root_minus_d=pow(III_C_inverse,1./(T)d);
        
    	T Psi_I=mu_10*III_C_root_minus_d+mu_01*I_C*sqr(III_C_root_minus_d);
    	T Psi_II=-.5*mu_01*sqr(III_C_root_minus_d);
	T Psi_III=-(mu_10/(T)d)*I_C*III_C_inverse*III_C_root_minus_d-(mu_01/(T)d)*(sqr(I_C)-II_C)*III_C_inverse*sqr(III_C_root_minus_d)+.25*kappa*log(III_C)/III_C;
        
        SYMMETRIC_MATRIX<T,3> hessian;
        hessian.x11=mu_01*sqr(III_C_root_minus_d);
        hessian.x21=0;
        hessian.x31=-(mu_10/(T)d)*III_C_inverse*III_C_root_minus_d-(2*mu_01/(T)d)*I_C*III_C_inverse*sqr(III_C_root_minus_d);
        hessian.x22=0;
        hessian.x32=(mu_01/(T)d)*III_C_inverse*sqr(III_C_root_minus_d);
        hessian.x33=(mu_10*(T)(d+1)/sqr((T)d))*I_C*sqr(III_C_inverse)*III_C_root_minus_d+(mu_01*(T)(d+2)/sqr((T)d))*(sqr(I_C)-II_C)*sqr(III_C_inverse)*sqr(III_C_root_minus_d)+.25*kappa*(1.-log(III_C))/sqr(III_C);

        SYMMETRIC_MATRIX<T,d> Alpha;
        for(int i=1;i<=d;i++)
            for(int j=1;j<=i;j++)
                Alpha(i,j)=2*Psi_I+4*Psi_II*(C(i,i)+C(j,j));
        
        SYMMETRIC_MATRIX<T,d> Beta;
        for(int i=1;i<=d;i++)
            for(int j=1;j<=i;j++)
                Beta(i,j)=4*Psi_II*Sigma(i,i)*Sigma(j,j)-(2*III_C*Psi_III)*Sigma_inverse(i,i)*Sigma_inverse(j,j);

        SYMMETRIC_MATRIX<T,d> Gamma;
        Gamma(1,1)=4.*(hessian.x11*Sigma(1,1)*Sigma(1,1)+hessian.x31*III_C*(Sigma(1,1)*Sigma_inverse(1,1)+Sigma(1,1)*Sigma_inverse(1,1))
		       +2.*hessian.x32*III_C*(Sigma_cube(1,1)*Sigma_inverse(1,1)+Sigma_cube(1,1)*Sigma_inverse(1,1))
		       +(hessian.x33*III_C+Psi_III)*III_C*Sigma_inverse(1,1)*Sigma_inverse(1,1));
        Gamma(2,1)=4.*(hessian.x11*Sigma(2,2)*Sigma(1,1)+hessian.x31*III_C*(Sigma(1,1)*Sigma_inverse(2,2)+Sigma(2,2)*Sigma_inverse(1,1))
		       +2.*hessian.x32*III_C*(Sigma_cube(1,1)*Sigma_inverse(2,2)+Sigma_cube(2,2)*Sigma_inverse(1,1))
		       +(hessian.x33*III_C+Psi_III)*III_C*Sigma_inverse(2,2)*Sigma_inverse(1,1));
        Gamma(2,2)=4.*(hessian.x11*Sigma(2,2)*Sigma(2,2)+hessian.x31*III_C*(Sigma(2,2)*Sigma_inverse(2,2)+Sigma(2,2)*Sigma_inverse(2,2))
		       +2.*hessian.x32*III_C*(Sigma_cube(2,2)*Sigma_inverse(2,2)+Sigma_cube(2,2)*Sigma_inverse(2,2))
		       +(hessian.x33*III_C+Psi_III)*III_C*Sigma_inverse(2,2)*Sigma_inverse(2,2));
        Gamma(3,1)=4.*(hessian.x11*Sigma(3,3)*Sigma(1,1)+hessian.x31*III_C*(Sigma(1,1)*Sigma_inverse(3,3)+Sigma(3,3)*Sigma_inverse(1,1))
		       +2.*hessian.x32*III_C*(Sigma_cube(1,1)*Sigma_inverse(3,3)+Sigma_cube(3,3)*Sigma_inverse(1,1))
		       +(hessian.x33*III_C+Psi_III)*III_C*Sigma_inverse(3,3)*Sigma_inverse(1,1));
        Gamma(3,2)=4.*(hessian.x11*Sigma(3,3)*Sigma(2,2)+hessian.x31*III_C*(Sigma(2,2)*Sigma_inverse(3,3)+Sigma(3,3)*Sigma_inverse(2,2))
		       +2.*hessian.x32*III_C*(Sigma_cube(2,2)*Sigma_inverse(3,3)+Sigma_cube(3,3)*Sigma_inverse(2,2))
		       +(hessian.x33*III_C+Psi_III)*III_C*Sigma_inverse(3,3)*Sigma_inverse(2,2));
        Gamma(3,3)=4.*(hessian.x11*Sigma(3,3)*Sigma(3,3)+hessian.x31*III_C*(Sigma(3,3)*Sigma_inverse(3,3)+Sigma(3,3)*Sigma_inverse(3,3))
		       +2.*hessian.x32*III_C*(Sigma_cube(3,3)*Sigma_inverse(3,3)+Sigma_cube(3,3)*Sigma_inverse(3,3))
		       +(hessian.x33*III_C+Psi_III)*III_C*Sigma_inverse(3,3)*Sigma_inverse(3,3));

        rotated_dPdF.a1111=Alpha.x11+Beta.x11+Gamma.x11;
    	rotated_dPdF.a2222=Alpha.x22+Beta.x22+Gamma.x22;
    	rotated_dPdF.a3333=Alpha.x33+Beta.x33+Gamma.x33;
        rotated_dPdF.a1122=Gamma.x21;
    	rotated_dPdF.a1133=Gamma.x31;
    	rotated_dPdF.a2233=Gamma.x32;
        rotated_dPdF.a1212=Alpha.x21;
    	rotated_dPdF.a1313=Alpha.x31;
    	rotated_dPdF.a2323=Alpha.x32;
        rotated_dPdF.a1221=Beta.x21;
    	rotated_dPdF.a1331=Beta.x31;
    	rotated_dPdF.a2332=Beta.x32;

	return rotated_dPdF;
    }

    static ROTATED_STRESS_DERIVATIVE<T,2> Augmented_Rotated_Stress_Derivative(const DIAGONAL_MATRIX<T,2>& Sigma,const T p,const T mu_10,const T mu_01,const T alpha)
    {
        ROTATED_STRESS_DERIVATIVE<T,2> rotated_dPdF;
        
        DIAGONAL_MATRIX<T,d> C=Sigma*Sigma,Sigma_cube=C*Sigma,Sigma_inverse=Sigma.Inverse();
        T I_C=C.Trace(),II_C=(C*C).Trace(),III_C=C.Determinant();

	T III_C_inverse=1./III_C,III_C_root_minus_d=pow(III_C_inverse,1./(T)d);
        
    	T Psi_I=mu_10*III_C_root_minus_d+mu_01*I_C*sqr(III_C_root_minus_d);
    	T Psi_II=-.5*mu_01*sqr(III_C_root_minus_d);
	T Psi_III=-(mu_10/(T)d)*I_C*III_C_inverse*III_C_root_minus_d-(mu_01/(T)d)*(sqr(I_C)-II_C)*III_C_inverse*sqr(III_C_root_minus_d)+.5*alpha*p/III_C;
        
        SYMMETRIC_MATRIX<T,3> hessian;
        hessian.x11=mu_01*sqr(III_C_root_minus_d);
        hessian.x21=0;
        hessian.x31=-(mu_10/(T)d)*III_C_inverse*III_C_root_minus_d-(2*mu_01/(T)d)*I_C*III_C_inverse*sqr(III_C_root_minus_d);
        hessian.x22=0;
        hessian.x32=(mu_01/(T)d)*III_C_inverse*sqr(III_C_root_minus_d);
        hessian.x33=(mu_10*(T)(d+1)/sqr((T)d))*I_C*sqr(III_C_inverse)*III_C_root_minus_d+(mu_01*(T)(d+2)/sqr((T)d))*(sqr(I_C)-II_C)*sqr(III_C_inverse)*sqr(III_C_root_minus_d)-.5*alpha*p/sqr(III_C);

        SYMMETRIC_MATRIX<T,d> Alpha;
        for(int i=1;i<=d;i++)
            for(int j=1;j<=i;j++)
                Alpha(i,j)=2*Psi_I+4*Psi_II*(C(i,i)+C(j,j));
        
        SYMMETRIC_MATRIX<T,d> Beta;
        for(int i=1;i<=d;i++)
            for(int j=1;j<=i;j++)
                Beta(i,j)=4*Psi_II*Sigma(i,i)*Sigma(j,j)-(2*III_C*Psi_III)*Sigma_inverse(i,i)*Sigma_inverse(j,j);

        SYMMETRIC_MATRIX<T,d> Gamma;
        Gamma(1,1)=4.*(hessian.x11*Sigma(1,1)*Sigma(1,1)+hessian.x31*III_C*(Sigma(1,1)*Sigma_inverse(1,1)+Sigma(1,1)*Sigma_inverse(1,1))
		       +2.*hessian.x32*III_C*(Sigma_cube(1,1)*Sigma_inverse(1,1)+Sigma_cube(1,1)*Sigma_inverse(1,1))
		       +(hessian.x33*III_C+Psi_III)*III_C*Sigma_inverse(1,1)*Sigma_inverse(1,1));
        Gamma(2,1)=4.*(hessian.x11*Sigma(2,2)*Sigma(1,1)+hessian.x31*III_C*(Sigma(1,1)*Sigma_inverse(2,2)+Sigma(2,2)*Sigma_inverse(1,1))
		       +2.*hessian.x32*III_C*(Sigma_cube(1,1)*Sigma_inverse(2,2)+Sigma_cube(2,2)*Sigma_inverse(1,1))
		       +(hessian.x33*III_C+Psi_III)*III_C*Sigma_inverse(2,2)*Sigma_inverse(1,1));
        Gamma(2,2)=4.*(hessian.x11*Sigma(2,2)*Sigma(2,2)+hessian.x31*III_C*(Sigma(2,2)*Sigma_inverse(2,2)+Sigma(2,2)*Sigma_inverse(2,2))
		       +2.*hessian.x32*III_C*(Sigma_cube(2,2)*Sigma_inverse(2,2)+Sigma_cube(2,2)*Sigma_inverse(2,2))
		       +(hessian.x33*III_C+Psi_III)*III_C*Sigma_inverse(2,2)*Sigma_inverse(2,2));

        rotated_dPdF.a1111=Alpha.x11+Beta.x11+Gamma.x11;
    	rotated_dPdF.a2222=Alpha.x22+Beta.x22+Gamma.x22;
        rotated_dPdF.a1122=Gamma.x21;
        rotated_dPdF.a1212=Alpha.x21;
        rotated_dPdF.a1221=Beta.x21;

	return rotated_dPdF;
    }

    static ROTATED_STRESS_DERIVATIVE<T,3> Augmented_Rotated_Stress_Derivative(const DIAGONAL_MATRIX<T,3>& Sigma,const T p,const T mu_10,const T mu_01,const T alpha)
    {
        ROTATED_STRESS_DERIVATIVE<T,3> rotated_dPdF;
        
        DIAGONAL_MATRIX<T,d> C=Sigma*Sigma,Sigma_cube=C*Sigma,Sigma_inverse=Sigma.Inverse();
        T I_C=C.Trace(),II_C=(C*C).Trace(),III_C=C.Determinant();
        
	T III_C_inverse=1./III_C,III_C_root_minus_d=pow(III_C_inverse,1./(T)d);
        
    	T Psi_I=mu_10*III_C_root_minus_d+mu_01*I_C*sqr(III_C_root_minus_d);
    	T Psi_II=-.5*mu_01*sqr(III_C_root_minus_d);
	T Psi_III=-(mu_10/(T)d)*I_C*III_C_inverse*III_C_root_minus_d-(mu_01/(T)d)*(sqr(I_C)-II_C)*III_C_inverse*sqr(III_C_root_minus_d)+.5*alpha*p/III_C;
        
        SYMMETRIC_MATRIX<T,3> hessian;
        hessian.x11=mu_01*sqr(III_C_root_minus_d);
        hessian.x21=0;
        hessian.x31=-(mu_10/(T)d)*III_C_inverse*III_C_root_minus_d-(2*mu_01/(T)d)*I_C*III_C_inverse*sqr(III_C_root_minus_d);
        hessian.x22=0;
        hessian.x32=(mu_01/(T)d)*III_C_inverse*sqr(III_C_root_minus_d);
        hessian.x33=(mu_10*(T)(d+1)/sqr((T)d))*I_C*sqr(III_C_inverse)*III_C_root_minus_d+(mu_01*(T)(d+2)/sqr((T)d))*(sqr(I_C)-II_C)*sqr(III_C_inverse)*sqr(III_C_root_minus_d)-.5*alpha*p/sqr(III_C);

        SYMMETRIC_MATRIX<T,d> Alpha;
        for(int i=1;i<=d;i++)
            for(int j=1;j<=i;j++)
                Alpha(i,j)=2*Psi_I+4*Psi_II*(C(i,i)+C(j,j));
        
        SYMMETRIC_MATRIX<T,d> Beta;
        for(int i=1;i<=d;i++)
            for(int j=1;j<=i;j++)
                Beta(i,j)=4*Psi_II*Sigma(i,i)*Sigma(j,j)-(2*III_C*Psi_III)*Sigma_inverse(i,i)*Sigma_inverse(j,j);

        SYMMETRIC_MATRIX<T,d> Gamma;
        Gamma(1,1)=4.*(hessian.x11*Sigma(1,1)*Sigma(1,1)+hessian.x31*III_C*(Sigma(1,1)*Sigma_inverse(1,1)+Sigma(1,1)*Sigma_inverse(1,1))
		       +2.*hessian.x32*III_C*(Sigma_cube(1,1)*Sigma_inverse(1,1)+Sigma_cube(1,1)*Sigma_inverse(1,1))
		       +(hessian.x33*III_C+Psi_III)*III_C*Sigma_inverse(1,1)*Sigma_inverse(1,1));
        Gamma(2,1)=4.*(hessian.x11*Sigma(2,2)*Sigma(1,1)+hessian.x31*III_C*(Sigma(1,1)*Sigma_inverse(2,2)+Sigma(2,2)*Sigma_inverse(1,1))
		       +2.*hessian.x32*III_C*(Sigma_cube(1,1)*Sigma_inverse(2,2)+Sigma_cube(2,2)*Sigma_inverse(1,1))
		       +(hessian.x33*III_C+Psi_III)*III_C*Sigma_inverse(2,2)*Sigma_inverse(1,1));
        Gamma(2,2)=4.*(hessian.x11*Sigma(2,2)*Sigma(2,2)+hessian.x31*III_C*(Sigma(2,2)*Sigma_inverse(2,2)+Sigma(2,2)*Sigma_inverse(2,2))
		       +2.*hessian.x32*III_C*(Sigma_cube(2,2)*Sigma_inverse(2,2)+Sigma_cube(2,2)*Sigma_inverse(2,2))
		       +(hessian.x33*III_C+Psi_III)*III_C*Sigma_inverse(2,2)*Sigma_inverse(2,2));
        Gamma(3,1)=4.*(hessian.x11*Sigma(3,3)*Sigma(1,1)+hessian.x31*III_C*(Sigma(1,1)*Sigma_inverse(3,3)+Sigma(3,3)*Sigma_inverse(1,1))
		       +2.*hessian.x32*III_C*(Sigma_cube(1,1)*Sigma_inverse(3,3)+Sigma_cube(3,3)*Sigma_inverse(1,1))
		       +(hessian.x33*III_C+Psi_III)*III_C*Sigma_inverse(3,3)*Sigma_inverse(1,1));
        Gamma(3,2)=4.*(hessian.x11*Sigma(3,3)*Sigma(2,2)+hessian.x31*III_C*(Sigma(2,2)*Sigma_inverse(3,3)+Sigma(3,3)*Sigma_inverse(2,2))
		       +2.*hessian.x32*III_C*(Sigma_cube(2,2)*Sigma_inverse(3,3)+Sigma_cube(3,3)*Sigma_inverse(2,2))
		       +(hessian.x33*III_C+Psi_III)*III_C*Sigma_inverse(3,3)*Sigma_inverse(2,2));
        Gamma(3,3)=4.*(hessian.x11*Sigma(3,3)*Sigma(3,3)+hessian.x31*III_C*(Sigma(3,3)*Sigma_inverse(3,3)+Sigma(3,3)*Sigma_inverse(3,3))
		       +2.*hessian.x32*III_C*(Sigma_cube(3,3)*Sigma_inverse(3,3)+Sigma_cube(3,3)*Sigma_inverse(3,3))
		       +(hessian.x33*III_C+Psi_III)*III_C*Sigma_inverse(3,3)*Sigma_inverse(3,3));

        rotated_dPdF.a1111=Alpha.x11+Beta.x11+Gamma.x11;
    	rotated_dPdF.a2222=Alpha.x22+Beta.x22+Gamma.x22;
    	rotated_dPdF.a3333=Alpha.x33+Beta.x33+Gamma.x33;
        rotated_dPdF.a1122=Gamma.x21;
    	rotated_dPdF.a1133=Gamma.x31;
    	rotated_dPdF.a2233=Gamma.x32;
        rotated_dPdF.a1212=Alpha.x21;
    	rotated_dPdF.a1313=Alpha.x31;
    	rotated_dPdF.a2323=Alpha.x32;
        rotated_dPdF.a1221=Beta.x21;
    	rotated_dPdF.a1331=Beta.x31;
    	rotated_dPdF.a2332=Beta.x32;

	return rotated_dPdF;
    }

//#################################################
    static T Deviatoric_Energy(const DIAGONAL_MATRIX<T,d>& Sigma,const T mu) {PHYSBAM_FATAL_ERROR();}
    static DIAGONAL_MATRIX<T,d> Deviatoric_P_hat(const DIAGONAL_MATRIX<T,d>& Sigma,const T mu) {PHYSBAM_FATAL_ERROR();}
    static ROTATED_STRESS_DERIVATIVE<T,d> Rotated_Stress_Derivative(const DIAGONAL_MATRIX<T,d>& Sigma,const T mu,const T kappa) {PHYSBAM_FATAL_ERROR();}
    static ROTATED_STRESS_DERIVATIVE<T,d> Augmented_Rotated_Stress_Derivative(const DIAGONAL_MATRIX<T,d>& Sigma,const T p,const T mu,const T alpha) {PHYSBAM_FATAL_ERROR();}
//#################################################
};

}

#endif

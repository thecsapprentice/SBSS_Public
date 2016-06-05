//#####################################################################
// Copyright 2013, Nathan Mitchell, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ANALYTIC_COLLISION
//#####################################################################
#ifndef __ANALYTIC_COLLISION_H__
#define __ANALYTIC_COLLISION_H__

#include "COLLISION_INTERFACE.h"

namespace PhysBAM{

    enum ANALYTIC_COLLISION_SHAPE { AC_SPHERE, AC_PLANE };

    template<class T, int d, ANALYTIC_COLLISION_SHAPE shape> 
        class ANALYTIC_COLLISION : public COLLISION_INTERFACE<T,d>{
    public:
        ANALYTIC_COLLISION(){ };
    private:
        virtual T Phi_Implementation( const T world_location[d]) const  {return T();};
        virtual void Normal_Implementation(const T world_location[d], T normal[d]) const {
            for(int w=1;w<=d;w++) normal[w-1] = VECTOR<T,d>::Axis_Vector(1)(w);
        };
    };

    template<class T, int d> 
        class ANALYTIC_COLLISION<T, d ,AC_SPHERE> : public COLLISION_INTERFACE<T,d>{
    private:
        VECTOR<T,d> _center;
        T _radius;
    public:
        ANALYTIC_COLLISION(const VECTOR<T,d>& center, T radius) {_center=center;_radius=radius;};
    private:
        virtual T Phi_Implementation( const T world_location[d]) const {
            const VECTOR<T,d>& _world_location=*(const VECTOR<T,d>*)(world_location);
            return (_world_location-_center).Magnitude() - _radius;
        };
        virtual void Normal_Implementation(const T world_location[d], T normal[d]) const {
            const VECTOR<T,d>& _world_location=*(const VECTOR<T,d>*)(world_location);
            VECTOR<T,d> normalv = (_world_location-_center).Normalized();
            for(int w=1;w<=d;w++) normal[w-1] = normalv;
        };
    };

    template<class T, int d> 
        class ANALYTIC_COLLISION<T, d ,AC_PLANE> : public COLLISION_INTERFACE<T,d>{
    private:
        VECTOR<T,d> _normal;
        VECTOR<T,d> _center;
    public:
        ANALYTIC_COLLISION(const VECTOR<T,d>& normal, const VECTOR<T,d>& center) {_center=center;
            _normal=normal.Normalized();};
    private:
        virtual T Phi_Implementation( const T world_location[d]) const {
            const VECTOR<T,d>& _world_location=*(const VECTOR<T,d>*)(world_location);
            return VECTOR<T,d>::Dot_Product((_world_location-_center),_normal); 
        };
        virtual  void Normal_Implementation(const T world_location[d], T normal[d]) const {
            for(int w=1;w<=d;w++) normal[w-1] = _normal(w);
        };
    };



}

#endif

//#####################################################################
// Copyright 2013, Nathan Mitchell, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class COLLISION_INTERFACE
//#####################################################################
#ifndef __COLLISION_INTERFACE_H__
#define __COLLISION_INTERFACE_H__

#include <PhysBAM_Tools/Vectors/VECTOR_FORWARD.h>

namespace PhysBAM{

    template<class T, int d> 
        class COLLISION_INTERFACE{
    public:
        T Phi( const VECTOR<T,d>& world_location) const{
            T _world_location[d];
            for(int w=1;w<=d;w++) _world_location[w-1] = world_location(w);
            return Phi_Implementation( _world_location );
        };
        VECTOR<T,d> Normal( const VECTOR<T,d>& world_location) const{
            T _world_location[d], _normal[d];
            VECTOR<T,d> normal;
            for(int w=1;w<=d;w++) _world_location[w-1] = world_location(w);
            Normal_Implementation( _world_location, _normal);
            for(int w=1;w<=d;w++) normal(w) = _normal[w-1];
            return normal;
        };
        VECTOR<T,d> ClosestPoint(const VECTOR<T,d>& world_location) const{
            return world_location - (Phi(world_location)*Normal(world_location)); 
        }
    private:
        virtual T Phi_Implementation( const T world_location[d] ) const = 0;
        virtual void Normal_Implementation( const T world_location[d], T normal[d] ) const = 0;
    };

}

#endif

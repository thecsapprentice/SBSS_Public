//#####################################################################
// Copyright 2014, Nathan Mitchell
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class C_SHAPE
//##################################################################### 
#ifndef __C_SHAPE__
#define __C_SHAPE__

#include <PhysBAM_Tools/Log/DEBUG_UTILITIES.h>
#include <PhysBAM_Tools/Math_Tools/constants.h>
#include <PhysBAM_Tools/Math_Tools/pow.h>
#include <PhysBAM_Tools/Math_Tools/RANGE.h>
namespace PhysBAM{

template<class T>
class C_SHAPE
{
    typedef VECTOR<T,2> TV;
    enum WORKAROUND {d=TV::m};
public:
    typedef TV VECTOR_T;

    TV center;
    T r1, r2;
    T omega;

 C_SHAPE() : center(TV()), r1(.9), r2(1), omega(T(one_sixth_pi/T(2)))
    {}
    
    C_SHAPE( const TV& center, const T r1, const T r2, const T omega )
        : center(center), r1(r1), r2(r2), omega(omega)
    {}
    
    T Signed_Distance( const TV& location ) const {
        TV loc_norm = location.Normalized();
        T distance_from_center = (location-center).Magnitude();
        T distance_to_r1 = abs(distance_from_center-r1);
        T distance_to_r2 = abs(distance_from_center-r2);
        T distance_from_core = min( distance_to_r1, distance_to_r2 );
        T signed_distance_from_core = (distance_from_center > r1 && distance_from_center < r2) ? -distance_from_core : distance_from_core;
        
        T cap_radius = (r2-r1)/2.0;
        int loc_top = asin( loc_norm.y ) >= 0 ? 1 : -1;        
        TV closest_cap_center( cos(loc_top*omega)*(r1+cap_radius), sin(loc_top*omega)*(r1+cap_radius) );
        T distance_from_cap = (location-closest_cap_center).Magnitude();
        T signed_distance_from_cap = distance_from_cap - cap_radius;
        

        T  loc_angle = acos( loc_norm.x );
        if( loc_angle <= omega )
            return signed_distance_from_cap;
        else
            return signed_distance_from_core;
        
    }
    
    TV Normal( const TV& location ){
        return TV();

    }

    RANGE<TV> Bounding_Box() const
    {return RANGE<TV>(center).Thickened(r2);}


};


}


#endif

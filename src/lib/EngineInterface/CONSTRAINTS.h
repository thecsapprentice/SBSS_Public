//#####################################################################
// Copyright 2010-2013, Eftychios Sifakis, Nathan Mitchell
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class CONSTRAINT_NODE
// Class CONSTRAINT_SEGMENT
//#####################################################################
#ifndef __CONSTRAINTS_H__
#define __CONSTRAINTS_H__

#include <iostream>
#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Log/DEBUG_UTILITIES.h>


namespace PhysBAM{

template<class T, int d>
struct CONSTRAINT_NODE
{
    CONSTRAINT_NODE() : type(KINEMATIC) {
        for(int i=0;i<d;i++)
            _spatial_location[i] = (T)(0.0);
    }

    CONSTRAINT_NODE( const CONSTRAINT_NODE& other ){
        type=other.type;
        grid_index() = other.grid_index();
        multilinear_coordinates() = other.multilinear_coordinates();
    }

    CONSTRAINT_NODE& operator= (const CONSTRAINT_NODE& other){
        if(&other != this){
            type=other.type;
            grid_index() = other.grid_index();
            multilinear_coordinates() = other.multilinear_coordinates();
        }
        return *this;
    }

    typedef VECTOR<T,d> TV;
    typedef VECTOR<int,d> T_INDEX;
    enum {
        constraint_kinematic = 0x1,
        constraint_embedded  = 0x2,
        constraint_mesh      = 0x4, // Negation of this implies grid - N/A with non-embedded constraints
        constraint_pointer   = 0x8  // Negation of this implies absolute - N/A with embedded constraints
    };

    typedef enum {
        GRID_FIXED    = constraint_embedded,
        MESH_FIXED    = constraint_embedded|constraint_mesh,
        KINEMATIC     = constraint_kinematic,
        KINEMATIC_PTR = constraint_kinematic|constraint_pointer
    } CONSTRAINT_ENDPOINT_TYPE;
            
    CONSTRAINT_ENDPOINT_TYPE type;
            
    union{
        int _grid_index[d];
        int _mesh_index;
        T _spatial_location[d];
        int _spatial_location_ptr;
    };

    // union{
    //     float _spatial_velocity[d];
    //     int _spatial_velocity_ptr;
    // };

    union{
        T _multilinear_coordinates[d];
    };

    const T_INDEX& grid_index() const {return *reinterpret_cast<const T_INDEX*>(_grid_index);}
    T_INDEX& grid_index() {return *reinterpret_cast<T_INDEX*>(_grid_index);}
    const int& mesh_index() const {return _mesh_index;}
    int& mesh_index() {return _mesh_index;}
    const TV& spatial_location() const {return *reinterpret_cast<const TV*>(_spatial_location);}
    TV& spatial_location() {return *reinterpret_cast<TV*>(_spatial_location);}
    const TV& spatial_location_ptr(const ARRAY<TV> locations) const {return locations(_spatial_location_ptr);}
    TV& spatial_location_ptr(const ARRAY<TV> locations) {return locations(_spatial_location_ptr);}

    const TV& multilinear_coordinates() const {return *reinterpret_cast<const TV*>(_multilinear_coordinates);}
    TV& multilinear_coordinates() {return *reinterpret_cast<TV*>(_multilinear_coordinates);}

    // const TV& spatial_velocity() const {return *reinterpret_cast<const TV*>(_spatial_velocity);}
    // TV& spatial_velocity() {return *reinterpret_cast<TV*>(_spatial_velocity);}
    // const TV& spatial_velocity_ptr(const ARRAY<TV> velocities) const {return velocities(_spatial_velocity_ptr);}
    // TV& spatial_velocity_ptr(const ARRAY<TV> velocities) {return velocities(_spatial_velocity_ptr);}

};
    

template<class T, int d>
struct CONSTRAINT_SEGMENT
{
    CONSTRAINT_SEGMENT() : is_reference(false), spring_coefficient((T)(0.0)) {}

    CONSTRAINT_SEGMENT( const CONSTRAINT_SEGMENT& other ){
        endpoints[0] = other.endpoints[0];
        endpoints[1] = other.endpoints[1];
        is_reference = other.is_reference;
        spring_coefficient = other.spring_coefficient;
    }
    CONSTRAINT_SEGMENT& operator= (const CONSTRAINT_SEGMENT& other){
        if(&other != this){
            endpoints[0] = other.endpoints[0];
            endpoints[1] = other.endpoints[1];
            is_reference = other.is_reference;
            spring_coefficient = other.spring_coefficient;
        }
        return *this;
    }

    CONSTRAINT_NODE<T,d> endpoints[2];
    bool is_reference;
    union{
        T spring_coefficient;
        int spring_coefficient_ptr;
    };
    
    bool isSinglePoint() const {return (

                                        (endpoints[0].type & CONSTRAINT_NODE<T,d>::constraint_embedded && (!(endpoints[1].type & CONSTRAINT_NODE<T,d>::constraint_embedded))) || 
                                        ((!(endpoints[0].type & CONSTRAINT_NODE<T,d>::constraint_embedded)) && (endpoints[1].type & CONSTRAINT_NODE<T,d>::constraint_embedded))

                                        ) ;}
    bool isDualPoint() const {return (((endpoints[0].type & CONSTRAINT_NODE<T,d>::constraint_embedded) && 
                                       ((endpoints[1].type & CONSTRAINT_NODE<T,d>::constraint_embedded)))); }


    const CONSTRAINT_NODE<T,d>& getEmbeddedPoint() const {PHYSBAM_ASSERT( isSinglePoint() ); 
        if(endpoints[0].type & CONSTRAINT_NODE<T,d>::constraint_embedded)
            return endpoints[0];
        if(endpoints[1].type & CONSTRAINT_NODE<T,d>::constraint_embedded)
            return endpoints[1];
        PHYSBAM_FATAL_ERROR("Should never reach here!");
    }
};


template<class T, int d>
std::ostream& operator<<( std::ostream& out, const CONSTRAINT_NODE<T,d>& cn){
    out << "TYPE = ";
    switch(cn.type){
    case CONSTRAINT_NODE<T,d>::GRID_FIXED:
        out << "GRID" << std::endl;
        out << "INDEX = " << cn.grid_index() << std::endl;
        out << "M_COORDS = " << cn.multilinear_coordinates() << std::endl;
        break;
    case CONSTRAINT_NODE<T,d>::MESH_FIXED:
        out << "MESH" << std::endl;
        out << "INDEX = " << cn.mesh_index() << std::endl;
        out << "M_COORDS = " << cn.multilinear_coordinates() << std::endl;
        break;
    case CONSTRAINT_NODE<T,d>::KINEMATIC:
        out << "KINEMATIC" << std::endl;
        out << "LOCATION = " << cn.spatial_location() << std::endl;
        break;
    case CONSTRAINT_NODE<T,d>::KINEMATIC_PTR:
        out << "KINEMATIC_PTR" << std::endl;
        out << "POINTER = " << cn._spatial_location_ptr << std::endl;
        break;
    default:
        out << "UNKNOWN" << std::endl;
        break;
    }
    return out;
}


template<class T, int d>
std::ostream& operator<<( std::ostream& out, const CONSTRAINT_SEGMENT<T,d>& cs){
    out << "Endpoint 1" << std::endl;
    out << cs.endpoints[0] << std::endl;
    out << "Endpoint 2" << std::endl;
    out << cs.endpoints[1] << std::endl;
    if(cs.is_reference)
        out << "Spring Constant PTR: " << cs.spring_coefficient_ptr << std::endl;
    else
        out << "Spring Constant: " << cs.spring_coefficient << std::endl;
    return out;
}

}
#endif 

//#####################################################################
// Copyright 2014, Nathan Mitchell
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class HYBRID_ARRAY
//#####################################################################
#ifndef _HYBRID_ARRAY_H_
#define _HYBRID_ARRAY_H_

#include <PhysBAM_Tools/Data_Structures/PAIR.h>
#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Arrays_Computations/ARRAY_COPY.h>
#include <PhysBAM_Tools/Arrays_Computations/ARRAY_MIN_MAX.h>
#include <PhysBAM_Tools/Arrays_Computations/SUMMATIONS.h>
#include <iostream>

namespace PhysBAM{

    template <class T, int d>
    class HYBRID_ARRAY 
    {
    public:

        template<class T2> struct REBIND{typedef HYBRID_ARRAY<T2,d> TYPE;};
        typedef int ID;
        typedef VECTOR<int,d> TV_INT;
        enum WORKAROUND1 {dimension=d};
        typedef T ELEMENT;
        typedef PAIR<TV_INT,ID> INDEX;

        ARRAY<T,ID> flat_array;
        ARRAY<T,TV_INT> nd_array;

        HYBRID_ARRAY() :
            flat_array(0, true), nd_array(RANGE<TV_INT>(), true)
        {};
        HYBRID_ARRAY( const RANGE<TV_INT>& domain_input, const ID m_input,
                      const bool initialize_using_default_constructor=true ) : 
            flat_array(m_input, initialize_using_default_constructor),
            nd_array(domain_input, initialize_using_default_constructor)
        { }

        HYBRID_ARRAY( const HYBRID_ARRAY& array )
        {
            flat_array( array.flat_array );
            nd_array( array.nd_array );
        }

        HYBRID_ARRAY& operator=(const HYBRID_ARRAY& source ){
            if( &source == this )
                return *this;
            
            flat_array = source.flat_array;
            nd_array = source.nd_array;
            return *this;
        }

        void Resize(const RANGE<TV_INT>& domain_input, const ID m_input,
                    const bool initialize_new_elements=true, const bool copy_existing_elements=true,
                    const T& initialization_value=T()){
            flat_array.Resize( m_input, initialize_new_elements, copy_existing_elements, initialization_value);
            nd_array.Resize( domain_input, initialize_new_elements, copy_existing_elements, initialization_value);
        }

        void Resize(const ID m_input,
                    const bool initialize_new_elements=true, const bool copy_existing_elements=true,
                    const T& initialization_value=T()){
            flat_array.Resize( m_input, initialize_new_elements, copy_existing_elements, initialization_value);
        }

        void Resize(const RANGE<TV_INT>& domain_input,
                    const bool initialize_new_elements=true, const bool copy_existing_elements=true,
                    const T& initialization_value=T()){
            nd_array.Resize( domain_input, initialize_new_elements, copy_existing_elements, initialization_value);
        }

        RANGE<TV_INT> Domain_Indices() const
        {return nd_array.Domain_Indices();}

        ID Size() const
        {return nd_array.Size().Product() + flat_array.Size();}

        ID Flat_Size() const
        {return flat_array.Size();}

        TV_INT Domain_Size() const
        {return nd_array.Size();}

        void Fill(const T& constant)
        {
            ARRAYS_COMPUTATIONS::Fill(flat_array,constant);
            nd_array.Fill( constant );
        }

        bool Valid_Index( const INDEX& index ){            
            return (index.y == 0 && nd_array.Valid_Index( index.x )) ||
                (index.y != 0 && flat_array.Valid_Index(index.y));
        }

        T& operator()(const INDEX& index ) {
            if(index.y == 0)
                return nd_array(index.x);
            else
                return flat_array(index.y);
        } 

        const T& operator()(const INDEX& index ) const {
            if(index.y == 0)
                return nd_array(index.x);
            else
                return flat_array(index.y);
        }

        INDEX Append( const T& element){
            ID i = flat_array.Append(element);
            return INDEX(TV_INT(),i);            
        }

        T Sum() const
        {
            return ARRAYS_COMPUTATIONS::Sum(flat_array) + nd_array.Sum();
        }
      
    };

    template<int d>
    std::ostream& operator<<(std::ostream& out, const PAIR<VECTOR<int,d>,int>& index ){
        out << "{ " << index.x << " , " << index.y << " }";
        return out;
    }

}

#endif

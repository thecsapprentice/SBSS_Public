//#####################################################################
// Copyright 2014, Nathan Mitchell
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class HYBRID_ARRAY_ITERATOR
//#####################################################################
#ifndef __HYBRID_ARRAY_ITERATOR__
#define __HYBRID_ARRAY_ITERATOR__

namespace PhysBAM{

template<class T, int d> class HYBRID_ARRAY;

template<int d> 
class HYBRID_ARRAY_ITERATOR
{
    typedef VECTOR<int,d> TV_INT;
    typedef RANGE<TV_INT> T_RANGE;
    typedef PAIR<TV_INT,int> HINDEX;
    enum workaround {dimension=d};

    const T_RANGE& domain_range;
    const int flat_range; 
    TV_INT domain_index;
    int flat_index;
    HINDEX index;

public:
    template<class T>
    HYBRID_ARRAY_ITERATOR(const HYBRID_ARRAY<T,d>& array_input)
        :domain_range(array_input.Domain_Indices()), flat_range(array_input.Flat_Size())
    {
        Reset();
    }

    HYBRID_ARRAY_ITERATOR(const T_RANGE& domain_range_input, const int flat_range_input)
        :domain_range(domain_range_input), flat_range(flat_range_input)
    {
        Reset();
    }


    void Reset()
    {
        domain_index=domain_range.min_corner;
        flat_index=0;
        index = HINDEX(domain_index,flat_index);
    }

    bool Valid() const
    {
        return (domain_index.x<=domain_range.max_corner.x && flat_index == 0) ||
               (domain_index == TV_INT() && flat_index <= flat_range && flat_index > 0);
    }

    void Next()
    {
        if(flat_index == 0){
            for(int i=dimension;i>=1;i--)
                if(domain_index(i)<domain_range.max_corner(i) || i==1){
                    domain_index(i)++; 
                    if(domain_index.x > domain_range.max_corner.x){
                        domain_index = TV_INT();
                        flat_index++;
                    }
                    index = HINDEX(domain_index,flat_index);
                    return;
                }
                else 
                    domain_index(i)=domain_range.min_corner(i);
        }
        else{
            flat_index++;
            index = HINDEX(domain_index,flat_index);
            return;
        }
        index = HINDEX(domain_index,flat_index);
    }

    HINDEX& Index()
    {        
        return index;
    }
};

//#####################################################################
}
#endif

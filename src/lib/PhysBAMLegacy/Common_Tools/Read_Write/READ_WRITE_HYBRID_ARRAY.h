//#####################################################################
// Copyright 2014, Nathan Mitchell.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class READ_WRITE_HYBRID_ARRAY
//#####################################################################
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
#ifndef __READ_WRITE_HYBRID_ARRAY__
#define __READ_WRITE_HYBRID_ARRAY__

#include <PhysBAM_Tools/Read_Write/Utilities/READ_WRITE.h>
#include <Common_Tools/Arrays/HYBRID_ARRAY.h>
namespace PhysBAM{

    template<class RW,class T,int d>
    class Read_Write<HYBRID_ARRAY<T,d>,RW>
{
public:
    static void Read(std::istream& input,HYBRID_ARRAY<T,d>& object)
    {Read_Binary<RW>(input,object.nd_array,object.flat_array);}

    static void Write(std::ostream& output,const HYBRID_ARRAY<T,d>& object)
    {Write_Binary<RW>(output,object.nd_array,object.flat_array);}
};
}
#endif
#endif

//#####################################################################
// Copyright 2013, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class POWER
//#####################################################################
#ifndef __POWER__
#define __POWER__

#include <PhysBAM_Tools/Utilities/STATIC_ASSERT.h>
namespace PhysBAM{

template<unsigned base, unsigned exponent> struct POWER;
template<unsigned base> struct POWER<base,0> { enum WORKAROUND {value = 1}; };
template<unsigned base,unsigned exponent> struct POWER { enum WORKAROUND {value=base*POWER<base,exponent-1>::value}; };

}
#endif

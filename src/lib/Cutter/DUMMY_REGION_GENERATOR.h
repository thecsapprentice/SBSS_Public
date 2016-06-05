#ifndef __DUMMY_REGION_GENERATOR_H__
#define __DUMMY_REGION_GENERATOR_H__


#include "REGION_GENERATOR.h"


namespace PhysBAM{

    template<class T, int d>
        class DummyRegionGenerator : public REGION_GENERATOR<T,d>
        {
            private:
            typedef REGION_GENERATOR<T,d> BASE;

            typedef typename BASE::T_INDEX T_INDEX;
            typedef typename BASE::TV TV;
            typedef typename BASE::T_CELL T_CELL;
            typedef typename BASE::T_FACE T_FACE;
            typedef typename BASE::T_GRID T_GRID;

 
            virtual void ResolveCuts();
            virtual bool IsMaterialContinous(const T_INDEX& indexA, const T_CELL& cellA, const T_FACE& faceA,
                                             const T_INDEX& indexB, const T_CELL& cellB, const T_FACE& faceB);                        
            
        public:
            DummyRegionGenerator( const T_GRID& grid_input);
    
            
        };

}


#endif

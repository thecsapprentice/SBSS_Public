#include "DUMMY_REGION_GENERATOR.h"

using namespace PhysBAM;



template<class T, int d> 
DummyRegionGenerator<T,d>::DummyRegionGenerator(const T_GRID& grid_input) : BASE( grid_input)
{

} 


template<class T, int d> void
DummyRegionGenerator<T,d>::ResolveCuts()
{




}


template<class T, int d> bool
DummyRegionGenerator<T,d>::IsMaterialContinous(const T_INDEX& indexA, const T_CELL& cellA, const T_FACE& faceA,
                                               const T_INDEX& indexB, const T_CELL& cellB, const T_FACE& faceB)
{



    return true;
}


template class DummyRegionGenerator<float, 2>;
template class DummyRegionGenerator<float, 3>;
template class DummyRegionGenerator<double, 2>;
template class DummyRegionGenerator<double, 3>;

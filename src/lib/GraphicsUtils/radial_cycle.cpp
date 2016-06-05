#include "bmesh.h"
#include <algorithm>

using namespace BMESH;

template<class ANCHOR, class LEAF, int cycle_tag>
typename ELEMENT_CYCLE<ANCHOR,LEAF,cycle_tag>::T_LINK* 
ELEMENT_CYCLE<ANCHOR,LEAF,cycle_tag>::get_cycle_link( LEAF l, ANCHOR a ) {
    T_LINK* rl;
    if( l->edge == a ) {// Are we the edge of this loop
        rl = &(l->radial_cycle);
    }
    else{ // Oops. We aren't associated with this loop. Big Problem.
        throw BMeshError( "Edge-Loop mismatch during RADIAL_CYCLE::get_cycle_link" );
    }
    return rl;
}

template<class ANCHOR, class LEAF, int cycle_tag>
LEAF
ELEMENT_CYCLE<ANCHOR,LEAF,cycle_tag>::next( ANCHOR a, LEAF l, bool nothrow){
    LEAF next_loop = NULL;
    if( a == l->edge )
        next_loop = l->radial_cycle.next;
    else
        throw BMeshError( "Edge-Loop mismatch during RADIAL_CYCLE::next" );   
    if( next_loop == NULL && !nothrow)
        throw BMeshError( "Broken cycle detected during RADIAL_CYCLE::next" );   
    return next_loop;
}

template<class ANCHOR, class LEAF, int cycle_tag>
LEAF
ELEMENT_CYCLE<ANCHOR,LEAF,cycle_tag>::prev( ANCHOR a, LEAF l, bool nothrow){
    LEAF prev_loop = NULL;
    if( a == l->edge )
        prev_loop = l->radial_cycle.prev;
    else
        throw BMeshError( "Edge-Loop mismatch during RADIAL_CYCLE::prev" );   
    if( prev_loop == NULL && !nothrow)
        throw BMeshError( "Broken cycle detected during RADIAL_CYCLE::prev" );   
    return prev_loop;
}

template<class ANCHOR, class LEAF, int cycle_tag>
bool
ELEMENT_CYCLE<ANCHOR,LEAF,cycle_tag>::validate( ANCHOR a ){

    if( a->base == NULL ) // this may be a problem. but from the anchor's perpective, its fine
        return true;

    // Now that we have a base, traverse the cycle. We should not encounter any dead ends, and each
    // component should point to the anchor.
    LEAF current = a->base;
    do{
        if( !(current->edge == a) )
            throw BMeshError( "Broken link between Edge/Loop during RADIAL_CYCLE::validate" );

        current = next( a, current, false );
        if( current == NULL )
            throw BMeshError( "Cycle dead-ends during RADIAL_CYCLE::validate" );
    }while( current != a->base ); // Continue until we reach the base again
       
    return true;
}

// Radial Cycle
template struct ELEMENT_CYCLE< BMeshPolicy<float, 2>::T_EDGE_PTR,BMeshPolicy<float, 2>::T_LOOP_PTR, RADIAL_CYCLE_T>;
template struct ELEMENT_CYCLE< BMeshPolicy<float, 3>::T_EDGE_PTR,BMeshPolicy<float, 3>::T_LOOP_PTR, RADIAL_CYCLE_T>;
template struct ELEMENT_CYCLE< BMeshPolicy<double, 2>::T_EDGE_PTR,BMeshPolicy<double, 2>::T_LOOP_PTR, RADIAL_CYCLE_T>;
template struct ELEMENT_CYCLE< BMeshPolicy<double, 3>::T_EDGE_PTR,BMeshPolicy<double, 3>::T_LOOP_PTR, RADIAL_CYCLE_T>;

#include "bmesh.h"
#include <algorithm>

using namespace BMESH;

template<class ANCHOR, class LEAF, int cycle_tag>
typename ELEMENT_CYCLE<ANCHOR,LEAF,cycle_tag>::T_LINK* 
ELEMENT_CYCLE<ANCHOR,LEAF,cycle_tag>::get_cycle_link( LEAF l, ANCHOR a ) {
    T_LINK* dl;
    if( l->vertex1 == a ) {// Are we the first vertex?
        dl = &(l->v1_disk_cycle);
    }
    else if( l->vertex2 == a ) {// Are we the second vertex?
        dl = &(l->v2_disk_cycle);
        }
    else{ // Oops. We aren't associated with this vertex. Big Problem.
        throw BMeshError( "Vertex-Edge mismatch during DISK_CYCLE::get_cycle_link" );
    }
    return dl;
}

template<class ANCHOR, class LEAF, int cycle_tag>
LEAF
ELEMENT_CYCLE<ANCHOR,LEAF,cycle_tag>::next( ANCHOR a, LEAF l, bool nothrow){
    LEAF next_edge;
    if( !(a == l->vertex1 || a == l->vertex2) )
        throw BMeshError( "Vertex-Edge mismatch during DISK_CYCLE::next" );
    if( a == l->vertex1 )
        next_edge = l->v1_disk_cycle.next;
    if( a == l->vertex2 )
        next_edge = l->v2_disk_cycle.next;
    if( next_edge == NULL && !nothrow)
        throw BMeshError( "Broken cycle detected during DISK_CYCLE::next" );   
    return next_edge;
}

template<class ANCHOR, class LEAF, int cycle_tag>
LEAF
ELEMENT_CYCLE<ANCHOR,LEAF,cycle_tag>::prev( ANCHOR a, LEAF l, bool nothrow){
    LEAF prev_edge;
    if( !(a == l->vertex1 || a == l->vertex2) )
        throw BMeshError( "Vertex-Edge mismatch during DISK_CYCLE::next" );
    if( a == l->vertex1 )
        prev_edge = l->v1_disk_cycle.prev;
    if( a == l->vertex2 )
        prev_edge = l->v2_disk_cycle.prev;
    if( prev_edge == NULL && !nothrow)
        throw BMeshError( "Broken cycle detected during DISK_CYCLE::prev" );   
    return prev_edge;
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
        if( !(current->vertex1 == a || current->vertex2 == a) )
            throw BMeshError( "Broken link between Vertex/Edge during DISK_CYCLE::validate" );
        if( current->vertex1 == current->vertex2 )
            throw BMeshError( "Edge points to vertex twice during DISK_CYCLE::validate" );

        current = next( a, current, false );
        if( current == NULL )
            throw BMeshError( "Cycle dead-ends during DISK_CYCLE::validate" );
    }while( current != a->base ); // Continue until we reach the base again
       
    return true;
}


// Disk Cycle
template struct ELEMENT_CYCLE< BMeshPolicy<float, 2>::T_VERT_PTR,BMeshPolicy<float, 2>::T_EDGE_PTR, DISK_CYCLE_T>;
template struct ELEMENT_CYCLE< BMeshPolicy<float, 3>::T_VERT_PTR,BMeshPolicy<float, 3>::T_EDGE_PTR, DISK_CYCLE_T>;
template struct ELEMENT_CYCLE< BMeshPolicy<double, 2>::T_VERT_PTR,BMeshPolicy<double, 2>::T_EDGE_PTR, DISK_CYCLE_T>;
template struct ELEMENT_CYCLE< BMeshPolicy<double, 3>::T_VERT_PTR,BMeshPolicy<double, 3>::T_EDGE_PTR, DISK_CYCLE_T>;

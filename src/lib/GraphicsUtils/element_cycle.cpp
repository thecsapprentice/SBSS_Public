#include "bmesh.h"
#include <algorithm>

using namespace BMESH;

template<class ANCHOR, class LEAF, int cycle_tag>
int
ELEMENT_CYCLE<ANCHOR,LEAF,cycle_tag>::size( ANCHOR a ){
    LEAF start = a->base;
    if( start == NULL )
        return 0;
    
    int size = 1;
    LEAF current = next( a, start, false );
    while( current != start ){
        size++;
        current = next( a, current, false );
    }
    return size;
}


template<class ANCHOR, class LEAF, int cycle_tag>
void
ELEMENT_CYCLE<ANCHOR,LEAF,cycle_tag>::append( ANCHOR a, LEAF l ){

    if( a->base != NULL ){ // Are we connected to any edge yet?
        T_LINK *el_new, *el_start, *el_end;

        // We need to insert this loop into the cycle based at f
        el_new = get_cycle_link( l, a );
        // This is our base handle into the cycle
        el_start = get_cycle_link( a->base, a );
        // This is the next cycle link from the base, if it exists...
        el_end = el_start->prev ? get_cycle_link( el_start->prev, a ) : NULL;

        // Stitch in new cycle link
        //        dl_start.p <-> n.dl_end
        //                    |
        //                    v
        // dl_start.p <-> n.dl_new.p <-> n.dl_end
        el_new->next = a->base;
        el_new->prev = el_start->prev;
        el_start->prev = l;
        if( el_end != NULL )
            el_end->next = l;
    }
    else{ // Nope!
        T_LINK* el = get_cycle_link( l, a );
        a->base = l; // Link the vertex with the edge;
        el->next = el->prev = l; // Build a loop cycle of one edge around the vertex;
    }
}

template<class ANCHOR, class LEAF, int cycle_tag>
void
ELEMENT_CYCLE<ANCHOR,LEAF,cycle_tag>::remove( ANCHOR a, LEAF l ){

    T_LINK* el_this_edge, *el_neighbor_edge;

    // Get our cycle link to this vertex
    el_this_edge = get_cycle_link( l, a ); 
    if( el_this_edge->prev ) { // If we are pointing to another edge in the prev direction, relink it to our next
        el_neighbor_edge = get_cycle_link( el_this_edge->prev, a );
        el_neighbor_edge->next = el_this_edge->next;
    }
    if( el_this_edge->next ) { // If we are pointing to another edge in the next direction, relink it to our prev
        el_neighbor_edge = get_cycle_link( el_this_edge->next, a );
        el_neighbor_edge->prev = el_this_edge->prev;
    }

    // if we are the edge for our cycle's base, switch it to our
    //        next edge if we aren't the only edge in the cycle.
    if( a->base == l )
        a->base = ( l != el_this_edge->next )? el_this_edge->next : NULL;
    
    // Now we are disconnected, so erase our links to eliminate dangling references...
    el_this_edge->next = el_this_edge->prev = NULL;
}


template<class ANCHOR, class LEAF, int cycle_tag>
bool
ELEMENT_CYCLE<ANCHOR,LEAF,cycle_tag>::contains( ANCHOR a, LEAF l){
    LEAF start = a->base;
    if( start == NULL ) // If the anchor doesn't have a base, then there is no cycle to be part of
        return false;
    
    LEAF current = next( a, start, false ); // Traverse the whole cycle
    while( current != start ){
        if( current = l ) // If we discover ourselves, stop
            return true; 
        current = next( a, current, false );        
    }
    return false; // Otherwise, we are not in this cycle after all
}


// Disk Cycle
template struct ELEMENT_CYCLE< BMeshPolicy<float, 2>::T_VERT_PTR,BMeshPolicy<float, 2>::T_EDGE_PTR, DISK_CYCLE_T>;
template struct ELEMENT_CYCLE< BMeshPolicy<float, 3>::T_VERT_PTR,BMeshPolicy<float, 3>::T_EDGE_PTR, DISK_CYCLE_T>;
template struct ELEMENT_CYCLE< BMeshPolicy<double, 2>::T_VERT_PTR,BMeshPolicy<double, 2>::T_EDGE_PTR, DISK_CYCLE_T>;
template struct ELEMENT_CYCLE< BMeshPolicy<double, 3>::T_VERT_PTR,BMeshPolicy<double, 3>::T_EDGE_PTR, DISK_CYCLE_T>;

// Loop Cycle
template struct ELEMENT_CYCLE< BMeshPolicy<float, 2>::T_FACE_PTR,BMeshPolicy<float, 2>::T_LOOP_PTR, LOOP_CYCLE_T>;
template struct ELEMENT_CYCLE< BMeshPolicy<float, 3>::T_FACE_PTR,BMeshPolicy<float, 3>::T_LOOP_PTR, LOOP_CYCLE_T>;
template struct ELEMENT_CYCLE< BMeshPolicy<double, 2>::T_FACE_PTR,BMeshPolicy<double, 2>::T_LOOP_PTR, LOOP_CYCLE_T>;
template struct ELEMENT_CYCLE< BMeshPolicy<double, 3>::T_FACE_PTR,BMeshPolicy<double, 3>::T_LOOP_PTR, LOOP_CYCLE_T>;

// Radial Cycle
template struct ELEMENT_CYCLE< BMeshPolicy<float, 2>::T_EDGE_PTR,BMeshPolicy<float, 2>::T_LOOP_PTR, RADIAL_CYCLE_T>;
template struct ELEMENT_CYCLE< BMeshPolicy<float, 3>::T_EDGE_PTR,BMeshPolicy<float, 3>::T_LOOP_PTR, RADIAL_CYCLE_T>;
template struct ELEMENT_CYCLE< BMeshPolicy<double, 2>::T_EDGE_PTR,BMeshPolicy<double, 2>::T_LOOP_PTR, RADIAL_CYCLE_T>;
template struct ELEMENT_CYCLE< BMeshPolicy<double, 3>::T_EDGE_PTR,BMeshPolicy<double, 3>::T_LOOP_PTR, RADIAL_CYCLE_T>;

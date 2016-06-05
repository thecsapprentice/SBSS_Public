#include "bmesh.h"
#include <algorithm>

using namespace BMESH;


template<class T, int d>
bool
BMESH_VERT<T,d>::validate(){
    if( base == NULL )
        return true;
    else
        return DISK_CYCLE::validate( this );    
}


template<class T, int d>
bool
BMESH_EDGE<T,d>::validate(){
    if( vertex1 == NULL || vertex2 == NULL )
        throw BMeshError( "Dangling edge detected during EDGE validation" );
    if( !DISK_CYCLE::contains( vertex1, this ) )
        throw BMeshError( "Edge is not contained by its first vertex during EDGE validation" );
    if( !DISK_CYCLE::contains( vertex2, this ) )
        throw BMeshError( "Edge is not contained by its second vertex during EDGE validation" );
    if( base == NULL ) 
        return true;
    else
        return RADIAL_CYCLE::validate( this );          
}


template<class T, int d>
bool
BMESH_LOOP<T,d>::validate(){
    if( vertex == NULL || edge == NULL || face == NULL )
        throw BMeshError( "Disconnected loop detected during LOOP validation." );
    if( !LOOP_CYCLE::contains( face, this ) )
        throw BMeshError( "Loop is not contained by its face during LOOP validation." );
    if( !RADIAL_CYCLE::contains( edge, this ) )
        throw BMeshError( "Loop is not contained by its edge during LOOP validation." );
    if( !(edge->vertex1 == vertex || edge->vertex2 == vertex) )
        throw BMeshError( "Loop does not match edge/vertex relationship during LOOP validation.");
    return true;
}


template<class T, int d>
bool
BMESH_FACE<T,d>::validate(){

    if( LOOP_CYCLE::size( this ) != len )
        throw BMeshError( "Internal face length does not match size of associated loop cycle during FACE validation.");
    if( base == NULL )
        return true;
    else
        return LOOP_CYCLE::validate( this );
}

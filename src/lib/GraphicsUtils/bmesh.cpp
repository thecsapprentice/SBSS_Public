#include "bmesh.h"
#include <algorithm>

using namespace BMESH;

template<class T, int d>
BMesh<T,d>::BMesh(){

}

template<class T, int d>
BMesh<T,d>::~BMesh(){
    for( T_C_FACE_ITR iter = _faces.begin(); iter != _faces.end(); iter++)
        delete *iter;
    _faces.clear();        

    for( T_C_LOOP_ITR iter = _loops.begin(); iter != _loops.end(); iter++)
        delete *iter;
    _loops.clear();

    for( T_C_EDGE_ITR iter = _edges.begin(); iter != _edges.end(); iter++)
        delete *iter;
    _edges.clear();

    for( T_C_VERT_ITR iter = _verts.begin(); iter != _verts.end(); iter++)
        delete *iter;
    _verts.clear();
}


template<class T, int d>
typename BMesh<T,d>::IDX
BMesh<T,d>::vertex_make(){
    T_VERT_PTR vert = new T_VERT;
    vert->base = NULL;
    vert->position.assign(d, T(0));
    _verts.push_back( vert );

    validate(); // Validate whole BMesh structure before returning

    return IDX(_verts.size()-1);
}

 
 template<class T, int d>
 typename BMesh<T,d>::IDX
 BMesh<T,d>::edge_make(IDX vertex1, IDX vertex2){

     if( vertex1 >= _verts.size() )
         throw BMeshError( "Vertex 1 out of bounds during edge_make" );

     if( vertex2 >= _verts.size() )
         throw BMeshError( "Vertex 2 out of bounds during edge_make" );

     T_EDGE_PTR edge = new T_EDGE;
     edge->vertex1 = _verts.at(vertex1);
     edge->vertex2 = _verts.at(vertex2);
     edge->base = NULL;

     // Link the new edge into the cycles for its vertex
     DISK_CYCLE::append( edge->vertex1, edge );
     DISK_CYCLE::append( edge->vertex2, edge );

     _edges.push_back( edge );

    validate(); // Validate whole BMesh structure before returning

     return IDX(_edges.size()-1);
 }

template<class T, int d>
typename BMesh<T,d>::IDX
BMesh<T,d>::loop_make_internal( IDX vertex, IDX edge, IDX face ){

    if( vertex >= _verts.size() )
         throw BMeshError( "Vertex out of bounds during loop_make" );
    if( edge >= _edges.size() )
         throw BMeshError( "Edge out of bounds during loop_make" );
    if( face >= _faces.size() )
         throw BMeshError( "Face out of bounds during loop_make" );

    T_LOOP_PTR loop = new T_LOOP;
    loop->vertex = _verts.at(vertex);
    loop->edge = _edges.at(edge);
    loop->face = _faces.at(face);
    
    _loops.push_back( loop );

    return IDX(_loops.size()-1);
}


template<class T, int d>
 typename BMesh<T,d>::IDX
 BMesh<T,d>::face_make(const T_INDEX_COLLECTION& verts, const T_INDEX_COLLECTION& edges){
    
    if( verts.size() != edges.size() )
        throw BMeshError( "Number of vertices does not match number of edges in face_make" );

    C_IDX_ITR viter, eiter;
    for( viter = verts.begin(), eiter = edges.begin();
         viter != verts.end();
         viter++, eiter++){
        if( *viter >= _verts.size() )
            throw BMeshError( "Vertex out of bounds during face_make" );
        if( *eiter >= _edges.size() )
            throw BMeshError( "Edge out of bounds during face_make" );
        if( _edges.at( *eiter )->vertex1 != _verts.at( *viter ) &&
            _edges.at( *eiter )->vertex2 != _verts.at( *viter ) )
            throw BMeshError( "Edge does not contain vertex during face_make" );
    }

    T_FACE_PTR face = new T_FACE;
    face->len = 0;
    face->base = NULL;
    _faces.push_back( face );
    IDX face_id = IDX( _faces.size()-1 );


    // Set up loop cycle by creating first loop then add more
    viter = verts.begin(); eiter = edges.begin();
    IDX loop_id = _loop_make( *viter, *eiter, face_id );
    T_VERT_PTR vert = _verts.at(*viter);
    T_EDGE_PTR edge = _edges.at(*eiter);
    T_LOOP_PTR loop = _loops.at(loop_id);
    face->base = loop; // Important! Assign the base before appending loop into the cycle.
    LOOP_CYCLE::append( face, loop );    
    RADIAL_CYCLE::append( edge, loop );
    viter++; eiter++;

    // Now create and add the rest of the loops
    for( ; viter != verts.end();
           viter++, eiter++){      
        loop_id = _loop_make( *viter, *eiter, face_id );
        vert = _verts.at(*viter);
        edge = _edges.at(*eiter);
        loop = _loops.at(loop_id);
        LOOP_CYCLE::append( face, loop );
        RADIAL_CYCLE::append( edge, loop );
    }
    face->len = LOOP_CYCLE::size(face);

    validate(); // Validate whole BMesh structure before returning

    return face_id;
}




template<class T, int d>
void
BMesh<T,d>::face_edge_kill(IDX face){
    
    

}
 




template<class T, int d>
bool
BMesh<T,d>::validate(){
    for( T_C_VERT_ITR iter = _verts.begin(); iter != _verts.end(); iter++)
        (*iter)->validate();    

    for( T_C_EDGE_ITR iter = _edges.begin(); iter != _edges.end(); iter++)
        (*iter)->validate();

    for( T_C_LOOP_ITR iter = _loops.begin(); iter != _loops.end(); iter++)
        (*iter)->validate();

    for( T_C_FACE_ITR iter = _faces.begin(); iter != _faces.end(); iter++)
        (*iter)->validate();
}



template class BMesh<float,2>;
template class BMesh<float,3>;
template class BMesh<double,2>;
template class BMesh<double,3>;

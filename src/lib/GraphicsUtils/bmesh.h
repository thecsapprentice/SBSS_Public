#ifndef __BMESH_H__
#define __BMESH_H__

#include <vector>
#include <list>
#include <exception>
#include <string>

namespace BMESH {

// Declarations
template<class T, int d>
struct BMESH_VERT;
template<class T, int d>
struct BMESH_LOOP;
template<class T, int d>
struct BMESH_EDGE;
template<class T, int d>
struct BMESH_FACE;
template<class T, int d>
struct BMESH;

// Exception Definitions

class BMeshError : public std::exception {
public:
    BMeshError() throw() : _msg("Unknown BMesh Error.") {}
    BMeshError( std::string msg ) throw() : _msg(msg) {}
    virtual ~BMeshError() throw() {};
    virtual const char* what() const throw() {
        return _msg.c_str();
    }
private:
    std::string _msg;
};


// Type Definitions

#define DISK_CYCLE_T 0
#define LOOP_CYCLE_T 1
#define RADIAL_CYCLE_T 2

template<class E>
struct CYCLELINK{
    CYCLELINK(): next(NULL), prev(NULL) {}
    E next;
    E prev;
};

template<class T, int d>
struct BMeshPolicy{
    typedef BMESH_FACE<T,d> T_FACE;
    typedef BMESH_LOOP<T,d> T_LOOP;
    typedef BMESH_EDGE<T,d> T_EDGE;
    typedef BMESH_VERT<T,d> T_VERT;
    typedef T_VERT* T_VERT_PTR;
    typedef T_EDGE* T_EDGE_PTR;
    typedef T_LOOP* T_LOOP_PTR;
    typedef T_FACE* T_FACE_PTR;
    typedef CYCLELINK< T_EDGE_PTR > T_DISK_LINK;
    typedef CYCLELINK< T_LOOP_PTR > T_LOOP_LINK;
    typedef CYCLELINK< T_LOOP_PTR > T_RADIAL_LINK;
};


// Functions to manage geometric cycles
// There are three types of cycles in a BMesh:
//
// Disk Cycle
// These are cycles of edges around a vertex
// Anchor: Vertex
// Cycle Member: Edge
//
//        v --- edge
//        | \
//        |  \
//        |   edge
//        edge
//
//
// Loop Cycle
// These are cycles of loops around a face
// Anchor: Face
// Cycle Member: Loop
//
//         loop->
//      +---------+
//    ^ |         | |
//    | |  face   | v loop
// loop |         | 
//      +---------+
//        <- loop
//
//
//
// Radial Cycle
// These are cycles of faces(loops) around an edge
// Anchor: Edge
// Cycle Member: Loop
//
//         
//      +------------+------------+
//      |          ^ | l          | 
//      |  face    | e |   face   | 
//      |          l | v          | 
//      +------------+------------+
//      
//
template<class ANCHOR, class LEAF, int tag>
struct ELEMENT_CYCLE {
    typedef CYCLELINK< LEAF > T_LINK;
    static T_LINK* get_cycle_link( LEAF l, ANCHOR a );
    static int size( ANCHOR a );
    static void append( ANCHOR a, LEAF l );
    static void remove( ANCHOR a, LEAF l );
    static bool contains( ANCHOR a, LEAF l );
    static LEAF next( ANCHOR a, LEAF l, bool nothrow);
    static LEAF prev( ANCHOR a, LEAF l, bool nothrow);
    static bool validate( ANCHOR a);
};


template<class T, int d>
struct BMESH_VERT {
    typedef BMeshPolicy<T,d> BP;
    typedef typename BP::T_EDGE_PTR T_EDGE_PTR;

    std::vector<T> position;
    T_EDGE_PTR base; // Base into this vertex's disk cycle

    bool validate();
};

template<class T, int d>
struct BMESH_EDGE {
    typedef BMeshPolicy<T,d> BP;
    typedef typename BP::T_VERT_PTR T_VERT_PTR;
    typedef typename BP::T_LOOP_PTR T_LOOP_PTR;
    typedef typename BP::T_DISK_LINK T_DISK_LINK;

    T_VERT_PTR vertex1;
    T_DISK_LINK v1_disk_cycle;

    T_VERT_PTR vertex2;
    T_DISK_LINK v2_disk_cycle;

    T_LOOP_PTR base; // Base into this edge's radial cycle

    bool validate();
};

template<class T, int d>
struct BMESH_LOOP  {
    typedef BMeshPolicy<T,d> BP;
    typedef typename BP::T_VERT_PTR T_VERT_PTR;
    typedef typename BP::T_EDGE_PTR T_EDGE_PTR;
    typedef typename BP::T_FACE_PTR T_FACE_PTR;
    typedef typename BP::T_LOOP_LINK T_LOOP_LINK;
    typedef typename BP::T_RADIAL_LINK T_RADIAL_LINK;


    T_VERT_PTR vertex;
    T_EDGE_PTR edge;
    T_FACE_PTR face;

    T_LOOP_LINK loop_cycle;
    T_RADIAL_LINK radial_cycle;

    //Pointer to vertex specific data
    void *vertex_data;

    bool validate();
};


template<class T, int d>
struct BMESH_FACE {
    typedef BMeshPolicy<T,d> BP; 
    typedef typename BP::T_LOOP_PTR T_LOOP_PTR;

    // Number of loops in this face
    int len;
    
    // First loop of this face
    T_LOOP_PTR base; // Base into this face's loop cycle

    //Pointer to face specific data
    void *face_data;

    bool validate();
};

template<class T, int d>
class BMesh {  
public:
    typedef BMESH_FACE<T,d> T_FACE;
    typedef BMESH_LOOP<T,d> T_LOOP;
    typedef BMESH_EDGE<T,d> T_EDGE;
    typedef BMESH_VERT<T,d> T_VERT;
    typedef T_VERT* T_VERT_PTR;
    typedef T_EDGE* T_EDGE_PTR;
    typedef T_LOOP* T_LOOP_PTR;
    typedef T_FACE* T_FACE_PTR;

    typedef unsigned IDX;
    typedef std::vector<IDX> T_INDEX_COLLECTION;
    typedef std::vector<T_VERT_PTR> T_VERT_COLLECTION;
    typedef std::vector<T_EDGE_PTR> T_EDGE_COLLECTION;
    typedef std::vector<T_LOOP_PTR> T_LOOP_COLLECTION;
    typedef std::vector<T_FACE_PTR> T_FACE_COLLECTION;

    typedef T_INDEX_COLLECTION::iterator IDX_ITR;
    typedef T_INDEX_COLLECTION::const_iterator C_IDX_ITR;
    typedef typename T_VERT_COLLECTION::iterator T_VERT_ITR;
    typedef typename T_EDGE_COLLECTION::iterator T_EDGE_ITR;
    typedef typename T_LOOP_COLLECTION::iterator T_LOOP_ITR;
    typedef typename T_FACE_COLLECTION::iterator T_FACE_ITR;
    typedef typename T_VERT_COLLECTION::const_iterator T_C_VERT_ITR;
    typedef typename T_EDGE_COLLECTION::const_iterator T_C_EDGE_ITR;
    typedef typename T_LOOP_COLLECTION::const_iterator T_C_LOOP_ITR;
    typedef typename T_FACE_COLLECTION::const_iterator T_C_FACE_ITR;

    typedef ELEMENT_CYCLE<T_VERT_PTR, T_EDGE_PTR, DISK_CYCLE_T> DISK_CYCLE;
    typedef ELEMENT_CYCLE<T_FACE_PTR, T_LOOP_PTR, LOOP_CYCLE_T> LOOP_CYCLE;
    typedef ELEMENT_CYCLE<T_EDGE_PTR, T_LOOP_PTR, RADIAL_CYCLE_T> RADIAL_CYCLE;



private:

    T_FACE_COLLECTION _faces;
    T_LOOP_COLLECTION _loops;
    T_EDGE_COLLECTION _edges;
    T_VERT_COLLECTION _verts;

    IDX _loop_make( IDX vertex, IDX edge, IDX face );

public:

    BMesh();
    ~BMesh();

    // Make Methods
    IDX vertex_make();
    IDX edge_make(IDX vertex1, IDX vertex2);
    IDX face_make(const T_INDEX_COLLECTION& verts,
                  const T_INDEX_COLLECTION& edges);

    // Delete Methods

    // Remove all edges used by face and any additional faces that use those edges
    void face_edge_kill(IDX face); 

    // Remove all vertices used by face and any additional faces that use those vertices
    void face_vertex_kill(IDX face); 

    // Remove all loops used by face and any faces which use those loops
    void face_loop_kill(IDX face);

    // Remove this edge and any higher level elements that depend on it (loops, faces)
    void edge_kill( IDX edge );

    // Remove this vertex and any higher level elements that depend on it (edges, loops, faces)
    void vertex_kill( IDX vertex );

    // Highlevel Mesh Manipulations
    void split_edge_make_vertex();
    void join_edge_kill_vertex();
    void split_face_make_edge();
    void join_face_kill_edge();
    void reverse_face_loop();


    // Validation
    bool validate();
};

}

#endif

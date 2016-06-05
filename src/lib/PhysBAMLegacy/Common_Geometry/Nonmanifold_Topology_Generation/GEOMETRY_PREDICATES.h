#ifndef __GEOMETRY_PREDICATES_H__
#define __GEOMETRY_PREDICATES_H__

#include <Common_Geometry/Nonmanifold_Topology_Generation/MATERIAL_PREDICATE_TESSELATED_VOLUME.h>
#include <PhysBAM_Geometry/Basic_Geometry_Intersections/BOX_SEGMENT_2D_INTERSECTION.h>
#include <PhysBAM_Geometry/Basic_Geometry_Intersections/BOX_TRIANGLE_3D_INTERSECTION.h>
#include <PhysBAM_Geometry/Collisions_And_Interactions/ROBUST_SIMPLEX_INTERACTIONS.h>

#include <PhysBAM_Geometry/Topology_Based_Geometry/TETRAHEDRALIZED_VOLUME.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_AREA.h>
#include <PhysBAM_Geometry/Basic_Geometry/TRIANGLE_2D.h>
#include <PhysBAM_Geometry/Basic_Geometry/TETRAHEDRON.h>
#include <PhysBAM_Geometry/Basic_Geometry/TRIANGLE_3D.h>
#include <PhysBAM_Geometry/Basic_Geometry/SEGMENT_2D.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/SEGMENTED_CURVE_2D.h>

#include <exception>

namespace PhysBAM{
//==================================================================================================
//                                            Inside Tests
//
//==================================================================================================
    
    template<class T>
    struct GEOMETRY_PREDICATES {
        
        class GEOMETRY_PREDICATE_EXCEPTION : public std::exception {
        public:
            GEOMETRY_PREDICATE_EXCEPTION( bool i, bool r, std::string m ) :
                msg(m), intersects(i), robust(r) {}
            
            virtual ~GEOMETRY_PREDICATE_EXCEPTION() throw() {}
            
            virtual const char* what() const throw(){
                return msg.c_str();
            }
            
            std::string msg;
            bool intersects;
            bool robust;                                
        };
        

        static bool
        Inside( const TETRAHEDRON<T>& T1, const VECTOR<T,3>& P ) throw(GEOMETRY_PREDICATE_EXCEPTION);
        
        static bool
        Inside( const TRIANGLE_2D<T>& T1, const VECTOR<T,2>& P ) throw(GEOMETRY_PREDICATE_EXCEPTION);
        
        template<int d> static int
        InsideOutside( const RANGE<VECTOR<T,d> >& box, const typename MATERIAL_PREDICATE_TESSELATED_VOLUME_POLICY<T,d>::T_ELEMENT& T1 ) throw(GEOMETRY_PREDICATE_EXCEPTION);
        
//==================================================================================================
//                                            Intersection Tests
//
//        Referenced from: http://www.geometrictools.com/Documentation/MethodOfSeparatingAxes.pdf
//==================================================================================================
            
        template<int d> static int 
        WhichSide( const VECTOR< VECTOR<T,d>, d+1>& points, VECTOR<T,d> D, VECTOR<T,d> P ) throw(GEOMETRY_PREDICATE_EXCEPTION);
        
        static bool 
        TestIntersection( const TRIANGLE_2D<T>& T1, const SEGMENT_2D<T>& T2 ) throw(GEOMETRY_PREDICATE_EXCEPTION);
        
        static bool 
        TestIntersection( const TETRAHEDRON<T>& T1, const TRIANGLE_3D<T>& T2 ) throw(GEOMETRY_PREDICATE_EXCEPTION);
        
        static bool 
        TestIntersection( const TRIANGLE_2D<T>& T1, const TRIANGLE_2D<T>& T2 ) throw(GEOMETRY_PREDICATE_EXCEPTION);
            
        static bool 
        TestIntersection( const TETRAHEDRON<T>& T1, const TETRAHEDRON<T>& T2 ) throw(GEOMETRY_PREDICATE_EXCEPTION);
            
        static bool 
        TestIntersection( const RANGE<VECTOR<T,2> >& box, 
                          const typename MATERIAL_PREDICATE_TESSELATED_VOLUME_POLICY<T,2>::T_ELEMENT& T1 ) throw(GEOMETRY_PREDICATE_EXCEPTION);
            
        static bool 
        TestIntersection( const RANGE<VECTOR<T,3> >& box, 
                          const typename MATERIAL_PREDICATE_TESSELATED_VOLUME_POLICY<T,3>::T_ELEMENT& T1 ) throw(GEOMETRY_PREDICATE_EXCEPTION);

        static bool 
        TestIntersection( const VECTOR< VECTOR<T,3>, 8 >& hex, 
                          const typename MATERIAL_PREDICATE_TESSELATED_VOLUME_POLICY<T,3>::T_ELEMENT& T1 ) throw(GEOMETRY_PREDICATE_EXCEPTION);
            
        static bool 
        TestIntersection( const RANGE<VECTOR<T,2> >& box, 
                          const typename MATERIAL_PREDICATE_TESSELATED_VOLUME_POLICY<T,2>::T_BOUNDARY_ELEMENT& T1 ) throw(GEOMETRY_PREDICATE_EXCEPTION);
            
        static bool 
        TestIntersection( const RANGE<VECTOR<T,3> >& box, 
                          const typename MATERIAL_PREDICATE_TESSELATED_VOLUME_POLICY<T,3>::T_BOUNDARY_ELEMENT& T1 ) throw(GEOMETRY_PREDICATE_EXCEPTION);
            
//=====================================================================================
//
//                              Predicate Method Helpers
//
//=====================================================================================
            
        static bool
        Element_Neighbors( const TETRAHEDRON_MESH& mesh, int T1, int T2, VECTOR<int,3>& edge_nodes );
            
        static bool
        Element_Neighbors( const TRIANGLE_MESH& mesh, int T1, int T2, VECTOR<int,2>& edge_nodes );
            
        static bool
        Element_Neighbors( const SEGMENT_MESH& mesh, int T1, int T2, VECTOR<int,1>& edge_nodes );
            
        static HASHTABLE<int, ARRAY< VECTOR<int,2> > >
        GenerateBoundaryPrimatives( TRIANGLE_MESH& mesh, ARRAY<int> elements);
            
        static HASHTABLE<int, ARRAY< VECTOR<int,3> > > 
        GenerateBoundaryPrimatives( TETRAHEDRON_MESH& mesh, ARRAY<int> elements);
            
        static T 
        DistanceToElement( const TRIANGLE_3D<T>& E, const VECTOR<T,3>& p ) throw(GEOMETRY_PREDICATE_EXCEPTION);
            
        static T 
        DistanceToElement( const SEGMENT_2D<T>& E, const VECTOR<T,2>& p ) throw(GEOMETRY_PREDICATE_EXCEPTION);
            
    };
}

#endif

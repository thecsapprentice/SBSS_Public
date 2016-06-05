#include <Common_Geometry/Nonmanifold_Topology_Generation/GEOMETRY_PREDICATES.h>

using namespace PhysBAM;

//==================================================================================================
//                                            Inside Tests
//
//==================================================================================================


template<class T> bool
GEOMETRY_PREDICATES<T>::Inside( const TETRAHEDRON<T>& T1, const VECTOR<T,3>& P ) throw(GEOMETRY_PREDICATE_EXCEPTION) {
    
    bool perturb=false;
    ROBUST_SIMPLEX_INTERACTIONS<VECTOR<double,3> > robust_simplex_interactions;
    robust_simplex_interactions.Set_Tolerance((double)1e-15);
    robust_simplex_interactions.Flush_Intersection_Cache();
    
    bool inside = true;
    VECTOR<TRIANGLE_3D<T>,4> triangles;
    triangles(1)=T1.triangle1;triangles(2)=T1.triangle2;
    triangles(3)=T1.triangle3;triangles(4)=T1.triangle4;
    VECTOR<double,3> p;
    for(int k=1;k<=3;++k) p(k)=(double)P(k);
    
    for(int i=1;i<=4;++i){
        VECTOR<VECTOR<double,3>,3> tri;
        for(int v=1;v<=3;++v) tri(1)(v)=(double)triangles(i).x1(v);
        for(int v=1;v<=3;++v) tri(2)(v)=(double)triangles(i).x2(v);
        for(int v=1;v<=3;++v) tri(3)(v)=(double)triangles(i).x3(v);
        
        PAIR<double,bool> volume1=robust_simplex_interactions.Signed_Volume_Times_Six(VECTOR<VECTOR<double,3>,4>(p,
                                                                                                                 tri(1),
                                                                                                                 tri(2),
                                                                                                                 tri(3)));
        if(!volume1.y) perturb=true;
        if(volume1.x < 0 ) inside=false;
    }

    if( perturb )
        throw GEOMETRY_PREDICATE_EXCEPTION( inside, false, "Non-Robust input detected during intersection tests. Perturbation may be needed." );


    //PHYSBAM_ASSERT( T1.Inside(P) == inside );

    return inside;
}

template<class T> bool
GEOMETRY_PREDICATES<T>::Inside( const TRIANGLE_2D<T>& T1, const VECTOR<T,2>& P ) throw(GEOMETRY_PREDICATE_EXCEPTION) {

    bool perturb=false;
    ROBUST_SIMPLEX_INTERACTIONS<VECTOR<double,2> > robust_simplex_interactions;
    robust_simplex_interactions.Set_Tolerance((double)1e-15);
    robust_simplex_interactions.Flush_Intersection_Cache();
    
    VECTOR<VECTOR<double,2>, 1> p;
    for(int k=1;k<=2;++k) p(1)(k)=(double)P(k);
    VECTOR<VECTOR<double,2>,3> tri;
    for(int v=1;v<=2;++v) tri(1)(v)=(double)T1.X(1)(v);
    for(int v=1;v<=2;++v) tri(2)(v)=(double)T1.X(2)(v);
    for(int v=1;v<=2;++v) tri(3)(v)=(double)T1.X(3)(v);

    VECTOR<bool,2> intersects;
    intersects(1) = robust_simplex_interactions.Intersection(tri, p, &intersects(2));

    if(!intersects(2)) perturb=true;

    if( perturb )
        throw GEOMETRY_PREDICATE_EXCEPTION( intersects(1), false, "Non-Robust input detected during intersection tests. Perturbation may be needed." );


    //PHYSBAM_ASSERT( T1.Inside(P) == intersects(1) );

    return intersects(1);
}


template<typename T> template<int d> int
GEOMETRY_PREDICATES<T>::InsideOutside( const RANGE<VECTOR<T,d> >& box,
               const typename MATERIAL_PREDICATE_TESSELATED_VOLUME_POLICY<T,d>::T_ELEMENT& T1 ) throw(GEOMETRY_PREDICATE_EXCEPTION)
{
    int inside = 0;
    int outside = 0;
        
    for( int i = 1; i <= T1.X.m; i++){
        if( box.Boundary( T1.X(i), 0 ) )
            continue;
        if( box.Inside( T1.X(i) ) )
            inside ++;
        else
            outside ++;
    }

    if( inside && outside == 0)
        return 1;
    if( outside && inside == 0)
        return -1;
    return 0;
}


//==================================================================================================
//                                            Intersection Tests
//
//        Referenced from: http://www.geometrictools.com/Documentation/MethodOfSeparatingAxes.pdf
//==================================================================================================

template<typename T> template<int d> int
GEOMETRY_PREDICATES<T>::WhichSide( const VECTOR< VECTOR<T,d>, d+1>& points, VECTOR<T,d> D, VECTOR<T,d> P ) throw(GEOMETRY_PREDICATE_EXCEPTION){
    int positive = 0; int negative = 0;
    for( int i = 1; i <= points.m; i++){
        T t = VECTOR<T,d>::Dot_Product(D, points(i)-P);
        if( t > 0 ) positive++;
        else if( t < 0 ) negative++;
        if(positive && negative) return 0;
    }
    return (positive ? 1 : -1);
}

template<class T> bool 
GEOMETRY_PREDICATES<T>::TestIntersection( const TRIANGLE_2D<T>& T1, 
                  const SEGMENT_2D<T>& T2 ) throw(GEOMETRY_PREDICATE_EXCEPTION)
{

    VECTOR< VECTOR<double, 2>, 3 > tri;
    VECTOR< VECTOR<double, 2>, 2 > seg;
    for( int i=1; i<=3; i++) for(int v=1; v<=2; v++) tri(i)(v) = (double)T1.X(i)(v);
    seg(1)(1) = (double)T2.x1(1);     seg(1)(2) = (double)T2.x1(2);
    seg(2)(1) = (double)T2.x2(1);     seg(2)(2) = (double)T2.x2(2);

    bool perturb=false,intersection=false;
    ROBUST_SIMPLEX_INTERACTIONS<VECTOR<double,2> > robust_simplex_interactions;
    robust_simplex_interactions.Set_Tolerance((double)1e-15);
    robust_simplex_interactions.Flush_Intersection_Cache();

    VECTOR<bool,2> intersect = robust_simplex_interactions.Intersection_Test(tri, seg);

    if(intersect(1) && !intersect(2)) perturb=true;
    else if(intersect(1) && intersect(2)) perturb=false;
    if(intersect(1) && intersect(2)) intersection=true;

    if( perturb )
        throw GEOMETRY_PREDICATE_EXCEPTION( intersection, false, "Non-Robust input detected during intersection tests. Perturbation may be needed." );

    return intersection;
}


template<class T> bool 
GEOMETRY_PREDICATES<T>::TestIntersection( const TETRAHEDRON<T>& T1, 
                  const TRIANGLE_3D<T>& T2 ) throw(GEOMETRY_PREDICATE_EXCEPTION)
{

    VECTOR< VECTOR<double, 3>, 4 > tet;
    VECTOR< VECTOR<double, 3>, 3 > tri;
    for( int i=1; i<=4; i++) for(int v=1; v<=3; v++) tet(i)(v) = (double)T1.X(i)(v);
    for( int i=1; i<=3; i++) for(int v=1; v<=3; v++) tri(i)(v) = (double)T2.X(i)(v);

    bool perturb=false,intersection=false;
    ROBUST_SIMPLEX_INTERACTIONS<VECTOR<double,3> > robust_simplex_interactions;
    robust_simplex_interactions.Set_Tolerance((double)1e-15);
    robust_simplex_interactions.Flush_Intersection_Cache();

    VECTOR<bool,2> intersect;
    intersect(1) = robust_simplex_interactions.Intersection(tet, tri, &intersect(2));

    if(intersect(1) && !intersect(2)) perturb=true;
    else if(intersect(1) && intersect(2)) perturb=false;
    if(intersect(1) && intersect(2)) intersection=true;

    if( perturb )
        throw GEOMETRY_PREDICATE_EXCEPTION( intersection, false, "Non-Robust input detected during intersection tests. Perturbation may be needed." );

    return intersection;
}



template<class T> bool 
GEOMETRY_PREDICATES<T>::TestIntersection( const TRIANGLE_2D<T>& T1, 
                  const TRIANGLE_2D<T>& T2 ) throw(GEOMETRY_PREDICATE_EXCEPTION)
{
    /*
    typedef VECTOR<T,2> TV;

    for( int i0=0, i1=T1.X.m-1; i0 < T1.X.m; i1=i0, i0++){  
        TV t = (T1.X(i0+1) - T1.X(i1+1));
        TV D( t.y, -t.x);
        //LOG::cout << "D: " << D << std::endl;
        if( WhichSide( T2.X, D, T1.X(i0+1)) > 0 ){
            //LOG::cout << "Triangle 2 can be separated from Triangle 1 with this line: " << D << " " << T1.X(i0+1) << std::endl;
            return false;
        }
            
    }

    for( int i0=0, i1=T2.X.m-1; i0 < T2.X.m; i1=i0, i0++){  
        TV t = (T2.X(i0+1) - T2.X(i1+1));
        TV D( t.y, -t.x);
        //LOG::cout << "D: " << D << std::endl;
        if( WhichSide( T1.X, D, T2.X(i0+1)) > 0 ){
            //LOG::cout << "Triangle 1 can be separated from Triangle 2 with this line: " << D << " " << T1.X(i0+1) << std::endl;
            return false;
        }
    }
    return true;
    */

    bool perturb=false,intersection=false;
    ROBUST_SIMPLEX_INTERACTIONS<VECTOR<double,2> > robust_simplex_interactions;
    robust_simplex_interactions.Set_Tolerance((double)1e-15);
    robust_simplex_interactions.Flush_Intersection_Cache();

    VECTOR< TRIANGLE_2D<T>, 2 > exterior( T2, T1 );
    VECTOR< TRIANGLE_2D<T>, 2 > interior( T1, T2 );

    for( int e = 1; e<=2; e++){
        for(int i=1;i<=interior(e).X.m;++i){VECTOR<VECTOR<double,2>,3> T2_X;VECTOR<double,2> X;
            for(int j=1;j<=3;++j) for(int k=1;k<=2;++k) T2_X(j)(k)=(double)exterior(e).X(j)(k);
            for(int j=1;j<=2;++j) X(j)=(double)interior(e).X(i)(j);
            
            VECTOR<bool,2> intersect=robust_simplex_interactions.Intersection_Test(T2_X,VECTOR<VECTOR<double,2>,1>(X));
            if(intersect(1) && !intersect(2)) perturb=true;
            else if(intersect(1) && intersect(2)) perturb=false;
            if(intersect(1) && intersect(2)) intersection=true;}

        if( perturb )
            PHYSBAM_FATAL_ERROR( "Non-Robust input detected during intersection tests. Perturbation may be needed." );
        if( intersection )
            return true;

        for(int i=1;i<=interior(e).X.m;++i){int next=(i==3)?1:i+1;
            VECTOR<VECTOR<double,2>,3> T2_X;VECTOR<double,2> X,X_next;
            for(int j=1;j<=3;++j) for(int k=1;k<=2;++k) T2_X(j)(k)=(double)exterior(e).X(j)(k);
            for(int j=1;j<=2;++j){X(j)=(double)interior(e).X(i)(j);X_next(j)=(double)interior(e).X(next)(j);}
            
            VECTOR<bool,2> intersect=robust_simplex_interactions.Intersection_Test(T2_X,VECTOR<VECTOR<double,2>,2>(X,X_next));
            if(intersect(1) && !intersect(2)) perturb=true;
            else if(intersect(1) && intersect(2)) perturb=false;
            if(intersect(1) && intersect(2)) intersection=true;}
        
        if( perturb )
            PHYSBAM_FATAL_ERROR( "Non-Robust input detected during intersection tests. Perturbation may be needed." );
        if( intersection )
            return true;
    }

    return intersection;
}


template<class T> bool 
GEOMETRY_PREDICATES<T>::TestIntersection( const TETRAHEDRON<T>& T1, 
                  const TETRAHEDRON<T>& T2 ) throw(GEOMETRY_PREDICATE_EXCEPTION)
{
    /*

    typedef VECTOR<T,3> TV;
    typedef VECTOR<TRIANGLE_3D<T>, 4> FACES;
    typedef VECTOR<SEGMENT_3D<T>, 6> EDGES;
    FACES F1( T1.triangle1,T1.triangle2,T1.triangle3,T1.triangle4 );
    FACES F2( T2.triangle1,T2.triangle2,T2.triangle3,T2.triangle4 );
    EDGES E1( SEGMENT_3D<T>(T1.X(1),T1.X(2)), SEGMENT_3D<T>(T1.X(1),T1.X(3)),
              SEGMENT_3D<T>(T1.X(1),T1.X(4)), SEGMENT_3D<T>(T1.X(2),T1.X(3)),
              SEGMENT_3D<T>(T1.X(2),T1.X(4)), SEGMENT_3D<T>(T1.X(3),T1.X(4)));
    EDGES E2( SEGMENT_3D<T>(T2.X(1),T2.X(2)), SEGMENT_3D<T>(T2.X(1),T2.X(3)),
              SEGMENT_3D<T>(T2.X(1),T2.X(4)), SEGMENT_3D<T>(T2.X(2),T2.X(3)),
              SEGMENT_3D<T>(T2.X(2),T2.X(4)), SEGMENT_3D<T>(T2.X(3),T2.X(4)));

    for( int i = 0; i < F1.m; i++){
        TV D = F1(i+1).normal;
        if( WhichSide(T2.X, D, F1(i+1).x1) > 0 )
            return false;
    }

    for( int i = 0; i < F2.m; i++){
        TV D = F2(i+1).normal;
        if( WhichSide(T1.X, D, F2(i+1).x1) > 0 )
            return false;
    }

    for( int i=0; i < E1.m; i++){
        for( int j=0; j< E2.m; j++){
            TV D = TV::Cross_Product( (E1(i+1).x2 - E1(i+1).x1), (E2(j+1).x2 - E2(j+1).x1 ) );
            int side1 = WhichSide( T1.X, D, E1(i+1).x1 );
            if(side1 == 0) continue;
            int side2 = WhichSide( T2.X, D, E1(i+1).x1 );
            if(side2 == 0) continue;
            
            if( side1*side2 < 0 ) return false;
        }
    }
    return true;
    
    */

    bool perturb=false,intersection=false;
    ROBUST_SIMPLEX_INTERACTIONS<VECTOR<double,3> > robust_simplex_interactions;
    robust_simplex_interactions.Set_Tolerance((double)1e-15);
    robust_simplex_interactions.Flush_Intersection_Cache();
    VECTOR<TRIANGLE_3D<T>,4> triangles;

    VECTOR< TETRAHEDRON<T>, 2 > exterior( T1, T2 );
    VECTOR< TETRAHEDRON<T>, 2 > interior( T2, T1 );

    for( int e = 1; e<=2; e++){
        // Test for Triangle Intersections
        {
            
            triangles(1)=interior(e).triangle1;triangles(2)=interior(e).triangle2;triangles(3)=interior(e).triangle3;triangles(4)=interior(e).triangle4;
            for(int i=1;i<=4;++i){
                VECTOR<VECTOR<double,3>,4> T2_X;VECTOR<VECTOR<double,3>,3> tri;
                for(int j=1;j<=4;++j) for(int k=1;k<=3;++k) T2_X(j)(k)=(double)exterior(e).X(j)(k);
                for(int j=1;j<=3;++j) tri(1)(j)=(double)triangles(i).x1(j);
                for(int j=1;j<=3;++j) tri(2)(j)=(double)triangles(i).x2(j);
                for(int j=1;j<=3;++j) tri(3)(j)=(double)triangles(i).x3(j);
                
                VECTOR<bool,2> intersect;     
                intersect(1)=robust_simplex_interactions.Intersection(T2_X,tri,&intersect(2));
                if(intersect(1) && !intersect(2)) perturb=true;
                else if(intersect(1) && intersect(2)) perturb=(perturb || false);
                if(intersect(1) && intersect(2)) intersection=true;
            }
        }
        
        if( perturb )
            PHYSBAM_FATAL_ERROR( "Non-Robust input detected during intersection tests. Perturbation may be needed." );
        if( intersection )
            return true;

        // Test For Edge Intersections
        {
            ARRAY<VECTOR<VECTOR<T,3>,2> > edges;
            for(int i=1;i<=4;++i){VECTOR<VECTOR<T,3>,2> e1(triangles(i).x1,triangles(i).x2),e2(triangles(i).x2,triangles(i).x3),e3(triangles(i).x3,triangles(i).x1);
                edges.Append(e1);edges.Append(e2);edges.Append(e3);}
            for(int i=1;i<=edges.m;++i){VECTOR<VECTOR<double,3>,4> T2_X;VECTOR<VECTOR<double,3>,2> edge;
                for(int j=1;j<=4;++j) for(int k=1;k<=3;++k) T2_X(j)(k)=(double)exterior(e).X(j)(k);
                for(int j=1;j<=2;++j) for(int k=1;k<=3;++k) edge(j)(k)=(double)edges(i)(j)(k);
                
                VECTOR<bool,2> intersect;
                intersect(1)=robust_simplex_interactions.Intersection(T2_X,edge,&intersect(2));
                if(intersect(1) && !intersect(2)) perturb=true;
                else if(intersect(1) && intersect(2)) perturb=false;
                if(intersect(1) && intersect(2)) intersection=true;}
        }
        if( perturb )
            PHYSBAM_FATAL_ERROR( "Non-Robust input detected during intersection tests. Perturbation may be needed." );
        if( intersection )
            return true;
    }
/*
    // Test For Full Inculsions ( T1 inside T2 || T2 inside T1 )
    {

        bool inside;

        for( int e = 1; e<=2; e++){
            inside = true;
            VECTOR<TRIANGLE_3D<T>,4> triangles;
            triangles(1)=exterior(e).triangle1;triangles(2)=exterior(e).triangle2;
            triangles(3)=exterior(e).triangle3;triangles(4)=exterior(e).triangle4;
            for(int i=1;i<=4;++i){
                for(int j=1;j<=4;++j){
                    VECTOR<VECTOR<double,3>,3> tri;
                    VECTOR<double,3> p;
                    for(int k=1;k<=3;++k) p(k)=(double)interior(e).X(j)(k);
                    for(int j=1;j<=3;++j) tri(1)(j)=(double)triangles(i).x1(j);
                    for(int j=1;j<=3;++j) tri(2)(j)=(double)triangles(i).x2(j);
                    for(int j=1;j<=3;++j) tri(3)(j)=(double)triangles(i).x3(j);

                    PAIR<double,bool> volume1=robust_simplex_interactions.Signed_Volume_Times_Six(VECTOR<VECTOR<double,3>,4>(p,
                                                                                                                   tri(1),
                                                                                                                   tri(2),
                                                                                                                   tri(3)));
                    if(!volume1.y) perturb=true;
                    if(volume1.x < 0 ) inside=false;
                }
            }
            if(inside)
                break;
        }
        if(inside) intersection=true;
    }
*/
    
    return intersection;

}



template<class T> bool 
GEOMETRY_PREDICATES<T>::TestIntersection( const RANGE<VECTOR<T,2> >& box, 
                  const typename MATERIAL_PREDICATE_TESSELATED_VOLUME_POLICY<T,2>::T_ELEMENT& T1 ) throw(GEOMETRY_PREDICATE_EXCEPTION)
{
    // Intersect box by dividing it into two triangles and testing each in turn
    typedef VECTOR<int,2> TI;
    ARRAY< VECTOR<T,2>, TI > corners;
    box.Corners( corners );
    const TRIANGLE_2D<T> T_A( corners(TI(0,0)), corners(TI(1,0)), corners(TI(1,1)) );
    const TRIANGLE_2D<T> T_B( corners(TI(0,0)), corners(TI(1,1)), corners(TI(0,1)) );   

    return TestIntersection( T_A, T1 ) || TestIntersection( T_B, T1 );
}

template<class T> bool 
GEOMETRY_PREDICATES<T>::TestIntersection( const RANGE<VECTOR<T,3> >& box, 
                  const typename MATERIAL_PREDICATE_TESSELATED_VOLUME_POLICY<T,3>::T_ELEMENT& T1 ) throw(GEOMETRY_PREDICATE_EXCEPTION)
{

    VECTOR< VECTOR<double, 3>, 2 > b;
    VECTOR< VECTOR<double, 3>, 4 > tet;
    for( int i=1; i<=4; i++) for(int v=1; v<=3; v++) tet(i)(v) = (double)T1.X(i)(v);
    for(int v=1; v<=3; v++){
        b(1)(v) = (double)box.min_corner(v);
        b(2)(v) = (double)box.max_corner(v);
    }

    bool perturb=false,intersection=false;
    ROBUST_SIMPLEX_INTERACTIONS<VECTOR<double,3> > robust_simplex_interactions;
    robust_simplex_interactions.Set_Tolerance((double)1e-15);
    robust_simplex_interactions.Flush_Intersection_Cache();

    VECTOR<bool,2> intersect;
    intersect(1) = robust_simplex_interactions.Intersection(b, tet, &intersect(2));

    if(intersect(1) && !intersect(2)) perturb=true;
    else if(intersect(1) && intersect(2)) perturb=false;
    if(intersect(1) && intersect(2)) intersection=true;

    if( perturb )
        throw GEOMETRY_PREDICATE_EXCEPTION( intersection, false, "Non-Robust input detected during intersection tests. Perturbation may be needed." );


    return intersection;
}


template<class T> bool 
GEOMETRY_PREDICATES<T>::TestIntersection( const VECTOR< VECTOR<T,3>, 8 >& hex, 
                  const typename MATERIAL_PREDICATE_TESSELATED_VOLUME_POLICY<T,3>::T_ELEMENT& T1 ) throw(GEOMETRY_PREDICATE_EXCEPTION)
{
    // Intersect box by dividing it into five tetrahedra and testing each in turn
    //typedef VECTOR<int,3> TI;
    //ARRAY< VECTOR<T,3> , TI > corners;
    //box.Corners( corners );
    const TETRAHEDRON<T> T_A( hex(2), hex(1), hex(6), hex(4) );
    const TETRAHEDRON<T> T_B( hex(1), hex(5), hex(6), hex(7) );
    const TETRAHEDRON<T> T_C( hex(1), hex(6), hex(4), hex(7) );
    const TETRAHEDRON<T> T_D( hex(1), hex(4), hex(3), hex(7) );
    const TETRAHEDRON<T> T_E( hex(6), hex(4), hex(7), hex(8) );
    bool intersects = (TestIntersection( T_A, T1 ) || TestIntersection( T_B, T1 ) ||
                       TestIntersection( T_C, T1 ) || TestIntersection( T_D, T1 ) ||
                       TestIntersection( T_E, T1 ) );

    return intersects;
}


template<class T> bool 
GEOMETRY_PREDICATES<T>::TestIntersection( const RANGE<VECTOR<T,2> >& box, 
                  const typename MATERIAL_PREDICATE_TESSELATED_VOLUME_POLICY<T,2>::T_BOUNDARY_ELEMENT& T1 ) throw(GEOMETRY_PREDICATE_EXCEPTION)
{
    //bool old_intersect =  INTERSECTION::Intersects( box, T1 );

    // Intersect box by dividing it into two triangles and testing each in turn
    typedef VECTOR<int,2> TI;
    ARRAY< VECTOR<T,2>, TI > corners;
    box.Corners( corners );
    const TRIANGLE_2D<T> T_A( corners(TI(0,0)), corners(TI(1,0)), corners(TI(1,1)) );
    const TRIANGLE_2D<T> T_B( corners(TI(0,0)), corners(TI(1,1)), corners(TI(0,1)) );  

    bool intersects = (TestIntersection( T_A, T1 ) || TestIntersection( T_B, T1 ));
    //PHYSBAM_ASSERT( old_intersect == intersects );

    return intersects;
}

template<class T> bool 
GEOMETRY_PREDICATES<T>::TestIntersection( const RANGE<VECTOR<T,3> >& box, 
                  const typename MATERIAL_PREDICATE_TESSELATED_VOLUME_POLICY<T,3>::T_BOUNDARY_ELEMENT& T1 ) throw(GEOMETRY_PREDICATE_EXCEPTION)
{
    //bool old_intersect =  INTERSECTION::Intersects( box, T1 );

    // Intersect box by dividing it into five tetrahedra and testing each in turn
    typedef VECTOR<int,3> TI;
    ARRAY< VECTOR<T,3> , TI > corners;
    box.Corners( corners );
    const TETRAHEDRON<T> T_A( corners(TI(0,0,1)), corners(TI(0,0,0)), corners(TI(1,0,1)), corners(TI(0,1,1)) );
    const TETRAHEDRON<T> T_B( corners(TI(0,0,0)), corners(TI(1,0,0)), corners(TI(1,0,1)), corners(TI(1,1,0)) );
    const TETRAHEDRON<T> T_C( corners(TI(0,0,0)), corners(TI(1,0,1)), corners(TI(0,1,1)), corners(TI(1,1,0)) );
    const TETRAHEDRON<T> T_D( corners(TI(0,0,0)), corners(TI(0,1,1)), corners(TI(0,1,0)), corners(TI(1,1,0)) );
    const TETRAHEDRON<T> T_E( corners(TI(1,0,1)), corners(TI(0,1,1)), corners(TI(1,1,0)), corners(TI(1,1,1)) );
    bool intersects = (TestIntersection( T_A, T1 ) || TestIntersection( T_B, T1 ) ||
                       TestIntersection( T_C, T1 ) || TestIntersection( T_D, T1 ) ||
                       TestIntersection( T_E, T1 ) );
    //PHYSBAM_ASSERT( old_intersect == intersects );

    return intersects;
}

//=====================================================================================
//
//                              Predicate Method Helpers
//
//=====================================================================================

template< class T> bool 
GEOMETRY_PREDICATES<T>::Element_Neighbors( const TETRAHEDRON_MESH& mesh, int T1, int T2, VECTOR<int,3>& edge_nodes ){
    bool neighbors = mesh.Face_Neighbors( T1, T2 );
    if(!neighbors)
        return false;
    int i,j,k,l;mesh.elements(T1).Get(i,j,k,l);

    if(mesh.Triangle_In_Tetrahedron(j,k,l,T2)){
        edge_nodes = VECTOR<int,3>(j,k,l); return true;}
    if(mesh.Triangle_In_Tetrahedron(i,k,l,T2)){
        edge_nodes = VECTOR<int,3>(i,k,l); return true;}
    if(mesh.Triangle_In_Tetrahedron(i,j,l,T2)){
        edge_nodes = VECTOR<int,3>(i,j,l); return true;}
    if(mesh.Triangle_In_Tetrahedron(i,j,k,T2)){
        edge_nodes = VECTOR<int,3>(i,j,k); return true;}
    PHYSBAM_FATAL_ERROR();
}

template< class T> bool 
GEOMETRY_PREDICATES<T>::Element_Neighbors( const TRIANGLE_MESH& mesh, int T1, int T2, VECTOR<int,2>& edge_nodes ){
    bool neighbors = mesh.Edge_Neighbors( T1, T2 );
    if(!neighbors)
        return false;
    int i,j,k;mesh.elements(T1).Get(i,j,k);

    if(mesh.Segment_In_Triangle(i,j,T2)){
        edge_nodes = VECTOR<int,2>(i,j); return true;}
    if(mesh.Segment_In_Triangle(j,k,T2)){
        edge_nodes = VECTOR<int,2>(j,k); return true;}
    if(mesh.Segment_In_Triangle(k,i,T2)){
        edge_nodes = VECTOR<int,2>(k,i); return true;}
    PHYSBAM_FATAL_ERROR();
}

template< class T> bool 
GEOMETRY_PREDICATES<T>::Element_Neighbors( const SEGMENT_MESH& mesh, int T1, int T2, VECTOR<int,1>& edge_nodes ){
    bool neighbors = mesh.Segments_Adjacent( T1, T2 );
    if(!neighbors)
        return false;
    int i,j;mesh.elements(T1).Get(i,j);

    if(mesh.Node_In_Segment(i,T2)){
        edge_nodes = VECTOR<int,1>(i); return true;}
    if(mesh.Node_In_Segment(j,T2)){
        edge_nodes = VECTOR<int,1>(j); return true;}
    PHYSBAM_FATAL_ERROR();
}


template< class T> HASHTABLE<int, ARRAY< VECTOR<int,2> > > 
GEOMETRY_PREDICATES<T>::GenerateBoundaryPrimatives( TRIANGLE_MESH& mesh, ARRAY<int> elements){

    HASHTABLE<int, ARRAY< VECTOR<int, 2> > > boundary_primatives;
    ARRAY<int> taf;

    bool incident_elements_defined=mesh.incident_elements!=0;
    if(!mesh.incident_elements) mesh.Initialize_Incident_Elements();
    mesh.Initialize_Boundary_Nodes();
    mesh.Initialize_Node_On_Boundary();

    for(int t=1;t<=elements.m;t++){
        int e_id = elements(t);
        boundary_primatives.Set( e_id, ARRAY< VECTOR<int, 2> >() );
        int i,j,k;mesh.elements(e_id).Get(i,j,k);
        //if( (*mesh.node_on_boundary)(i) ) LOG::cout << "Node " << i << " is boundary." <<std::endl;
        //if( (*mesh.node_on_boundary)(j) ) LOG::cout << "Node " << j << " is boundary." <<std::endl;
        //if( (*mesh.node_on_boundary)(k) ) LOG::cout << "Node " << k << " is boundary." <<std::endl;
        taf.Clean_Memory();
        if(mesh.Triangles_Across_Edge(e_id,i,k,taf) == 0)
            boundary_primatives.Get(e_id).Append(VECTOR<int,2>(i,k));
        taf.Clean_Memory();
        if(mesh.Triangles_Across_Edge(e_id,j,i,taf) == 0)
            boundary_primatives.Get(e_id).Append(VECTOR<int,2>(j,i));
        taf.Clean_Memory();
        if(mesh.Triangles_Across_Edge(e_id,k,j,taf) == 0)
            boundary_primatives.Get(e_id).Append(VECTOR<int,2>(k,j));
        //LOG::cout << boundary_primatives.Get(e_id).m << std::endl;
        //LOG::cout << taf.m << std::endl;
        //LOG::cout << "--------" << std::endl;
    }
    if(!incident_elements_defined){delete mesh.incident_elements;mesh.incident_elements=0;}
    return boundary_primatives;
}

template< class T> HASHTABLE<int, ARRAY< VECTOR<int,3> > > 
GEOMETRY_PREDICATES<T>::GenerateBoundaryPrimatives( TETRAHEDRON_MESH& mesh, ARRAY<int> elements){

    HASHTABLE<int, ARRAY< VECTOR<int, 3> > > boundary_primatives;
    ARRAY<int> taf;

    bool incident_elements_defined=mesh.incident_elements!=0;
    if(!mesh.incident_elements) mesh.Initialize_Incident_Elements();
    for(int t=1;t<=elements.m;t++){
        boundary_primatives.Set( elements(t), ARRAY< VECTOR<int, 3> >() );
        int i,j,k,l;mesh.elements(elements(t)).Get(i,j,k,l);
        taf.Clean_Memory();
        if(mesh.Tetrahedrons_Across_Face(elements(t),i,k,j,taf) == 0)
            boundary_primatives.Get(elements(t)).Append(VECTOR<int,3>(i,k,j));
        taf.Clean_Memory();
        if(mesh.Tetrahedrons_Across_Face(elements(t),i,j,l,taf) == 0)
            boundary_primatives.Get(elements(t)).Append(VECTOR<int,3>(i,j,l));
        taf.Clean_Memory();
        if(mesh.Tetrahedrons_Across_Face(elements(t),i,l,k,taf) == 0)
            boundary_primatives.Get(elements(t)).Append(VECTOR<int,3>(i,l,k));
        taf.Clean_Memory();
        if(mesh.Tetrahedrons_Across_Face(elements(t),j,k,l,taf) == 0)
            boundary_primatives.Get(elements(t)).Append(VECTOR<int,3>(j,k,l));
    }
    if(!incident_elements_defined){delete mesh.incident_elements;mesh.incident_elements=0;}
    return boundary_primatives;
}

template< class T> T 
GEOMETRY_PREDICATES<T>::DistanceToElement( const TRIANGLE_3D<T>& E, const VECTOR<T,3>& p ) throw(GEOMETRY_PREDICATE_EXCEPTION) {
    return E.Distance_To_Triangle(p);
}

template< class T> T 
GEOMETRY_PREDICATES<T>::DistanceToElement( const SEGMENT_2D<T>& E, const VECTOR<T,2>& p ) throw(GEOMETRY_PREDICATE_EXCEPTION) {
    //LOG::cout << "Computing Distance from " << p << " to line: " << std::endl;
    //LOG::cout << "\t" << E.x1 << "   " << E.x2 << std::endl;
    return E.Distance_From_Point_To_Segment(p);
}

template struct GEOMETRY_PREDICATES<float>;
template struct GEOMETRY_PREDICATES<double>;


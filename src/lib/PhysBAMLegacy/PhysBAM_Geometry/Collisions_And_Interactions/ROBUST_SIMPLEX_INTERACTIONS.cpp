//#####################################################################
// Copyright 2006-2007, Andrew Selle, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ROBUST_SIMPLEX_INTERACTIONS
//##################################################################### 
#include <PhysBAM_Geometry/Collisions_And_Interactions/ROBUST_SIMPLEX_INTERACTIONS.h>
#include <PhysBAM_Tools/Data_Structures/PAIR.h>
#include <PhysBAM_Tools/Math_Tools/Robust_Arithmetic.h>
#include <PhysBAM_Tools/Math_Tools/permutation.h>
#include <PhysBAM_Tools/Matrices/MATRIX_3X3.h>
#include <PhysBAM_Tools/Log/DEBUG_UTILITIES.h>
#include <PhysBAM_Tools/Log/LOG.h>

using namespace PhysBAM;
//#####################################################################
// Function Triple_Product
//#####################################################################
template<class T> VECTOR<T,2> ROBUST_SIMPLEX_INTERACTIONS<VECTOR<T,3> >::
Triple_Product(const VECTOR<TV,3>& locations)
{
    VECTOR<T,2> result; // (non_negative_sum,negative_sum)
    for(int i=1;i<=6;i++) if(triple_products.Get(permute_three(locations,i),result)){if(!permutation_of_three_is_even(i)) exchange(result.x,result.y);return result;}
    T product;
    product=locations.x.x*locations.y.y*locations.z.z;if(product>=0) result.x+=product;else result.y-=product;
    product=locations.y.x*locations.z.y*locations.x.z;if(product>=0) result.x+=product;else result.y-=product;
    product=locations.z.x*locations.x.y*locations.y.z;if(product>=0) result.x+=product;else result.y-=product;
    product=locations.x.x*locations.z.y*locations.y.z;if(product<0) result.x-=product;else result.y+=product;
    product=locations.y.x*locations.x.y*locations.z.z;if(product<0) result.x-=product;else result.y+=product;
    product=locations.z.x*locations.y.y*locations.x.z;if(product<0) result.x-=product;else result.y+=product;
    triple_products.Insert(locations,result);
    return result;
}
//#####################################################################
// Function Precompute_Triple_Products
//#####################################################################
namespace{
template<class T,int Np> void
Precompute_Triple_Products(const VECTOR<VECTOR<T,3>,Np>& locations,VECTOR<VECTOR<VECTOR<VECTOR<T,2>,Np>,Np>,Np>& triple_products)
{
    for(int i=1;i<=Np-2;i++)
        for(int j=i+1;j<=Np-1;j++)
            for(int k=j+1;k<=Np;k++){
                VECTOR<T,2> result;T product;
                product=locations(i).x*locations(j).y*locations(k).z;if(product>=0) result.x+=product;else result.y-=product;
                product=locations(j).x*locations(k).y*locations(i).z;if(product>=0) result.x+=product;else result.y-=product;
                product=locations(k).x*locations(i).y*locations(j).z;if(product>=0) result.x+=product;else result.y-=product;
                product=locations(i).x*locations(k).y*locations(j).z;if(product<0) result.x-=product;else result.y+=product;
                product=locations(j).x*locations(i).y*locations(k).z;if(product<0) result.x-=product;else result.y+=product;
                product=locations(k).x*locations(j).y*locations(i).z;if(product<0) result.x-=product;else result.y+=product;
                triple_products(i)(j)(k)=result;}
     return;
}
}
//#####################################################################
// Function Signed_Volume_Times_Six
//#####################################################################
template<class T> PAIR<T,bool> ROBUST_SIMPLEX_INTERACTIONS<VECTOR<T,3> >::
Signed_Volume_Times_Six(const VECTOR<TV,4>& locations)
{
    PAIR<T,bool> result;
    for(int i=1;i<=24;i++) if(signed_volumes_times_six.Get(permute_four(locations,i),result)){if(!permutation_of_four_is_even(i)) result.x=-result.x;return result;}
    VECTOR<T,2> determinant=Triple_Product(VECTOR<TV,3>(locations[2],locations[3],locations[4]))+Triple_Product(VECTOR<TV,3>(locations[1],locations[4],locations[3]))+
        Triple_Product(VECTOR<TV,3>(locations[1],locations[2],locations[4]))+Triple_Product(VECTOR<TV,3>(locations[1],locations[3],locations[2]));
    result.x=TV::Triple_Product(locations(2)-locations(1),locations(3)-locations(1),locations(4)-locations(1));
    if(abs(determinant.x-determinant.y)<tolerance*(determinant.x+determinant.y)) result.y=false;
    else{result.y=true;
        bool positive1=(determinant.x-determinant.y>0),positive2=(result.x>0);
        if(positive1^positive2) PHYSBAM_FATAL_ERROR();}
    signed_volumes_times_six.Insert(locations,result);
    return result;
}
//#####################################################################
// Function Signed_Volume_From_Precomputed_Triples
//#####################################################################
namespace{
template<class T,int Np> PAIR<T,bool>
    Signed_Volume_From_Precomputed_Triples(const VECTOR<VECTOR<T,3>,Np>& locations,const VECTOR<VECTOR<VECTOR<VECTOR<T,2>,Np>,Np>,Np>& triple_products,const T tolerance,int i,int j,int k,int l)
{
    // Start with bubble sort
    T sign=1.0;
    if(i>j) {exchange(i,j);sign=-sign;}if(j>k) {exchange(j,k);sign=-sign;}if(k>l) {exchange(k,l);sign=-sign;}
    if(i>j) {exchange(i,j);sign=-sign;}if(j>k) {exchange(j,k);sign=-sign;}if(i>j) {exchange(i,j);sign=-sign;}

    PAIR<T,bool> result;
    VECTOR<T,2> determinant=triple_products(j)(k)(l)+triple_products(i)(k)(l).Reversed()+triple_products(i)(j)(l)+triple_products(i)(j)(k).Reversed();
    if(sign<0) determinant=determinant.Reversed();
    result.x=sign*VECTOR<T,3>::Triple_Product(locations(j)-locations(i),locations(k)-locations(i),locations(l)-locations(i));
    if(abs(determinant.x-determinant.y)<tolerance*(determinant.x+determinant.y)) result.y=false;
    else{result.y=true;
        bool positive1=(determinant.x-determinant.y>0),positive2=(result.x>0);
        if(positive1^positive2) PHYSBAM_FATAL_ERROR();}
    return result;
}
}
//#####################################################################
// Function Triangle_Segment_Intersection_Weights
//#####################################################################
template<class T> void ROBUST_SIMPLEX_INTERACTIONS<VECTOR<T,3> >::    
Triangle_Segment_Intersection_Weights(const VECTOR<TV,3>& triangle,const VECTOR<TV,2>& segment,VECTOR<T,2>& triangle_weights,T& segment_weight,bool *is_robust_input)
{
    MATRIX<T,3> matrix(triangle(1)-triangle(3),triangle(2)-triangle(3),segment(2)-segment(1));
    VECTOR<T,3> weights=matrix.Robust_Solve_Linear_System(segment(2)-triangle(3));
    for(int i=1;i<=2;i++) triangle_weights(i)=weights(i);segment_weight=weights(3);
}
//#####################################################################
// Function Intersection
//#####################################################################
template<class T> bool ROBUST_SIMPLEX_INTERACTIONS<VECTOR<T,3> >::
Intersection(const VECTOR<TV,4>& tetrahedron,const VECTOR<TV,2>& segment,bool *is_robust_input)
{
    bool is_robust=true;
    // Separator plane must contain two tetrahedron vertices and one segment vertex
    for(int i=1;i<=3;i++) for(int j=i+1;j<=4;j++) for(int k=1;k<=2;k++){
        PAIR<T,bool> volume1=Signed_Volume_Times_Six(VECTOR<TV,4>(tetrahedron(i),tetrahedron(j),segment(k),segment(3-k)));
        if(!volume1.y){is_robust=false;continue;}
        for(int l=1;l<=4;l++) if(l!=i && l!=j){
            PAIR<T,bool> volume2=Signed_Volume_Times_Six(VECTOR<TV,4>(tetrahedron(i),tetrahedron(j),segment(k),tetrahedron(l)));
            if(!volume2.y){is_robust=false;goto Next_Plane;}
            if((volume1.x>0)^(volume2.x<0)) goto Next_Plane;}
        if(is_robust_input) *is_robust_input=true;return false;
        Next_Plane:;}
    // No separator plane found, objects are intersecting
    if(is_robust_input) *is_robust_input=is_robust;
    return true;
}
//#####################################################################
// Function Intersection
//#####################################################################
template<class T> bool ROBUST_SIMPLEX_INTERACTIONS<VECTOR<T,3> >::
Intersection(const VECTOR<TV,4>& tetrahedron,const VECTOR<TV,3>& triangle,bool *is_robust_input)
{
    bool is_robust=true;
    PAIR<T,bool> volume1,volume2;
    // Case 1 : Entire tetrahedron on same half-space of triangle plane
    volume1=Signed_Volume_Times_Six(VECTOR<TV,4>(triangle(1),triangle(2),triangle(3),tetrahedron(1)));
    if(!volume1.y) is_robust=false;
    else for(int i=2;i<=4;i++){
        volume2=Signed_Volume_Times_Six(VECTOR<TV,4>(triangle(1),triangle(2),triangle(3),tetrahedron(i)));if(!volume2.y){is_robust=false;break;}
        if((volume1.x>0)^(volume2.x>0)) break;        
        if(i==4){
            if(is_robust_input) *is_robust_input=true;            
            //LOG::cout << "INTERSECTION FAILED IN CASE 1" << std::endl;
            return false;
        }
        }
    // Case 2 : Separator plane contains two triangle vertices
    for(int i=1;i<=2;i++) for(int j=i+1;j<=3;j++) for(int k=1;k<=4;k++){
        volume1=Signed_Volume_Times_Six(VECTOR<TV,4>(triangle(i),triangle(j),tetrahedron(k),triangle(6-i-j)));if(!volume1.y){is_robust=false;continue;}
        for(int l=1;l<=4;l++) if(l!=k){
            volume2=Signed_Volume_Times_Six(VECTOR<TV,4>(triangle(i),triangle(j),tetrahedron(k),tetrahedron(l)));
            if(!volume2.y){is_robust=false;goto Next_Plane1;}
            if((volume1.x>0)^(volume2.x<0)) goto Next_Plane1;}
        if(is_robust_input) *is_robust_input=true;
        //LOG::cout << "INTERSECTION FAILED IN CASE 2" << std::endl;
        return false;
        Next_Plane1:;}
    // Case 3 : Separator plane contains two tetraherdon vertices
    for(int i=1;i<=3;i++) for(int j=i+1;j<=4;j++) for(int k=1;k<=3;k++){
        VECTOR<int,2> indices=VECTOR<int,4>(1,2,3,4).Remove_Index(j).Remove_Index(i);
        volume1=Signed_Volume_Times_Six(VECTOR<TV,4>(tetrahedron(i),tetrahedron(j),triangle(k),tetrahedron(indices(1))));if(!volume1.y){is_robust=false;continue;}
        volume2=Signed_Volume_Times_Six(VECTOR<TV,4>(tetrahedron(i),tetrahedron(j),triangle(k),tetrahedron(indices(2))));if(!volume2.y){is_robust=false;continue;}
        if((volume1.x>0)^(volume2.x>0)) continue;
        for(int l=1;l<=3;l++) if(l!=k){
            volume2=Signed_Volume_Times_Six(VECTOR<TV,4>(tetrahedron(i),tetrahedron(j),triangle(k),triangle(l)));if(!volume2.y){is_robust=false;goto Next_Plane2;}
            if((volume1.x>0)^(volume2.x<0)) goto Next_Plane2;}
        if(is_robust_input) *is_robust_input=true;
        //LOG::cout << "INTERSECTION FAILED IN CASE 3" << std::endl;
        return false;
        Next_Plane2:;}
    // No separator plane found, objects are intersecting
    if(is_robust_input) *is_robust_input=is_robust;
    return true;
}
//#####################################################################
// Function Intersection
//#####################################################################
template<class T> bool ROBUST_SIMPLEX_INTERACTIONS<VECTOR<T,3> >::
Intersection(const VECTOR<TV,2>& box,const VECTOR<TV,4>& tetrahedron,bool *is_robust_input)
{
    VECTOR<TV,12> locations;
    locations(1)=TV(box(1).x,box(1).y,box(1).z);locations(2)=TV(box(1).x,box(1).y,box(2).z);locations(3)=TV(box(1).x,box(2).y,box(1).z);locations(4)=TV(box(1).x,box(2).y,box(2).z);
    locations(5)=TV(box(2).x,box(1).y,box(1).z);locations(6)=TV(box(2).x,box(1).y,box(2).z);locations(7)=TV(box(2).x,box(2).y,box(1).z);locations(8)=TV(box(2).x,box(2).y,box(2).z);
    locations(9)=tetrahedron(1);locations(10)=tetrahedron(2);locations(11)=tetrahedron(3);locations(12)=tetrahedron(4);
    VECTOR<VECTOR<VECTOR<VECTOR<T,2>,12>,12>,12> triple_products;
    Precompute_Triple_Products(locations,triple_products);

    bool is_robust=true;

    // Case 1 : Separator plane contains one hexahedron edge and one tetrahedron vertex
    for(int i=1;i<=7;i++) for(int j=i+1;j<=8;j++) if(((i-1)^(j-1))==1||((i-1)^(j-1))==2||((i-1)^(j-1))==4) // This if-clause checks if the vertices form an edge; 12 possibilities
        for(int k=1;k<=4;k++){
            ARRAY<PAIR<T,bool> > hex_volumes;
            for(int l=1;l<=8;l++) if(l!=i&&l!=j)
                hex_volumes.Append(Signed_Volume_From_Precomputed_Triples(locations,triple_products,tolerance,i,j,k+8,l));
            ARRAY<PAIR<T,bool> > tet_volumes;
            for(int l=1;l<=4;l++) if(l!=k)
                tet_volumes.Append(Signed_Volume_From_Precomputed_Triples(locations,triple_products,tolerance,i,j,k+8,l+8));

            bool hex_robust_positive=false,hex_robust_negative=false,hex_nonrobust=false;
            for(int l=1;l<=6;l++){
                if(hex_volumes(l).y)
                    if(hex_volumes(l).x>0) hex_robust_positive=true; else hex_robust_negative=true;
                else
                    hex_nonrobust=true;}
            if(hex_robust_positive&&hex_robust_negative) continue;
            bool tet_robust_positive=false,tet_robust_negative=false,tet_nonrobust=false;
            for(int l=1;l<=3;l++){
                if(tet_volumes(l).y)
                    if(tet_volumes(l).x>0) tet_robust_positive=true; else tet_robust_negative=true;
                else
                    tet_nonrobust=true;}
            if(tet_robust_positive&&tet_robust_negative) continue;
            if(hex_nonrobust||tet_nonrobust){is_robust=false;continue;} 

            if((hex_volumes(1).x>0)^(tet_volumes(1).x<0)) continue;
            if(is_robust_input) *is_robust_input=true;
            return false;}

    // Case 2 : Separator plane contains one hexahedron vertex and one tetrahedron edge
    for(int i=1;i<=8;i++) for(int j=1;j<=3;j++) for(int k=j+1;k<=4;k++){ // This if-clause checks if the vertices form an edge; 12 possibilities
        ARRAY<PAIR<T,bool> > hex_volumes;
        for(int l=1;l<=8;l++) if(l!=i)
            hex_volumes.Append(Signed_Volume_From_Precomputed_Triples(locations,triple_products,tolerance,i,j+8,k+8,l));
        ARRAY<PAIR<T,bool> > tet_volumes;
        for(int l=1;l<=4;l++) if(l!=j&&l!=k)
            tet_volumes.Append(Signed_Volume_From_Precomputed_Triples(locations,triple_products,tolerance,i,j+8,k+8,l+8));

        bool hex_robust_positive=false,hex_robust_negative=false,hex_nonrobust=false;
        for(int l=1;l<=7;l++){
            if(hex_volumes(l).y)
                if(hex_volumes(l).x>0) hex_robust_positive=true; else hex_robust_negative=true;
            else
                hex_nonrobust=true;}
        if(hex_robust_positive&&hex_robust_negative) continue;
        bool tet_robust_positive=false,tet_robust_negative=false,tet_nonrobust=false;
        for(int l=1;l<=2;l++){
            if(tet_volumes(l).y)
                if(tet_volumes(l).x>0) tet_robust_positive=true; else tet_robust_negative=true;
            else
                tet_nonrobust=true;}
        if(tet_robust_positive&&tet_robust_negative) continue;            
        if(hex_nonrobust||tet_nonrobust){is_robust=false;continue;} 

        if((hex_volumes(1).x>0)^(tet_volumes(1).x<0)) continue;
        if(is_robust_input) *is_robust_input=true;
        return false;}

    // No separator plane found, objects are intersecting
    if(is_robust_input) *is_robust_input=is_robust;
    return true;
}
//#####################################################################
// Function Cross_Product
//#####################################################################
template<class T> VECTOR<T,2> ROBUST_SIMPLEX_INTERACTIONS<VECTOR<T,2> >::
Cross_Product(const VECTOR<TV,2>& locations)
{
    VECTOR<T,2> result; // (non_negative_sum,negative_sum)
    for(int i=1;i<=2;i++) if(cross_products.Get(permute_two(locations,i),result)){if(!permutation_of_two_is_even(i)) exchange(result.x,result.y);return result;}
    T product;
    product=locations.x.x*locations.y.y;if(product>=0) result.x+=product;else result.y-=product;
    product=locations.x.y*locations.y.x;if(product<0) result.x-=product;else result.y+=product;
    cross_products.Insert(locations,result);
    return result;
}
//#####################################################################
// Function Signed_Area_Times_Two
//#####################################################################
template<class T> PAIR<T,bool> ROBUST_SIMPLEX_INTERACTIONS<VECTOR<T,2> >::
Signed_Area_Times_Two(const VECTOR<TV,3>& locations)
{
    PAIR<T,bool> result;
    for(int i=1;i<=6;i++) if(signed_areas_times_two.Get(permute_three(locations,i),result)){
        if(!permutation_of_three_is_even(i)) result.x=-result.x;return result;}
    VECTOR<T,2> determinant=Cross_Product(VECTOR<TV,2>(locations[2],locations[3]))+Cross_Product(VECTOR<TV,2>(locations[3],locations[1]))+Cross_Product(VECTOR<TV,2>(locations[1],locations[2]));
    result.x=TV::Cross_Product(locations[2]-locations[1],locations[3]-locations[1]).x;
    if(abs(determinant.x-determinant.y)<tolerance*(determinant.x+determinant.y)) result.y=false;
    else{result.y=true;
        bool positive1=(determinant.x-determinant.y>0),positive2=(result.x>0);
        if(positive1^positive2) PHYSBAM_FATAL_ERROR();}
    signed_areas_times_two.Insert(locations,result);
    return result;
}
//#####################################################################
// Function Intersection
//#####################################################################
template<class T> bool ROBUST_SIMPLEX_INTERACTIONS<VECTOR<T,2> >::
Intersection(const VECTOR<TV,3>& triangle,const VECTOR<TV,2>& segment,bool *is_robust_input)
{
    bool is_robust=true;
    PAIR<T,bool> area1,area2;
    // Case 1 : Entire triangle on same half-space of segment line
    area1=Signed_Area_Times_Two(VECTOR<TV,3>(segment(1),segment(2),triangle(1)));
    if(!area1.y) is_robust=false;
    else for(int i=2;i<=3;i++){
        area2=Signed_Area_Times_Two(VECTOR<TV,3>(segment(1),segment(2),triangle(i)));if(!area2.y){is_robust=false;break;}
        if((area1.x>0)^(area2.x>0)) break;
        if(i==4){if(is_robust_input) *is_robust_input=true;return false;}}
    // Case 2 : try combinations of one point on triangle and one point on segment
    for(int segment_i=1;segment_i<=2;segment_i++) for(int triangle_i=1;triangle_i<=3;triangle_i++){
        VECTOR<int,2> other_triangle_indices(triangle_i%3+1,(triangle_i%3+1)%3+1);
        int other_segment_i=3-segment_i;
        area1=Signed_Area_Times_Two(VECTOR<TV,3>(segment(segment_i),triangle(triangle_i),segment(other_segment_i)));if(!area1.y){is_robust=false;continue;}
        for(int dummy=1;dummy<=2;dummy++){int other_triangle_i=other_triangle_indices(dummy);
            area2=Signed_Area_Times_Two(VECTOR<TV,3>(segment(segment_i),triangle(triangle_i),triangle(other_triangle_i)));if(!area2.y){is_robust=false;goto Next_Line;}
            if((area1.x>0)^(area2.x<0)) goto Next_Line;}
        if(is_robust_input) *is_robust_input=true;return false;
      Next_Line:;}
    if(is_robust_input) *is_robust_input=is_robust;
    return true;
}
//#####################################################################
// Function Intersection
//#####################################################################
template<class T> bool ROBUST_SIMPLEX_INTERACTIONS<VECTOR<T,2> >::
Intersection(const VECTOR<TV,3>& triangle,const VECTOR<TV,1>& point,bool *is_robust_input)
{
    bool is_robust=true;
    PAIR<T,bool> area1,area2;
    // Iterate over all edges
    for(int triangle_1=1;triangle_1<=2;triangle_1++) for(int triangle_2=triangle_1+1;triangle_2<=3;triangle_2++){
        area1=Signed_Area_Times_Two(VECTOR<TV,3>(triangle(triangle_1),triangle(triangle_2),point(1)));if(!area1.y){is_robust=false;continue;}
        int other_triangle_index=VECTOR<int,3>(1,2,3).Remove_Index(triangle_2).Remove_Index(triangle_1)[1];
        area2=Signed_Area_Times_Two(VECTOR<TV,3>(triangle(triangle_1),triangle(triangle_2),triangle(other_triangle_index)));if(!area2.y){is_robust=false;continue;}
        if((area1.x>0)^(area2.x<0)) continue;
        if(is_robust_input) *is_robust_input=true;return false;}
    if(is_robust_input) *is_robust_input=is_robust;
    return true;
}
//#####################################################################
// Function Intersection_Test
//#####################################################################
template<class T> VECTOR<bool,2> ROBUST_SIMPLEX_INTERACTIONS<VECTOR<T,2> >::
Intersection_Test(const VECTOR<TV,3>& triangle,const VECTOR<TV,2>& segment)
{
    bool robust;
    bool value=Intersection(triangle,segment,&robust);
    return VECTOR<bool,2>(value,robust);
}
//#####################################################################
// Function Intersection_Test
//#####################################################################
template<class T> VECTOR<bool,2> ROBUST_SIMPLEX_INTERACTIONS<VECTOR<T,2> >::
Intersection_Test(const VECTOR<TV,3>& triangle,const VECTOR<TV,1>& point)
{
    bool robust;
    bool value=Intersection(triangle,point,&robust);
    return VECTOR<bool,2>(value,robust);
}
//####################################################################
template class ROBUST_SIMPLEX_INTERACTIONS<VECTOR<float,2> >;
template class ROBUST_SIMPLEX_INTERACTIONS<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class ROBUST_SIMPLEX_INTERACTIONS<VECTOR<double,2> >;
template class ROBUST_SIMPLEX_INTERACTIONS<VECTOR<double,3> >;
#endif

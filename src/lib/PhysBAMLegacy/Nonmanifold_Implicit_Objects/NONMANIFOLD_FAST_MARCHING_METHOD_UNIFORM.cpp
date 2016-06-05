//#####################################################################
// Copyright 2005-2007, Eran Guendelman, Geoffrey Irving, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Tools/Parallel_Computation/DOMAIN_ITERATOR_THREADED.h>
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#endif
#include <PhysBAM_Geometry/Basic_Geometry/RAY.h>
#include <PhysBAM_Geometry/Grids_Uniform_Collisions/GRID_BASED_COLLISION_GEOMETRY_UNIFORM.h>
#include <Nonmanifold_Implicit_Objects/NONMANIFOLD_FAST_MARCHING_METHOD_UNIFORM.h>
#include <Nonmanifold_Implicit_Objects/NONMANIFOLD_LEVELSET_2D.h>
#include <Nonmanifold_Implicit_Objects/NONMANIFOLD_LEVELSET_3D.h>
#include <PhysBAM_Geometry/Level_Sets/FAST_MARCHING.h>
#include <PhysBAM_Tools/Math_Tools/Hash.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T, int d> NONMANIFOLD_FAST_MARCHING_METHOD_UNIFORM<T,d>::
NONMANIFOLD_FAST_MARCHING_METHOD_UNIFORM(const T_NONMANIFOLD_LEVELSET& nonmanifold_levelset)
:nonmanifold_levelset(nonmanifold_levelset),stopping_distance(0),cell_mesh(nonmanifold_levelset.mesh),
 domain_range(nonmanifold_levelset.mesh.Node_Domain()),
 flat_range(nonmanifold_levelset.mesh.Node_Mesh_Counts())
{
//    cell_grid=nonmanifold_levelset.grid.Is_MAC_Grid()?nonmanifold_levelset.grid:nonmanifold_levelset.grid.Get_MAC_Grid_At_Regular_Positions();
//    RANGE<TV_INT> domain_indices=cell_grid.Domain_Indices().Thickened(ghost_cells);dimension_start=domain_indices.Minimum_Corner();dimension_end=domain_indices.Maximum_Corner();
}
//#####################################################################
// Destructor
//#####################################################################
template<class T, int d> NONMANIFOLD_FAST_MARCHING_METHOD_UNIFORM<T,d>::
~NONMANIFOLD_FAST_MARCHING_METHOD_UNIFORM()
{}

//#####################################################################
// Function Fast_Marching_Method
//#####################################################################
template<class T, int d> void NONMANIFOLD_FAST_MARCHING_METHOD_UNIFORM<T,d>::
Fast_Marching_Method(HYBRID_ARRAY<T,d>& phi, HYBRID_ARRAY<bool,d>& done)
{
#if 1
    int heap_length=0;
    HYBRID_ARRAY<int,d> close_k(domain_range, flat_range); 
    ARRAY<HINDEX> heap(phi.Size(),false); // a generic version of heap
    // Need to check what this does?
    Initialize_Interface(phi,done,close_k,heap,heap_length);

    while(heap_length != 0){
        HINDEX hindex=heap(1); // smallest point is on top of heap
        if(stopping_distance && abs(phi(hindex)) > stopping_distance){ // exit early
            for(HYBRID_ARRAY_ITERATOR<d> iterator(domain_range, flat_range); iterator.Valid(); iterator.Next())
                if(!done(iterator.Index())){
                    phi(iterator.Index()) = LEVELSET_UTILITIES<T>::Sign(phi(iterator.Index()))*stopping_distance;}
            return;}
        done(hindex)=true;close_k(hindex)=0; // add to done, remove from close
        //LOG::cout << "\n========== MARK " << hindex << " AS DONE =========\n " << std::endl;

        FAST_MARCHING<T>::Down_Heap(phi,close_k,heap,heap_length);heap_length--; // remove point from heap

        const AA_NEIGHBORS& neighbors = cell_mesh.Node_Neighbors(hindex);            
        for(int neighbor_direction=1; neighbor_direction <= 2*d; neighbor_direction++){
            for( int i = 1; i <= neighbors(neighbor_direction).m; i++){
                if( !done(neighbors(neighbor_direction)(i)))
                    Update_Or_Add_Neighbor(phi, done, close_k, heap, heap_length,neighbors(neighbor_direction)(i));
            }                
        }
    } 
#else
    int heap_length=0;
    T_ARRAYS_INT close_k(cell_grid.Domain_Indices(ghost_cells+1)); // extra cell so that it is the same size as the done array for optimizations
    ARRAY<TV_INT> heap(cell_grid.Domain_Indices().Thickened(ghost_cells).Size(),false); // a generic version of heap((m+6)*(n+6)*(mn+6),false);
    Initialize_Interface(phi_ghost,done,close_k,heap,heap_length,add_seed_indices_for_ghost_cells);

    while(heap_length != 0){
        TV_INT index=heap(1); // smallest point is on top of heap
        if(stopping_distance && abs(phi_ghost(index)) > stopping_distance){ // exit early
            for(CELL_ITERATOR iterator(cell_grid);iterator.Valid();iterator.Next()) if(!done(iterator.Cell_Index())){
                phi_ghost(iterator.Cell_Index())=LEVELSET_UTILITIES<T>::Sign(phi_ghost(iterator.Cell_Index()))*stopping_distance;}
            return;}
        done(index)=true;close_k(index)=0; // add to done, remove from close
        FAST_MARCHING<T>::Down_Heap(phi_ghost,close_k,heap,heap_length);heap_length--; // remove point from heap

        if(nonmanifold_levelset.collision_body_list){
            PHYSBAM_NOT_IMPLEMENTED();
            for(int axis=1;axis<=T_GRID::dimension;axis++){TV_INT axis_vector=TV_INT::Axis_Vector(axis);
                if(index[axis] != dimension_start[axis] && !done(index-axis_vector) && Neighbor_Visible(axis,index-axis_vector))
                    Update_Or_Add_Neighbor(phi_ghost,done,close_k,heap,heap_length,index-axis_vector);
                if(index[axis] != dimension_end[axis] && !done(index+axis_vector) && Neighbor_Visible(axis,index))
                    Update_Or_Add_Neighbor(phi_ghost,done,close_k,heap,heap_length,index+axis_vector);}}
        else
            for(int axis=1;axis<=T_GRID::dimension;axis++){TV_INT axis_vector=TV_INT::Axis_Vector(axis);
                if(index[axis] != dimension_start[axis] && !done(index-axis_vector))
                    Update_Or_Add_Neighbor(phi_ghost,done,close_k,heap,heap_length,index-axis_vector);
                if(index[axis] != dimension_end[axis] && !done(index+axis_vector))
                    Update_Or_Add_Neighbor(phi_ghost,done,close_k,heap,heap_length,index+axis_vector);}}
#endif
}
//#####################################################################
// Function Update_Or_Add_Neighbor
//#####################################################################
template<class T, int d> inline void NONMANIFOLD_FAST_MARCHING_METHOD_UNIFORM<T,d>::
Update_Or_Add_Neighbor(HYBRID_ARRAY<T,d>& phi, HYBRID_ARRAY<bool,d>& done,
                                    HYBRID_ARRAY<int,d>& close_k,ARRAY<HINDEX>& heap, int& heap_length,
                                    const HINDEX& neighbor)
{
#if 1
    if(close_k(neighbor)){
        Update_Close_Point(phi,done,neighbor);
        FAST_MARCHING<T>::Up_Heap(phi,close_k,heap,close_k(neighbor));}
    else{close_k(neighbor)=1; // add to close 
        Update_Close_Point(phi,done,neighbor);
        heap_length++;heap(heap_length)=neighbor;
        FAST_MARCHING<T>::Up_Heap(phi,close_k,heap,heap_length);}
#else
    if(close_k(neighbor)){
        Update_Close_Point(phi_ghost,done,neighbor);
        FAST_MARCHING<T>::Up_Heap(phi_ghost,close_k,heap,close_k(neighbor));}
    else{close_k(neighbor)=1; // add to close 
        Update_Close_Point(phi_ghost,done,neighbor);
        heap_length++;heap(heap_length)=neighbor;
        FAST_MARCHING<T>::Up_Heap(phi_ghost,close_k,heap,heap_length);}
#endif
}
//#####################################################################
// Function Initialize_Interface
//#####################################################################
// pass heap_length by reference
template<class T, int d> void NONMANIFOLD_FAST_MARCHING_METHOD_UNIFORM<T,d>::
Initialize_Interface(HYBRID_ARRAY<T,d>& phi, HYBRID_ARRAY<bool,d>& done,
                     HYBRID_ARRAY<int,d>& close_k,ARRAY<HINDEX>& heap, int& heap_length)
{

#if 1
  
    for(HYBRID_ARRAY_ITERATOR<d> iterator(domain_range, flat_range); iterator.Valid(); iterator.Next() ){
        if(done(iterator.Index())){
            //LOG::cout << "Seed point " << iterator.Index() << " has value: " << phi(iterator.Index()) << std::endl;
            Add_To_Initial(done,close_k,iterator.Index());
        }
    }
    //LOG::cout << std::endl;

    // initialize close points
    for(HYBRID_ARRAY_ITERATOR<d> iterator(domain_range, flat_range); iterator.Valid(); iterator.Next() ){
        if(close_k(iterator.Index())){
            Update_Close_Point(phi,done,iterator.Index());
            heap_length++;
            heap(heap_length)=iterator.Index();
            FAST_MARCHING<T>::Up_Heap(phi,close_k,heap,heap_length);
        }
    }
    
    
#else
    T_NONMANIFOLD_LEVELSET nonmanifold_levelset_ghost(cell_grid,phi_ghost);

    for(CELL_ITERATOR iterator(cell_grid,ghost_cells);iterator.Valid();iterator.Next()) if(done(iterator.Cell_Index())) Add_To_Initial(done,close_k,iterator.Cell_Index());
    if(add_seed_indices_for_ghost_cells){RANGE<TV_INT> ghost_domain=cell_grid.Domain_Indices().Thickened(ghost_cells);
        for(CELL_ITERATOR iterator(cell_grid,ghost_cells,T_GRID::GHOST_REGION);iterator.Valid();iterator.Next()){TV_INT index=iterator.Cell_Index();
            for(int i=1;i<=T_GRID::number_of_neighbors_per_cell;i++){TV_INT neighbor_index(iterator.Cell_Neighbor(i));
                if(ghost_domain.Lazy_Inside(neighbor_index) && LEVELSET_UTILITIES<T>::Interface(phi_ghost(index),phi_ghost(neighbor_index))){
                    if(!done(index))Add_To_Initial(done,close_k,index);
                    if(!done(neighbor_index))Add_To_Initial(done,close_k,neighbor_index);}}}}

   // initialize close points
    for(CELL_ITERATOR iterator(cell_grid,ghost_cells);iterator.Valid();iterator.Next()) if(close_k(iterator.Cell_Index())){
        Update_Close_Point(phi_ghost,done,iterator.Cell_Index());
        heap_length++;heap(heap_length)=iterator.Cell_Index();
        FAST_MARCHING<T>::Up_Heap(phi_ghost,close_k,heap,heap_length);}
#endif
}
//#####################################################################
// Function Update_Close_Point
//##################################################################### 
// needs done=0 around the outside of the domain
template<class T, int d> void NONMANIFOLD_FAST_MARCHING_METHOD_UNIFORM<T,d>::
Update_Close_Point(HYBRID_ARRAY<T,d>& phi,HYBRID_ARRAY<bool,d>& done, const HINDEX& hindex)
{
#if 1  
    T value[d] = {}; // the phi value to use in the given direction
    int number_of_axis=0; // the number of axis that we want to use later
    T dx[d] = {};

    //LOG::cout << "Updating point: " << hindex << std::endl;

    const AA_NEIGHBORS& neighbors = cell_mesh.Node_Neighbors(hindex);            
    for(int neighbor_axis=0; neighbor_axis < d; neighbor_axis++){
        int high = ((neighbor_axis)*2)+2;
        int low = ((neighbor_axis)*2)+1;
        value[number_of_axis] = FLT_MAX;
        dx[number_of_axis] = cell_mesh.dx;
        bool available = false;
        
        for( int i = 1; i <= neighbors(high).m; i++){
            if(done(neighbors(high)(i))){
                //LOG::cout << "Neighbor " << neighbors(high)(i) << " on axis spoke " << high << " has value " << phi(neighbors(high)(i)) << std::endl;
                value[number_of_axis] = minmag( value[number_of_axis], phi(neighbors(high)(i)) );
                available = true;
            }
        }                

        for( int i = 1; i <= neighbors(low).m; i++){
            if(done(neighbors(low)(i))){
                //LOG::cout << "Neighbor " << neighbors(low)(i) << " on axis spoke " << low << " has value " << phi(neighbors(low)(i)) << std::endl;
                value[number_of_axis] = minmag( value[number_of_axis], phi(neighbors(low)(i)) );
                available = true;
            }
        }                

        if(available){
            number_of_axis++;
        }
    }
    
    //for(int i = 0; i < d; i++){
    //    LOG::cout << "axis: " << value[i] << "  " << (unsigned int)(Hash( value[i] )) << std::endl;
    //    LOG::cout << "dx: " << dx[i] << "  " << (unsigned int)(Hash( dx[i] )) << std::endl;
    //   LOG::cout << "number of axis: " << number_of_axis << std::endl;
    //}

    //LOG::cout << hindex.x << "  " << (unsigned int)(Hash( phi(hindex) )) << std::endl;

    if( number_of_axis ) // Only do this if we have a neighbor (NOT ALL NODES HAVE NEIGHBORS!!!)
        phi(hindex)=FAST_MARCHING<T>::template Solve_Close_Point<d>(phi(hindex),number_of_axis,value,dx);
    {
        //ARRAY<T,T_INDEX> phi_copy(domain_range);
        //for(HYBRID_ARRAY_ITERATOR<d> iterator(domain_range, flat_range); iterator.Valid(); iterator.Next())
        //phi_copy(iterator.Index().x) = phi(iterator.Index());
        
        //LOG::cout << hindex.x << "  " << (unsigned int)(Hash( phi(hindex) )) << std::endl;
    }

    //LOG::cout << "New value for point: " << phi(hindex) << std::endl << std::endl;

#else
    T value[T_GRID::dimension]={}; // the phi value to use in the given direction
    T dx[T_GRID::dimension]; // the edge length in the given direction
    int number_of_axis=0; // the number of axis that we want to use later

    // check each principal axis
    for(int axis=1;axis<=T_GRID::dimension;axis++){TV_INT axis_vector=TV_INT::Axis_Vector(axis),low=index-axis_vector,high=index+axis_vector;
        bool check_low=done(low),check_high=done(high);
        if(nonmanifold_levelset.collision_body_list){
            if(check_low && !Neighbor_Visible(axis,low)) check_low=false;
            if(check_high && !Neighbor_Visible(axis,index)) check_high=false;}
        dx[number_of_axis]=cell_grid.dX[axis];
        if(!check_low){
            if(check_high)value[number_of_axis]=phi_ghost(high);
            else number_of_axis--;}
        else if(!check_high) value[number_of_axis]=phi_ghost(low);
        else value[number_of_axis]=minmag(phi_ghost(low),phi_ghost(high));
        number_of_axis++;}

    phi_ghost(index)=FAST_MARCHING<T>::template Solve_Close_Point<T_GRID::dimension>(phi_ghost(index),number_of_axis,value,dx);
#endif
}
//#####################################################################
// Function Add_To_Initial
//##################################################################### 
template<class T, int d> void NONMANIFOLD_FAST_MARCHING_METHOD_UNIFORM<T,d>::
Add_To_Initial(HYBRID_ARRAY<bool,d>& done,HYBRID_ARRAY<int,d>& close_k,const HINDEX& hindex)
{
#if 1
    done(hindex)=true;close_k(hindex)=0; // add to done, remove from close 
/*
    if(nonmanifold_levelset.collision_body_list){
        for(int axis=1;axis<=T_GRID::dimension;axis++){TV_INT axis_vector=TV_INT::Axis_Vector(axis);
            if(!done(index-axis_vector) && Neighbor_Visible(axis,index-axis_vector)) close_k(index-axis_vector)=1;
            if(!done(index+axis_vector) && Neighbor_Visible(axis,index)) close_k(index+axis_vector)=1;}}
            else{*/

    // add neighbors to close if not done     
    const AA_NEIGHBORS& neighbors = cell_mesh.Node_Neighbors(hindex);            
    for(int neighbor_direction=1; neighbor_direction <= 2*d; neighbor_direction++){
        for( int i = 1; i <= neighbors(neighbor_direction).m; i++){
            if(!done(neighbors(neighbor_direction)(i)))
                close_k(neighbors(neighbor_direction)(i))=1;
        }                
    }
#else
    done(index)=true;close_k(index)=0; // add to done, remove from close 
    // add neighbors to close if not done 
    if(nonmanifold_levelset.collision_body_list){
        for(int axis=1;axis<=T_GRID::dimension;axis++){TV_INT axis_vector=TV_INT::Axis_Vector(axis);
            if(!done(index-axis_vector) && Neighbor_Visible(axis,index-axis_vector)) close_k(index-axis_vector)=1;
            if(!done(index+axis_vector) && Neighbor_Visible(axis,index)) close_k(index+axis_vector)=1;}}
    else{
        for(int axis=1;axis<=T_GRID::dimension;axis++){TV_INT axis_vector=TV_INT::Axis_Vector(axis);
            if(!done(index-axis_vector)) close_k(index-axis_vector)=1;
            if(!done(index+axis_vector)) close_k(index+axis_vector)=1;}}
#endif
}
//#####################################################################
//template class NONMANIFOLD_FAST_MARCHING_METHOD_UNIFORM<float, 1 >;
template class NONMANIFOLD_FAST_MARCHING_METHOD_UNIFORM<float, 2 >;
template class NONMANIFOLD_FAST_MARCHING_METHOD_UNIFORM<float, 3 >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
//template class NONMANIFOLD_FAST_MARCHING_METHOD_UNIFORM<double, 1 >;
template class NONMANIFOLD_FAST_MARCHING_METHOD_UNIFORM<double, 2 >;
template class NONMANIFOLD_FAST_MARCHING_METHOD_UNIFORM<double, 3 >;
#endif

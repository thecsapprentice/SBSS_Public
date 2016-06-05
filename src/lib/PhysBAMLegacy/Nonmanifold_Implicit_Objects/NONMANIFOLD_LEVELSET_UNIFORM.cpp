//#####################################################################
// Copyright 2014, Nathan Mitchell, Raj Setaluri.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class NONMANIFOLD_LEVELSET_UNIFORM
//#####################################################################
#include <Nonmanifold_Implicit_Objects/NONMANIFOLD_LEVELSET_UNIFORM.h>
#include <PhysBAM_Tools/Interpolation/LINEAR_INTERPOLATION.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T, int d> NONMANIFOLD_LEVELSET_UNIFORM<T,d>::
NONMANIFOLD_LEVELSET_UNIFORM(T_MESH& mesh_input,T_ARRAY_SCALAR& phi_input)
    :mesh(mesh_input),phi(phi_input),normals(0),dx(mesh.dx),one_over_dx((T)1./dx)
{}
//#####################################################################
// Destructor
//#####################################################################
template<class T, int d> NONMANIFOLD_LEVELSET_UNIFORM<T,d>::
~NONMANIFOLD_LEVELSET_UNIFORM()
{
    if(normals) delete normals;
}
//#####################################################################
// Function Lazy_Phi
//#####################################################################
template<class T, int d> T NONMANIFOLD_LEVELSET_UNIFORM<T,d>::
Lazy_Phi(const int cell_index,const TV& weights) const
{
    T phi_node[T_MESH::nodes_per_cell];
    for(int i=1;i<=T_MESH::nodes_per_cell;i++)
        phi_node[i-1]=phi(mesh.cells(cell_index).nodes(i));
    return LINEAR_INTERPOLATION<T,T>::Linear(phi_node,weights.Reversed()); // weights reversed due to LINEAR_INTERPOLATION interface
}
//#####################################################################
// Function Distance_To_Closest_Point
//#####################################################################
template<class T, int d> T NONMANIFOLD_LEVELSET_UNIFORM<T,d>::
Distance_To_Closest_Point(const int input_cell_index,const TV& input_weights,int& output_cell_index,TV& output_weights,TV& normal) const
{
    const T small_threshold=(T).001;
    const T zero_threshold=(T)1e-6;
    const int max_iterations=4;
    const T phi_threshold=(T)1e-4;
    const T distance_threshold=(T)0.0;

    const TV X_original=mesh.Location(input_cell_index,input_weights);
    const T phi_original=Lazy_Phi(input_cell_index,input_weights);
    int current_cell_index=input_cell_index;
    TV entry_weights=input_weights,exit_weights;

    normal=Lazy_Normal(input_cell_index,input_weights);
    if(phi_original<0) normal*=(T)-1.;
    TV one_over_normal;
    for(int v=1;v<=d;v++)
        if(fabs(normal(v))<zero_threshold) one_over_normal(v)=std::numeric_limits<T>::max();
        else if(fabs(normal(v))<small_threshold) one_over_normal(v)=(normal(v)>0)?std::numeric_limits<T>::max():-std::numeric_limits<T>::max();
        else one_over_normal(v)=(T)1./normal(v);
    T phi_entry,phi_exit;
    while(true){
        phi_entry=Lazy_Phi(current_cell_index,entry_weights);
        T time_of_flight,min_positive_time_of_flight=std::numeric_limits<T>::max();
        int axis=0,side=0;
        for(int v=1;v<=d;v++) for(int s=1;s<=2;s++){
            if(fabs(one_over_normal(v))<std::numeric_limits<T>::max())
                time_of_flight=one_over_normal(v)*((s==1)?(-entry_weights(v)):((T)1.-entry_weights(v)));
            else time_of_flight=one_over_normal(v);
            if(time_of_flight>0&&time_of_flight<min_positive_time_of_flight){
                min_positive_time_of_flight=time_of_flight;axis=v;side=s;}}
        PHYSBAM_ASSERT(axis!=0 && side!=0); // debug
        exit_weights=RANGE<TV>::Unit_Box().Clamp(entry_weights+min_positive_time_of_flight*normal);
        phi_exit=Lazy_Phi(current_cell_index,exit_weights);
        const bool same_sign=!((phi_entry>(T)0.) ^ (phi_exit>(T)0.));
        if(same_sign){ // same sign -- go to next cell
            if(((X_original-mesh.Location(current_cell_index,exit_weights)).Magnitude()-fabs(phi_original)) > distance_threshold){ // traveled too far
                //LOG::cout << "Should not happen often. "<< std::endl;
                const T entry_distance=(X_original-mesh.Location(current_cell_index,entry_weights)).Magnitude();
                PHYSBAM_ASSERT((entry_distance-fabs(phi_original))<distance_threshold); // debug
                //PHYSBAM_ASSERT((entry_distance<phi_original)); // debug
                output_weights=entry_weights+/*one_over_dx**/((T)fabs(phi_original)-entry_distance)*normal;
                PHYSBAM_ASSERT(output_weights.Min() >=0.0 && output_weights.Max() <= 1.0 );
                output_cell_index=current_cell_index;
                break;}
            const int flat_face_index=(axis-1)*2+side;
            const TV_INT face_neighbor=(side==1)?-TV_INT::Axis_Vector(axis):TV_INT::Axis_Vector(axis);
            entry_weights=RANGE<TV>::Unit_Box().Clamp(exit_weights-TV(face_neighbor));
            ARRAY<int> cell_neighbors=mesh.cells(current_cell_index).cell_neighbors(flat_face_index);
            const int transition_index=mesh.cells(current_cell_index).transition(flat_face_index).x;
            if(transition_index){ // has transitions...append transition neighbors as well
                const typename T_MESH::TRANSITION_FACE& t=mesh.transition_faces(transition_index);
                const ARRAY<PAIR<int,int> >* transition_neighbors=(side==1)?(&(t.left_cells)):(&(t.right_cells));
                for(int i=1;i<=transition_neighbors->m;i++) cell_neighbors.Append((*transition_neighbors)(i).x);}
            if(cell_neighbors.m==0) PHYSBAM_FATAL_ERROR("No face neighbors found in backtrace!!"); // no face neighbors!!
            T min_phi=std::numeric_limits<T>::max();
            int next_cell_index=0;
            for(int i=1;i<=cell_neighbors.m;i++){
                const int neighbor_cell_index=cell_neighbors(i);
                const T neighbor_phi=Lazy_Phi(neighbor_cell_index,entry_weights);
                if(neighbor_phi<min_phi){
                    min_phi=neighbor_phi;
                    next_cell_index=neighbor_cell_index;}}
            if(transition_index && (min_phi > 0)){ // transition and no negative neighbor -- quit here
                output_cell_index=current_cell_index;
                output_weights=exit_weights;
                break;}
            else {
                //LOG::cout << "Transitioning to another cell: " << current_cell_index << " --> " << next_cell_index << "  ("<< transition_index << ")"<< std::endl;
                current_cell_index=next_cell_index;  // otherwise go to next cell
            }
        }
        else{                       // different sign -- search along line within cell
            //LOG::cout << "Searching along line." << std::endl;
            T phi_left=phi_entry, phi_right=phi_exit;
            TV weights_left=entry_weights, weights_right=exit_weights;
            TV weights_guess=weights_left,d_weights;
            int iterations=0;
            while(iterations<max_iterations){iterations++;    
                d_weights=weights_right-weights_left;
                weights_guess=weights_left+(-phi_left/(phi_right-phi_left))*d_weights;
                const T phi_guess=Lazy_Phi(current_cell_index,weights_guess);
                //LOG::cout << "I " << iterations << ": "<< weights_guess << ",  " << phi_guess << ", " << phi_left << ", " << phi_right << std::endl;
                if(phi_guess<phi_threshold) break;
                const bool guesses_same_sign=!((phi_left>(T)0.0)^(phi_guess>(T)0.0));
                if(!guesses_same_sign){
                    phi_right=phi_guess;
                    weights_right=weights_guess;
                }
                else{
                    phi_left=phi_guess;
                    weights_left=weights_guess;
                }
            }
            output_cell_index=current_cell_index;
            output_weights=weights_guess;
            break;}}
    return (mesh.Location(output_cell_index,output_weights)-X_original).Magnitude();
}
//#####################################################################
#define INSTANTIATION_HELPER(T,d) \
    template class NONMANIFOLD_LEVELSET_UNIFORM< T , d >;

INSTANTIATION_HELPER(float,2);
INSTANTIATION_HELPER(float,3);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
INSTANTIATION_HELPER(double,2);
INSTANTIATION_HELPER(double,3);
#endif

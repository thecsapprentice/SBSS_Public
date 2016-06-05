#include "SKINNING_NONLINEAR_ELASTICITY.h"
#include <Common/RANGE_ITERATOR.h>

using namespace PhysBAM;

//#####################################################################
// Function Initialize_Blocks_Constraints
//#####################################################################
template<class T,int d,bool enable_constraints,bool enable_muscles> void SKINNING_NONLINEAR_ELASTICITY<T,d,enable_constraints,enable_muscles>::
Initialize_Blocks_Constraints(const int number_of_partitions)
{
    LOG::SCOPE scope("SKINNING_NONLINEAR_ELASTICITY::Initialize_Blocks_Constraints()");

#if 0
    ARRAY<int> constraint_base_offsets;
    ARRAY<int> constraint_base_lengths;
    ARRAY<BLOCKED_TYPE<TV,8> > constraint_multilinear_coordinates;
    ARRAY<BLOCKED_TYPE<T,8> > constraint_spring_constants;

    int number_of_active_blocks=BASE::cell_block_base_indices.m;
    int block_constraint_offset = 1;

    constraint_base_offsets.Resize(number_of_active_blocks);
    constraint_base_lengths.Resize(number_of_active_blocks);

    ARRAY<ARRAY<PAIR<TV,T> >,T_INDEX> cell_constraints;
    cell_constraints.Resize(padded_cell_domain,true,false);


    for(int c=1;c<=static_point_constraints.m;c++)
        switch(static_point_constraints(c).type){
        case POINT_CONSTRAINT<T,d>::EMBEDDED_POINT_TO_FIXED_POINT_SPRING_TYPE:
            {
                const EMBEDDED_POINT_TO_FIXED_POINT_SPRING_CONSTRAINT<T,d>& embedded_point_to_fixed_point_spring=static_point_constraints(c).embedded_point_to_fixed_point_spring;
                PAIR<TV,T> spring_constraint(embedded_point_to_fixed_point_spring.Multilinear_Coordinates(),embedded_point_to_fixed_point_spring.Spring_Coefficient());
                cell_constraints(embedded_point_to_fixed_point_spring.Cell_Index()).Append(spring_constraint);
            }
            break;
          case POINT_CONSTRAINT<T,d>::TWO_EMBEDDED_POINT_SPRING_TYPE:
              {
                  PHYSBAM_FATAL_ERROR("Cannot add Two point Embedded constraints to the static constraint list.");
              }
              break;
        default:
                break;
        }




    for( int block = 1; block <= number_of_active_blocks; block++)
        {
           const T_INDEX& base_index = BASE::cell_block_base_indices(block);
   
            //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            //
            //              Setup Blocked Constraint Data
            //
            //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

            //LOG::cout << "Computing Constraint Data for block " << block << std::endl;

            int block_constraint_count = 0;
            
            for(RANGE_ITERATOR<d> cell_iterator(RANGE<T_INDEX>(base_index,base_index+1));
                cell_iterator.Valid();
                cell_iterator.Next()){
                const T_INDEX& index = cell_iterator.Index();
                //LOG::cout << "Cell " << index << "has a constraint count of " << cell_constraints(index).m << std::endl;
                block_constraint_count = max(cell_constraints(index).Size(), block_constraint_count);
            }

            
            constraint_base_offsets(block) = block_constraint_offset;
            constraint_base_lengths(block) = block_constraint_count;
            block_constraint_offset += block_constraint_count;
            
            //LOG::cout << "Constraint Block " << block << " has max level of " << block_constraint_count << std::endl;
            for( int block_constraint = 0; block_constraint < block_constraint_count; block_constraint++)
                {
                    //LOG::cout << "Computing block level of " << block_constraint << std::endl;

                    constraint_multilinear_coordinates.Append(BLOCKED_TYPE<TV,8>());
                    constraint_spring_constants.Append(BLOCKED_TYPE<T,8>());

                    int cell=0;
                    for(RANGE_ITERATOR<d> cell_iterator(RANGE<T_INDEX>(base_index,base_index+1));
                        cell_iterator.Valid();
                        cell_iterator.Next(),cell++)
                        {
                            const T_INDEX& index = cell_iterator.Index();
                            if(cell_constraints(index).m > block_constraint)
                                {
                                    constraint_multilinear_coordinates.Last().Set(cell_constraints(index)(block_constraint+1).x, cell);
                                    constraint_spring_constants.Last().Set(cell_constraints(index)(block_constraint+1).y, cell);
                                 }
                            else
                                {
                                    constraint_multilinear_coordinates.Last().Set(TV(.5,.5,.5), cell);
                                    constraint_spring_constants.Last().Set(T(0), cell);
                                }

                            if( !(cell_type(index)==INTERIOR_CELL_TYPE ||
                                  (/*first_order && */cell_type(index)==BOUNDARY_CELL_TYPE)) )
                                {
                                    constraint_multilinear_coordinates.Last().Set(TV(.5,.5,.5), cell);
                                    constraint_spring_constants.Last().Set(T(0), cell);
                                }

                        }          
                }
            

            BASE::specialized_data.InitializeConstraintData(constraint_base_offsets, constraint_base_lengths, constraint_multilinear_coordinates, constraint_spring_constants);

        }
#endif
}

//#####################################################################
// Function Initialize_Blocks_Muscles
//#####################################################################
template<class T,int d,bool enable_constraints,bool enable_muscles> void SKINNING_NONLINEAR_ELASTICITY<T,d,enable_constraints,enable_muscles>::
Initialize_Blocks_Muscles(const int number_of_partitions)
{
    LOG::SCOPE scope("SKINNING_NONLINEAR_ELASTICITY::Initialize_Blocks_Muscles()");
#ifdef USE_SPECIALIZED_KERNELS


    const int DEFAULT_MUSCLE = 1;

    int number_of_active_blocks=BASE::cell_block_base_indices.m;
    int block_muscle_offset = 1;

    ARRAY<int> muscle_base_offsets;
    ARRAY<int> muscle_base_lengths;
    ARRAY<BLOCKED_TYPE<int,8> > muscle_id;
    ARRAY<BLOCKED_TYPE<TV,8> > muscle_fiber;
    ARRAY<BLOCKED_TYPE<TV,8> > muscle_F_fiber;
    ARRAY<BLOCKED_TYPE<T,8> > muscle_density;
    ARRAY<BLOCKED_TYPE<T,8> > muscle_c1;
    ARRAY<BLOCKED_TYPE<T,8> > muscle_c2;

    muscle_base_offsets.Resize(number_of_active_blocks);
    muscle_base_lengths.Resize(number_of_active_blocks);

    for( int block = 1; block <= number_of_active_blocks; block++)
        {
           const T_INDEX& base_index = BASE::cell_block_base_indices(block);
            

            //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            //
            //              Setup Blocked Muscle Data
            //
            //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


            //LOG::cout << "Computing Muscle Data for block " << block << std::endl;

            int block_muscle_count = 0;

            for(RANGE_ITERATOR<d> cell_iterator(RANGE<T_INDEX>(base_index,base_index+1));
                cell_iterator.Valid();
                cell_iterator.Next()){
                const T_INDEX& index = cell_iterator.Index();
                //LOG::cout << "Cell " << index << "has a muscle count of " << this->cell_muscles(index).m << std::endl;
                block_muscle_count = max(this->cell_muscles(index).Size(), block_muscle_count);
            }

            //LOG::cout << "Muscle Block " << block << " has max level of " << block_muscle_count << std::endl;

            muscle_base_offsets(block) = block_muscle_offset;
            muscle_base_lengths(block) = block_muscle_count;
            block_muscle_offset += block_muscle_count;

            for( int block_muscle = 0; block_muscle < block_muscle_count; block_muscle++)
                {
                    //LOG::cout << "Computing block level of " << block_muscle << std::endl;

                    muscle_id.Append(BLOCKED_TYPE<int,8>());
                    muscle_fiber.Append(BLOCKED_TYPE<TV,8>());
                    muscle_density.Append(BLOCKED_TYPE<T,8>());

                    muscle_F_fiber.Append(BLOCKED_TYPE<TV,8>());
                    muscle_c1.Append(BLOCKED_TYPE<T,8>());
                    muscle_c2.Append(BLOCKED_TYPE<T,8>());

                    int cell=0;
                    for(RANGE_ITERATOR<d> cell_iterator(RANGE<T_INDEX>(base_index,base_index+1));
                        cell_iterator.Valid();
                        cell_iterator.Next(),cell++)
                        {
                            const T_INDEX& index = cell_iterator.Index();
                            if(this->cell_muscles(index).m > block_muscle)
                                {
                                    muscle_id.Last().Set(this->cell_muscles(index)(block_muscle+1)-1, cell);
                                    muscle_fiber.Last().Set(this->cell_fibers(index)(block_muscle+1), cell); 
                                    muscle_density.Last().Set(this->cell_densities(index)(block_muscle+1), cell); 
                                 }
                            else
                                {
                                    muscle_id.Last().Set(DEFAULT_MUSCLE-1,cell);
                                    muscle_fiber.Last().Set(TV(1,0,0), cell); 
                                    muscle_density.Last().Set(T(0), cell); 
                                }

                            if( !(this->cell_type(index)==INTERIOR_CELL_TYPE ||
                                  (/*first_order && */this->cell_type(index)==BOUNDARY_CELL_TYPE)) )
                                {
                                    muscle_id.Last().Set(DEFAULT_MUSCLE-1,cell);
                                    muscle_fiber.Last().Set(TV(1,0,0), cell); 
                                    muscle_density.Last().Set(T(0), cell); 
                                }

                        }          
                }
        }

    BASE::specialized_data.InitializeMuscleData(muscle_base_offsets, muscle_base_lengths, muscle_id,
                                                muscle_fiber, muscle_density);

#endif
}

template class SKINNING_NONLINEAR_ELASTICITY<float,3,true,true>;
template class SKINNING_NONLINEAR_ELASTICITY<float,3,true,false>;
#ifndef USE_SPECIALIZED_KERNELS
template class SKINNING_NONLINEAR_ELASTICITY<double,3,true,true>;
template class SKINNING_NONLINEAR_ELASTICITY<double,3,true,false>;
#endif

#ifndef __SPECIALIZED_KERNELS_DATA_H__
#define __SPECIALIZED_KERNELS_DATA_H__

#include <Common/ALIGNED_ARRAY.h>

namespace PhysBAM
{
    #define BUNDLE_WIDTH 16 
    

    template<class T,int width> class BLOCKED_TYPE
    {
        enum {size=sizeof(T)/sizeof(float)};
        typedef typename T::SCALAR T_SCALAR;
    public:
        T_SCALAR data[size][width];
        void Get(T& obj,const int i) const
        {
            typedef T_SCALAR (&refArray)[size];
            refArray obj_ref=reinterpret_cast<refArray>(obj);
            for(int j=0;j<size;j++) obj_ref[j]=data[j][i];
        }
        void Set(const T& obj,const int i)
        {
            typedef const T_SCALAR (&refArray)[size];
            refArray obj_ref=reinterpret_cast<refArray>(obj);
            for(int j=0;j<size;j++) data[j][i]=obj_ref[j];
        }
        template<int upwidth>
            void CopyUp( BLOCKED_TYPE<T,upwidth>& new_block, int position) const{
            PHYSBAM_ASSERT(upwidth%width==0);
            PHYSBAM_ASSERT(position*width<=upwidth);
            for(int i=0;i<width;i++)
                for(int j=0;j<size;j++)
                    new_block.data[j][i+(position*width)] = data[j][i];
        }

    };
    
    template<int width> class BLOCKED_TYPE<float,width>
        {
        public:
            float data[width];
            void Get(float& obj,const int i) const
            {
                obj=data[i];
            }
            void Set(const float& obj,const int i)
            {
                data[i]=obj;
            }
            template<int upwidth>
                void CopyUp( BLOCKED_TYPE<float,upwidth>& new_block, int position) const{
                PHYSBAM_ASSERT(upwidth%width==0);
                PHYSBAM_ASSERT(position*width<=upwidth);
                for(int i=0;i<width;i++)
                    new_block.data[i+(position*width)] = data[i];
            }
        };
    
    template<int width> class BLOCKED_TYPE<int,width>
        {
        public:
            int data[width];
            void Get(int& obj,const int i) const
            {
                obj=data[i];
            }
            void Set(const int& obj,const int i)
            {
                data[i]=obj;
            }
            template<int upwidth>
                void CopyUp( BLOCKED_TYPE<int,upwidth>& new_block, int position) const{
                PHYSBAM_ASSERT(upwidth%width==0);
                PHYSBAM_ASSERT(position*width<=upwidth);
                for(int i=0;i<width;i++)
                    new_block.data[i+(position*width)] = data[i];
            }
        };
    
    template<class T, int d>
    class SPECIALIZED_KERNEL_DATA
    {
    private:
        bool materials_initialized;
        bool muscles_initialized;
        bool constraints_initialized;

    public:
        static const int VECTOR_WIDTH=BUNDLE_WIDTH;
        static const int ALIGNMENT=VECTOR_WIDTH * 4;
        STATIC_ASSERT(VECTOR_WIDTH%8==0);
        static const int VECTOR_MULT=BUNDLE_WIDTH/8;

        typedef VECTOR<T,d> TV;
        typedef VECTOR<int,d> T_INDEX;

        typedef BLOCKED_TYPE<T,VECTOR_WIDTH> BUNDLED_T;
        typedef BLOCKED_TYPE<TV,VECTOR_WIDTH> BUNDLED_TV;
        typedef BLOCKED_TYPE<int,VECTOR_WIDTH> BUNDLED_INT;
        typedef BLOCKED_TYPE<MATRIX<T,d>,VECTOR_WIDTH> BUNDLED_MATRIX;
        typedef BLOCKED_TYPE<DIAGONAL_MATRIX<T,d>,VECTOR_WIDTH> BUNDLED_DIAG_MATRIX;
        typedef BLOCKED_TYPE<ROTATED_STRESS_DERIVATIVE<T,d>,VECTOR_WIDTH> BUNDLED_RSD;

        typedef BLOCKED_TYPE<T,8> BLOCKED_T;
        typedef BLOCKED_TYPE<TV,8> BLOCKED_TV;
        typedef BLOCKED_TYPE<int,8> BLOCKED_INT;
        typedef BLOCKED_TYPE<MATRIX<T,d>,8> BLOCKED_MATRIX;
        typedef BLOCKED_TYPE<DIAGONAL_MATRIX<T,d>,8> BLOCKED_DIAG_MATRIX;
        typedef BLOCKED_TYPE<ROTATED_STRESS_DERIVATIVE<T,d>,8> BLOCKED_RSD;

        // Flattened material parameters

        ALIGNED_ARRAY<BUNDLED_T, int, ALIGNMENT > mu_bundled;
        ALIGNED_ARRAY<BUNDLED_T, int, ALIGNMENT > mu_stab_bundled;
        ALIGNED_ARRAY<BUNDLED_T, int, ALIGNMENT > kappa_bundled;
        ALIGNED_ARRAY<BUNDLED_T, int, ALIGNMENT > alpha_bundled;
        ALIGNED_ARRAY<BUNDLED_T, int, ALIGNMENT > alpha_sqr_over_kappa_bundled;
        ALIGNED_ARRAY<BUNDLED_T, int, ALIGNMENT > cutoff_bundled;
        ALIGNED_ARRAY<BUNDLED_T, int, ALIGNMENT > h_bundled;
        ALIGNED_ARRAY<BUNDLED_T, int, ALIGNMENT > one_over_h_bundled;
        ALIGNED_ARRAY<BUNDLED_T, int, ALIGNMENT > cell_volume_bundled;
        
        // Flattened elasticity data

        ALIGNED_ARRAY<BUNDLED_MATRIX, int, ALIGNMENT> U_bundled;
        ALIGNED_ARRAY<BUNDLED_MATRIX, int, ALIGNMENT > V_bundled;
        ALIGNED_ARRAY<BUNDLED_DIAG_MATRIX, int, ALIGNMENT > Sigma_bundled;
        ALIGNED_ARRAY<BUNDLED_DIAG_MATRIX, int, ALIGNMENT > Q_hat_bundled;
        ALIGNED_ARRAY<BUNDLED_DIAG_MATRIX, int, ALIGNMENT > P_hat_bundled;
        ALIGNED_ARRAY<BUNDLED_RSD, int, ALIGNMENT > dPdF_bundled;
        mutable ARRAY<VECTOR<T,81> > u_compact_array;
        mutable ARRAY<VECTOR<T,8> > p_compact_array;
        mutable ARRAY<VECTOR<T,81> > d_compact_array;

        // Flattened muscle data
        
        ALIGNED_ARRAY<int, int, ALIGNMENT> muscle_base_offsets;
        ALIGNED_ARRAY<int, int, ALIGNMENT> muscle_base_lengths;
        ALIGNED_ARRAY<BUNDLED_INT, int, ALIGNMENT > muscle_id;
        ALIGNED_ARRAY<BUNDLED_TV, int, ALIGNMENT > muscle_fiber;
        ALIGNED_ARRAY<BUNDLED_TV, int, ALIGNMENT > muscle_F_fiber;
        ALIGNED_ARRAY<BUNDLED_T, int, ALIGNMENT > muscle_density;
        ALIGNED_ARRAY<BUNDLED_T, int, ALIGNMENT > muscle_c1;
        ALIGNED_ARRAY<BUNDLED_T, int, ALIGNMENT > muscle_c2;
        
        // Flattened constraint data
        
        ALIGNED_ARRAY<int, int, ALIGNMENT> constraint_base_offsets;
        ALIGNED_ARRAY<int, int, ALIGNMENT> constraint_base_lengths;
        ALIGNED_ARRAY<BUNDLED_TV, int, ALIGNMENT > constraint_multilinear_coordinates;
        ALIGNED_ARRAY<BUNDLED_T, int, ALIGNMENT > constraint_spring_constants;       
        ALIGNED_ARRAY<BUNDLED_TV, int, ALIGNMENT > constraint_node_positions;
        ALIGNED_ARRAY<BUNDLED_INT, int, ALIGNMENT > spring_id;
        ALIGNED_ARRAY<BUNDLED_INT, int, ALIGNMENT > spring_id_X;
        ALIGNED_ARRAY<BUNDLED_INT, int, ALIGNMENT > spring_id_Y;
        ALIGNED_ARRAY<BUNDLED_INT, int, ALIGNMENT > spring_id_Z;
        ARRAY<T,int>* collision_spring_constants;
        ARRAY<TV,int>* collision_spring_locations;

        SPECIALIZED_KERNEL_DATA() : materials_initialized(false), muscles_initialized(false), constraints_initialized(false), collision_spring_constants(NULL),  collision_spring_locations(NULL) {};            

        void InitializeMaterialParameters(const ARRAY<BLOCKED_T>& mu_blocked, const ARRAY<BLOCKED_T>& mu_stab_blocked, const ARRAY<BLOCKED_T>& kappa_blocked, const ARRAY<BLOCKED_T>& alpha_blocked, const ARRAY<BLOCKED_T>& alpha_sqr_over_kappa_blocked, const T cutoff_blocked, const T h_blocked, const T one_over_h_blocked, const T cell_volume_blocked){
            LOG::SCOPE scope("SPECIALIZED_KERNEL_DATA::InitializeMaterialParameters");

            int bundle_multiplier = (VECTOR_WIDTH/8);
            int blocks = mu_blocked.m;
            int bundles = blocks%bundle_multiplier ? (blocks/bundle_multiplier)+1 : (blocks/bundle_multiplier);
            int filler_blocks = blocks%bundle_multiplier ? bundle_multiplier - (blocks%bundle_multiplier) : 0;
            
            LOG::cout << "Creating " << bundles << " bundles of blocks. " << std::endl;
            if(filler_blocks) LOG::cout << "We are missing blocks.  Need to add " << filler_blocks
                                        << " filler blocks to complete bundles" << std::endl;

            mu_bundled.Resize(bundles);
            mu_stab_bundled.Resize(bundles);
            kappa_bundled.Resize(bundles);
            alpha_bundled.Resize(bundles);
            alpha_sqr_over_kappa_bundled.Resize(bundles);

            U_bundled.Resize(bundles);
            V_bundled.Resize(bundles);
            Sigma_bundled.Resize(bundles);
            P_hat_bundled.Resize(bundles);
            Q_hat_bundled.Resize(bundles);
            dPdF_bundled.Resize(bundles);

            u_compact_array.Exact_Resize(bundles*bundle_multiplier);
            p_compact_array.Exact_Resize(bundles*bundle_multiplier);
            d_compact_array.Exact_Resize(bundles*bundle_multiplier);

            for( int filler_block=blocks+1; filler_block <= blocks+filler_blocks; filler_block++){
                int bundle = ((filler_block-1)/bundle_multiplier)+1;
                int bundle_position = ((filler_block-1)%bundle_multiplier);
                for(int i=0;i<VECTOR_WIDTH;i++){
                    mu_bundled(bundle).Set(T(),i);
                    mu_stab_bundled(bundle).Set(T(),i);
                    kappa_bundled(bundle).Set(T(),i);
                    alpha_bundled(bundle).Set(T(),i);
                    alpha_sqr_over_kappa_bundled(bundle).Set(T(),i);
                }
            }

            for( int block=1; block <= blocks; block++){
                int bundle = ((block-1)/bundle_multiplier)+1;
                int bundle_position = ((block-1)%bundle_multiplier);
                mu_blocked(block).CopyUp(mu_bundled(bundle), bundle_position);
                mu_stab_blocked(block).CopyUp(mu_stab_bundled(bundle), bundle_position);
                kappa_blocked(block).CopyUp(kappa_bundled(bundle), bundle_position);
                alpha_blocked(block).CopyUp(alpha_bundled(bundle), bundle_position);
                alpha_sqr_over_kappa_blocked(block).CopyUp(alpha_sqr_over_kappa_bundled(bundle), bundle_position);
            }


            cutoff_bundled.Resize( 1 );
            one_over_h_bundled.Resize( 1 );
            cell_volume_bundled.Resize( 1 );
            h_bundled.Resize( 1 );

            for(int cell=0;cell<VECTOR_WIDTH;cell++){
                cutoff_bundled(1).Set(cutoff_blocked,cell);
                h_bundled(1).Set(h_blocked,cell);
                one_over_h_bundled(1).Set(one_over_h_blocked,cell);
                cell_volume_bundled(1).Set(cell_volume_blocked,cell);}
               
            materials_initialized = true;
            CheckDataConsistency();
        }

        void InitializeMuscleData(const ARRAY<int>& muscle_base_offsets_blocked, const ARRAY<int>& muscle_base_lengths_blocked, const ARRAY<BLOCKED_INT>& muscle_id_blocked,
                                  const ARRAY<BLOCKED_TV>& muscle_fiber_blocked,   const ARRAY<BLOCKED_T>& muscle_density_blocked ){
            LOG::SCOPE scope("SPECIALIZED_KERNEL_DATA::InitializeMuscleData");

            const int bundle_multiplier = (VECTOR_WIDTH/8);
            int blocks = muscle_base_offsets_blocked.m;
            int bundles = blocks%bundle_multiplier ? (blocks/bundle_multiplier)+1 : (blocks/bundle_multiplier);
            int filler_blocks = blocks%bundle_multiplier ? bundle_multiplier - (blocks%bundle_multiplier) : 0;

            LOG::cout << "Creating " << bundles << " bundles of blocks. " << std::endl;
            if(filler_blocks) LOG::cout << "We are missing blocks.  Need to add " << filler_blocks
                                        << " filler blocks to complete bundles" << std::endl;

            muscle_base_offsets.Resize(bundles);
            muscle_base_lengths.Resize(bundles);

            // Perform bundling step
            int bundle_muscle_offset = 1;
            muscle_base_offsets.Fill(-1);
            muscle_base_lengths.Fill(-1);
            
            for( int block = 1; block <= blocks; block+=bundle_multiplier ){
                const int bundle = ((block-1)/bundle_multiplier)+1;
                const int bundle_position = ((block-1) % bundle_multiplier);

                int prep_blocks[bundle_multiplier];
                for(int i=0;i<bundle_multiplier;i++){
                    prep_blocks[i] = block+i;
                    if(prep_blocks[i] > blocks)
                        prep_blocks[i]=-1;
                }
                
                int bundle_muscle_count = muscle_base_lengths_blocked(prep_blocks[0]);
                for(int i=1;i<bundle_multiplier;i++)
                    if(prep_blocks[i] > 0)
                        bundle_muscle_count=max(bundle_muscle_count, muscle_base_lengths_blocked(prep_blocks[i]));
                


                muscle_base_offsets(bundle) = bundle_muscle_offset;
                muscle_base_lengths(bundle) = bundle_muscle_count;
                bundle_muscle_offset += bundle_muscle_count;
                
                
                int prep_block_offsets[bundle_multiplier];
                for( int bundle_muscle = 0; bundle_muscle < bundle_muscle_count; bundle_muscle++){
                    for(int i=0;i<bundle_multiplier;i++)
                        if(prep_blocks[i] > 0 && muscle_base_lengths_blocked(prep_blocks[i]) > bundle_muscle)
                            prep_block_offsets[i] = muscle_base_offsets_blocked(prep_blocks[i]) + bundle_muscle;
                        else
                            prep_block_offsets[i] = 0;

                    
                    muscle_id.Append(BUNDLED_INT());
                    muscle_fiber.Append(BUNDLED_TV());
                    muscle_F_fiber.Append(BUNDLED_TV());
                    muscle_density.Append(BUNDLED_T());
                    muscle_c1.Append(BUNDLED_T());
                    muscle_c2.Append(BUNDLED_T());

                    for(int i=0;i<VECTOR_WIDTH;i++){
                        muscle_id.Last().Set(0,i);
                        muscle_fiber.Last().Set(TV(),i);
                        muscle_density.Last().Set(T(0),i);
                    }
                    
                    for(int i=0; i<bundle_multiplier; i++)
                        if(prep_block_offsets[i] > 0){
                            muscle_id_blocked(prep_block_offsets[i]).CopyUp(muscle_id.Last(), i);
                            muscle_fiber_blocked(prep_block_offsets[i]).CopyUp(muscle_fiber.Last(), i);
                            muscle_density_blocked(prep_block_offsets[i]).CopyUp(muscle_density.Last(), i);
                        }
                }
            }             

            muscles_initialized = true;
            CheckDataConsistency();
        }

        void InitializeConstraintData(const ARRAY<int>& constraint_base_offsets_blocked, const ARRAY<int>& constraint_base_lengths_blocked, const ARRAY<BLOCKED_TV>& constraint_multilinear_coordinates_blocked, const ARRAY<BLOCKED_T>& constraint_spring_constants_blocked, const ARRAY<BLOCKED_TV>& constraint_node_positions_blocked,  const ARRAY<BLOCKED_INT>& spring_id_blocked, const ARRAY<BLOCKED_INT>& spring_id_X_blocked, const ARRAY<BLOCKED_INT>& spring_id_Y_blocked, const ARRAY<BLOCKED_INT>& spring_id_Z_blocked, ARRAY<T, int>* collision_spring_constants_input, ARRAY<TV,int>* collision_spring_locations_input){
            LOG::SCOPE scope("SPECIALIZED_KERNEL_DATA::InitializeConstraintData");
            const int bundle_multiplier = (VECTOR_WIDTH/8);
            int blocks = constraint_base_offsets_blocked.m;
            int bundles = blocks%bundle_multiplier ? (blocks/bundle_multiplier)+1 : (blocks/bundle_multiplier);
            int filler_blocks = blocks%bundle_multiplier ? bundle_multiplier - (blocks%bundle_multiplier) : 0;

            LOG::cout << "Creating " << bundles << " bundles of blocks. " << std::endl;
            if(filler_blocks) LOG::cout << "We are missing blocks.  Need to add " << filler_blocks
                                        << " filler blocks to complete bundles" << std::endl;

            constraint_base_offsets.Resize(bundles);
            constraint_base_lengths.Resize(bundles);

            // Perform bundling step
            int bundle_constraint_offset = 1;
            constraint_base_offsets.Fill(-1);
            constraint_base_lengths.Fill(-1);
            
            for( int block = 1; block <= blocks; block+=bundle_multiplier ){
                const int bundle = ((block-1)/bundle_multiplier)+1;
                const int bundle_position = ((block-1) % bundle_multiplier);

                int prep_blocks[bundle_multiplier];
                for(int i=0;i<bundle_multiplier;i++){
                    prep_blocks[i] = block+i;
                    if(prep_blocks[i] > blocks)
                        prep_blocks[i]=-1;
                }
                
                int bundle_constraint_count = constraint_base_lengths_blocked(prep_blocks[0]);
                for(int i=1;i<bundle_multiplier;i++)
                    if(prep_blocks[i] > 0)
                        bundle_constraint_count=max(bundle_constraint_count, constraint_base_lengths_blocked(prep_blocks[i]));
                


                constraint_base_offsets(bundle) = bundle_constraint_offset;
                constraint_base_lengths(bundle) = bundle_constraint_count;
                bundle_constraint_offset += bundle_constraint_count;
                
                
                int prep_block_offsets[bundle_multiplier];
                for( int bundle_constraint = 0; bundle_constraint < bundle_constraint_count; bundle_constraint++){
                    for(int i=0;i<bundle_multiplier;i++)
                        if(prep_blocks[i] > 0 && constraint_base_lengths_blocked(prep_blocks[i]) > bundle_constraint )
                            prep_block_offsets[i] = constraint_base_offsets_blocked(prep_blocks[i]) + bundle_constraint;
                        else
                            prep_block_offsets[i] = 0;

                    constraint_multilinear_coordinates.Append(BUNDLED_TV());
                    constraint_spring_constants.Append(BUNDLED_T());
                    constraint_node_positions.Append(BUNDLED_TV());
                    spring_id.Append(BUNDLED_INT());
                    spring_id_X.Append(BUNDLED_INT());
                    spring_id_Y.Append(BUNDLED_INT());
                    spring_id_Z.Append(BUNDLED_INT());

                    for(int i=0;i<VECTOR_WIDTH;i++){
                        constraint_multilinear_coordinates.Last().Set(TV::All_Ones_Vector()*0.5, i);
                        constraint_spring_constants.Last().Set(T(0), i);
                        constraint_node_positions.Last().Set(TV(), i);
                        spring_id.Last().Set(0, i);
                        spring_id_X.Last().Set(0, i);
                        spring_id_Y.Last().Set(1, i);
                        spring_id_Z.Last().Set(2, i);                        
                    }
                    
                    for(int i=0; i<bundle_multiplier; i++)
                        if(prep_block_offsets[i] > 0){
                            constraint_multilinear_coordinates_blocked(prep_block_offsets[i]).CopyUp(constraint_multilinear_coordinates.Last(), i);
                            constraint_spring_constants_blocked(prep_block_offsets[i]).CopyUp(constraint_spring_constants.Last(), i);
                            constraint_node_positions_blocked(prep_block_offsets[i]).CopyUp(constraint_node_positions.Last(), i);
                            spring_id_blocked(prep_block_offsets[i]).CopyUp(spring_id.Last(), i);
                            spring_id_X_blocked(prep_block_offsets[i]).CopyUp(spring_id_X.Last(), i);
                            spring_id_Y_blocked(prep_block_offsets[i]).CopyUp(spring_id_Y.Last(), i);
                            spring_id_Z_blocked(prep_block_offsets[i]).CopyUp(spring_id_Z.Last(), i);                            
                        }
                }
            }             

            collision_spring_constants = collision_spring_constants_input;
            collision_spring_locations = collision_spring_locations_input;

            constraints_initialized = true;
            CheckDataConsistency();
        }

        void CheckDataConsistency(){
            int bundles = -1;
            
            if(materials_initialized){
                bundles = mu_bundled.m;
                PHYSBAM_ASSERT(mu_stab_bundled.m == bundles);
                PHYSBAM_ASSERT(kappa_bundled.m == bundles);
                PHYSBAM_ASSERT(alpha_bundled.m == bundles);
                PHYSBAM_ASSERT(alpha_sqr_over_kappa_bundled.m == bundles);}

            if(muscles_initialized)
                if(bundles!=-1){
                    PHYSBAM_ASSERT(muscle_base_offsets.m == bundles);
                    PHYSBAM_ASSERT(muscle_base_lengths.m == bundles);}
                else{
                    bundles = muscle_base_offsets.m;
                    PHYSBAM_ASSERT(muscle_base_lengths.m == bundles);}

            if(constraints_initialized){
                if(bundles!=-1){
                    PHYSBAM_ASSERT(constraint_base_offsets.m == bundles);
                    PHYSBAM_ASSERT(constraint_base_lengths.m == bundles);}
                else{
                    bundles = constraint_base_offsets.m;
                    PHYSBAM_ASSERT(constraint_base_lengths.m == bundles);}
                collision_spring_constants != NULL;
                collision_spring_locations != NULL;
            }
                       
        }

        bool MaterialsInitialized() const { return materials_initialized; };
        bool MusclesInitialized() const { return muscles_initialized; };
        bool ConstraintsInitialized() const { return constraints_initialized; };
    };


}

#endif

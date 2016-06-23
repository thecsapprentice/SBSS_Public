#include <PhysBAM_Tools/Data_Structures/HASHTABLE_ITERATOR.h>
#include <PhysBAM_Tools/Krylov_Solvers/SYMMQMR.h>
#include <PhysBAM_Tools/Matrices/SYMMETRIC_MATRIX_2X2.h>
#include <PhysBAM_Tools/Matrices/SYMMETRIC_MATRIX_3X3.h>
#include <PhysBAM_Geometry/Geometry_Particles/GEOMETRY_PARTICLES.h>
#include <PhysBAM_Tools/Random_Numbers/RANDOM_NUMBERS.h>
#include <PhysBAM_Geometry/Basic_Geometry/BOX.h>
#include "NONLINEAR_ELASTICITY.h"
#include "CG_SYSTEM.h"
#include "CG_VECTOR.h"
#include <Common_Tools/Math_Tools/RANGE_ITERATOR.h>
#include "ROTATED_STRESS_DERIVATIVE.h"

//For the Haxks
#include <iostream>
#include <fstream>

struct NEOHOOKEAN_TAG;
struct COROTATED_TAG;
#if defined(BIPHASIC_SUPPORT)
struct BIPHASIC_TAG;
#endif
struct __m256;
struct __m128;
#ifdef ENABLE_MIC
struct __m512;
#endif

#include <Thread_Queueing/PTHREAD_QUEUE.h>
using namespace PhysBAM;
extern PTHREAD_QUEUE* pthread_queue;

#include <SIMD_Optimized_Kernels/Kernels/Grid_Compact/Grid_Compact.h>
#include <SIMD_Optimized_Kernels/Kernels/Grid_UnCompact/Grid_UnCompact.h>
#include <SIMD_Optimized_Kernels/Kernels/Block_Duplicate/Block_Duplicate.h>
#include <SIMD_Optimized_Kernels/Kernels/Block_DeDuplicate/Block_DeDuplicate.h>

#include <SIMD_Optimized_Kernels/Kernels/Update_Position_Based_State_Blocked/Update_Position_Based_State_Blocked.h>
#include <SIMD_Optimized_Kernels/Kernels/Update_Position_Based_State/Update_Position_Based_State.h>
#include <SIMD_Optimized_Kernels/Kernels/Muscle_Update_Position_Based_State_Blocked/Muscle_Update_Position_Based_State_Blocked.h>
#include <SIMD_Optimized_Kernels/Kernels/Muscle_Update_Position_Based_State/Muscle_Update_Position_Based_State.h>

#include <SIMD_Optimized_Kernels/Kernels/Add_Force_First_Order/Add_Force_First_Order.h>
#include <SIMD_Optimized_Kernels/Kernels/Muscle_Forces/Muscle_Forces.h>
#include <SIMD_Optimized_Kernels/Kernels/Collision_Forces/Collision_Forces.h>
#include <SIMD_Optimized_Kernels/Kernels/Collision_Force_Differential/Collision_Force_Differential.h>
#include <SIMD_Optimized_Kernels/Kernels/Add_Force_Differential/Add_Force_Differential.h>
#include <SIMD_Optimized_Kernels/Kernels/Force_Stabilization/Force_Stabilization.h>
#include <SIMD_Optimized_Kernels/Kernels/Muscle_Differential_Forces/Muscle_Differential_Forces.h>
#include <SIMD_Optimized_Kernels/Kernels/Constraint_Differential_Forces/Constraint_Differential_Forces.h>
#include <SIMD_Optimized_Kernels/Kernels/Unweighted_Gradient/Unweighted_Gradient.h>
#include <SIMD_Optimized_Kernels/Kernels/Unweighted_Accumulation/Unweighted_Accumulation.h>
#include <SIMD_Optimized_Kernels/Kernels/Muscle_Differential/Muscle_Differential.h>

#include "SPECIALIZED_KERNELS.h"

#ifdef USE_SPECIALIZED_KERNELS

#if defined(BIPHASIC_SUPPORT)
namespace nm_biphasic {
    float* biphasic_threshold=NULL ;
    float* biphasic_factor=NULL;
}
#endif

namespace {
    template <class T_MATERIAL> struct T_MATERIAL_DESC { };
    
    template<> struct T_MATERIAL_DESC< PhysBAM::NEOHOOKEAN<float, 3> >
    {
        typedef NEOHOOKEAN_TAG T_MATERIAL;
    };

    template<> struct T_MATERIAL_DESC< PhysBAM::COROTATED<float, 3> >
    {
        typedef COROTATED_TAG T_MATERIAL;
    };
#if defined(BIPHASIC_SUPPORT)
    template<> struct T_MATERIAL_DESC< PhysBAM::BIPHASIC<float, 3> >
    {
        typedef BIPHASIC_TAG T_MATERIAL;
    };
#endif

    template <class T_ARCH> struct T_ARCH_DESC { };

    template<> struct T_ARCH_DESC<float>
    {
        static const int T_DATA_SIZE = 16;
        static const int WIDTH = 1;
        typedef float T_DATA[T_DATA_SIZE];
        typedef int I_DATA[T_DATA_SIZE];
        typedef T_MATERIAL_DESC< NONLINEAR_ELASTICITY<float,3>::T_CONSTITUTIVE_MODEL >::T_MATERIAL T_MATERIAL;
    };
    
    template<> struct T_ARCH_DESC<__m128>
    {
        static const int T_DATA_SIZE = 16;
        static const int WIDTH = 4;
        typedef float T_DATA[T_DATA_SIZE];
        typedef int I_DATA[T_DATA_SIZE];
        typedef T_MATERIAL_DESC< NONLINEAR_ELASTICITY<float,3>::T_CONSTITUTIVE_MODEL >::T_MATERIAL T_MATERIAL;
    };

    
    template<> struct T_ARCH_DESC<__m256>
    {
        static const int T_DATA_SIZE = 16;
        static const int WIDTH = 8;
        typedef float T_DATA[T_DATA_SIZE];
        typedef int I_DATA[T_DATA_SIZE];
        typedef T_MATERIAL_DESC< NONLINEAR_ELASTICITY<float,3>::T_CONSTITUTIVE_MODEL >::T_MATERIAL T_MATERIAL;
    };

#ifdef ENABLE_MIC
    template<> struct T_ARCH_DESC<__m512>
    {
        static const int T_DATA_SIZE = 16;
        static const int WIDTH = 16;
        typedef float T_DATA[T_DATA_SIZE];
        typedef int I_DATA[T_DATA_SIZE];
        typedef T_MATERIAL_DESC< NONLINEAR_ELASTICITY<float,3>::T_CONSTITUTIVE_MODEL >::T_MATERIAL T_MATERIAL;
    };
#endif

    typedef float (&fArr_W)[16];
    typedef const float (&cfArr_W)[16];

    typedef int (&iArr_W)[16];
    typedef const int (&ciArr_W)[16];

    typedef float (&fArr_3_8_W)[3][8][16];
    typedef const float (&cfArr_3_8_W)[3][8][16];

    typedef float (&fArr_3_W)[3][16];
    typedef const float (&cfArr_3_W)[3][16];

    typedef float (&fArr_9_W)[9][16];
    typedef const float (&cfArr_9_W)[9][16];

    typedef float (&fArr_12_W)[12][16];
    typedef const float (&cfArr_12_W)[12][16];

#define KCAST(T_ARCH_TYPE, VAR) reinterpret_cast<T_ARCH_TYPE>(VAR)
        

}

#ifdef ENABLE_MIC
#define SELECTED_ARCH __m512
#else
//#define SELECTED_ARCH __m256
#define SELECTED_ARCH float
#endif
//********************************************************************************************
//
//                   Update Position Based State: Threaded Task
//
//********************************************************************************************
    template<class T, bool enable_constraints, bool enable_muscles>
        class Update_Position_Based_State_Task : public PhysBAM::PTHREAD_QUEUE::TASK
        {
        private:
            typedef VECTOR<T,3> TV;

            int _bundle_start;
            int _bundle_count;

            T (*_u_compact)[3][27];
            T (*_p_compact)[8];
            T (*_d_compact)[3][27];
            const ARRAY<T>& _muscle_fiber_max_stresses;
            const ARRAY<T>& _muscle_activations;
            SPECIALIZED_KERNEL_DATA<T,3>& _specialized_data;


        public:

//#####################################################################
// Constructor
//#####################################################################
            Update_Position_Based_State_Task(
              int bundle_start,
              int bundle_count,                              
              T (*u_compact)[3][27],
              T (*p_compact)[8],   
              T (*d_compact)[3][27],
              const ARRAY<T>& muscle_fiber_max_stresses,
              const ARRAY<T>& muscle_activations,
              SPECIALIZED_KERNEL_DATA<T,3>& specialized_data
                                             )
                : 
                _bundle_start(bundle_start),
                _bundle_count(bundle_count),
                _u_compact(u_compact),                             
                _p_compact(p_compact),  
                _d_compact(d_compact),                           
                _muscle_fiber_max_stresses(muscle_fiber_max_stresses),
                _muscle_activations(muscle_activations),
                _specialized_data(specialized_data)
            {
#if defined(BIPHASIC_SUPPORT)
                if( nm_biphasic::biphasic_threshold == NULL )
                    nm_biphasic::biphasic_threshold = new float[SPECIALIZED_KERNEL_DATA<T,3>::VECTOR_WIDTH];
                if( nm_biphasic::biphasic_factor == NULL )
                    nm_biphasic::biphasic_factor = new float[SPECIALIZED_KERNEL_DATA<T,3>::VECTOR_WIDTH];
#endif
            }


//#####################################################################
// Function: Run
//#####################################################################
void Run()
{
    typedef SELECTED_ARCH T_ARCH;

    STATIC_ASSERT(T_ARCH_DESC<T_ARCH>::T_DATA_SIZE == SPECIALIZED_KERNEL_DATA<T,3>::VECTOR_WIDTH);

#if defined(BIPHASIC_SUPPORT)
    std::ifstream biphasic_params( "biphasic_parameters" );
    float bpt, bpf;
    biphasic_params >> bpt;
    biphasic_params >> bpf;
    for( int i = 0; i <= SPECIALIZED_KERNEL_DATA<T,3>::VECTOR_WIDTH; i++){
        nm_biphasic::biphasic_threshold[i] = bpt;
        nm_biphasic::biphasic_factor[i] = bpf;
    }
    biphasic_params.close();
#endif

    for( int _bundle = 0; _bundle < _bundle_count; _bundle++)
        {
            int bundle = _bundle + _bundle_start + 1;
                       
            typedef T pressure_p[T_ARCH_DESC<T_ARCH>::T_DATA_SIZE];
            
//            __declspec(align(SPECIALIZED_KERNEL_DATA<T,3>::ALIGNMENT)) T due[3][8][T_ARCH_DESC<T_ARCH>::T_DATA_SIZE];
//            __declspec(align(SPECIALIZED_KERNEL_DATA<T,3>::ALIGNMENT)) pressure_p dpe; 
            __attribute__ ((aligned( BUNDLE_WIDTH * 4))) T due[3][8][T_ARCH_DESC<T_ARCH>::T_DATA_SIZE];
            __attribute__ ((aligned( BUNDLE_WIDTH * 4))) T dde[3][8][T_ARCH_DESC<T_ARCH>::T_DATA_SIZE];
            __attribute__ ((aligned( BUNDLE_WIDTH * 4))) pressure_p dpe;
            memset(dde,0,sizeof(T)*3*8*T_ARCH_DESC<T_ARCH>::T_DATA_SIZE);
                      
            int first_block = (bundle-1)*SPECIALIZED_KERNEL_DATA<T,3>::VECTOR_MULT;
            for( int i = 0; i < T_ARCH_DESC<T_ARCH>::T_DATA_SIZE; i += 8, first_block++){
                Block_Duplicate<T_ARCH>(_u_compact[first_block], KCAST(fArr_3_8_W, due[0][0][i]));
                for(int j=0; j<8; j++)
                    KCAST(fArr_W, dpe[i])[j] = _p_compact[first_block][j];
            }
            
            {
                //LOG::SCOPE scope("SPECIALIZED_KERNELS::UPBS");

            for( int i = 0; i < T_ARCH_DESC<T_ARCH>::T_DATA_SIZE; i += T_ARCH_DESC<T_ARCH>::WIDTH)
                ::Update_Position_Based_State<
                    T_ARCH_DESC<T_ARCH>::T_MATERIAL,
                    T_ARCH,
                    typename  T_ARCH_DESC<T_ARCH>::T_DATA,
                    typename T_ARCH_DESC<T_ARCH>::I_DATA>::Run(
                    KCAST(cfArr_3_8_W, due[0][0][i]),
                    KCAST(cfArr_W, dpe[i]),
                    KCAST(cfArr_W, (_specialized_data.mu_bundled)(bundle).data[i]),
                    KCAST(cfArr_W, (_specialized_data.mu_stab_bundled)(bundle).data[i]),
                    KCAST(cfArr_W, (_specialized_data.kappa_bundled)(bundle).data[i]),
                    KCAST(cfArr_W, (_specialized_data.alpha_bundled)(bundle).data[i]),
                    KCAST(cfArr_W, (_specialized_data.cutoff_bundled)(1).data[i]),
                    KCAST(cfArr_W, (_specialized_data.one_over_h_bundled)(1).data[i]),
                    KCAST(cfArr_W, (_specialized_data.cell_volume_bundled)(1).data[i]),
                    KCAST(fArr_9_W, (_specialized_data.U_bundled)(bundle).data[0][i]),
                    KCAST(fArr_9_W, (_specialized_data.V_bundled)(bundle).data[0][i]),
                    KCAST(fArr_3_W, (_specialized_data.Sigma_bundled)(bundle).data[0][i]),
                    KCAST(fArr_3_W, (_specialized_data.Q_hat_bundled)(bundle).data[0][i]),
                    KCAST(fArr_12_W, (_specialized_data.dPdF_bundled)(bundle).data[0][i]),
                    KCAST(fArr_3_8_W, dde[0][0][i]));
            }

                        

            if( enable_muscles ){
                //LOG::SCOPE scope("SPECIALIZED_KERNELS::M-UPBS");

                int muscle_bundle_offset = (_specialized_data.muscle_base_offsets)(bundle);
                for( int layer=0; layer < (_specialized_data.muscle_base_lengths)(bundle); layer++)
                    {
                        int muscle_bundle = muscle_bundle_offset + layer;
                        for( int i = 0; i < T_ARCH_DESC<T_ARCH>::T_DATA_SIZE; i += T_ARCH_DESC<T_ARCH>::WIDTH)
                        ::Muscle_Update_Position_Based_State<
                            T_ARCH,
                            typename  T_ARCH_DESC<T_ARCH>::T_DATA,
                            typename  T_ARCH_DESC<T_ARCH>::I_DATA>(
                            KCAST(cfArr_3_8_W, due[0][0][i]),
                            KCAST(ciArr_W, (_specialized_data.muscle_id)(muscle_bundle).data[i]),
                            KCAST(cfArr_3_W, (_specialized_data.muscle_fiber)(muscle_bundle).data[0][i]),
                            KCAST(cfArr_W, (_specialized_data.muscle_density)(muscle_bundle).data[i]),
                            KCAST(cfArr_W, (_specialized_data.one_over_h_bundled)(1).data[i]),
                            KCAST(fArr_W, (_specialized_data.muscle_c1)(muscle_bundle).data[i]),
                            KCAST(fArr_W, (_specialized_data.muscle_c2)(muscle_bundle).data[i]),
                            KCAST(fArr_3_W, (_specialized_data.muscle_F_fiber)(muscle_bundle).data[0][i]),
                            (_muscle_activations).Get_Array_Pointer(),
                            (_muscle_fiber_max_stresses).Get_Array_Pointer()
                                                                   );
                    }
            }


            first_block = (bundle-1)*SPECIALIZED_KERNEL_DATA<T,3>::VECTOR_MULT;
            //for( int i = 0; i < T_ARCH_DESC<T_ARCH>::T_DATA_SIZE; i += 8, first_block++)
            //    Block_DeDuplicate<T_ARCH>(KCAST(fArr_3_8_W, dde[0][0][i]), _d_compact[first_block] );
                      
        }
}

    };
    
//********************************************************************************************
//
//                   Update Position Based State: Master Function
//
//********************************************************************************************

template <typename T, bool enable_constraints, bool enable_muscles>  void PhysBAM::
Update_Position_Based_State_With_Specialized_Kernels(
    T (*u_compact)[3][27],
    T (*p_compact)[8], 
    T (*d_compact)[3][27],
    const int number_of_blocks,
    const int number_of_partitions,
    const ARRAY<T>& muscle_fiber_max_stresses,
    const ARRAY<T>& muscle_activations,
    SPECIALIZED_KERNEL_DATA<T,3>& specialized_data)
{
    //LOG::SCOPE scope("SPECIALIZED_KERNELS::Update_Position_Based_State");

#ifdef USE_THREADED_KERNELS
    int bundles = (int)(ceil(number_of_blocks / (T)SPECIALIZED_KERNEL_DATA<T,3>::VECTOR_MULT));
    Update_Position_Based_State_Task<T,enable_constraints,enable_muscles>** tasks = new Update_Position_Based_State_Task<T,enable_constraints,enable_muscles>*[number_of_partitions];

    for(int partition=0;
        partition<number_of_partitions;
        partition++){
        int bundle_begin=(partition)*(bundles/number_of_partitions)+min(partition,bundles%number_of_partitions);
        int bundle_end=(partition+1)*(bundles/number_of_partitions)+min((partition+1),bundles%number_of_partitions);
        tasks[partition] =
            new Update_Position_Based_State_Task<T,enable_constraints,enable_muscles>(                                              
                                               bundle_begin,
                                               bundle_end-bundle_begin,
                                               u_compact,
                                               p_compact,
                                               d_compact,
                                               muscle_fiber_max_stresses,
                                               muscle_activations,
                                               specialized_data
                                                                                                                                    );
	}

#pragma omp parallel for num_threads(number_of_partitions)
    for(int partition=0;partition<number_of_partitions;partition++){
        tasks[partition]->Run();
        delete tasks[partition];
    }
    delete[] tasks;

 

#else
    int bundles = (int)(ceil(number_of_blocks / (T)SPECIALIZED_KERNEL_DATA<T,3>::VECTOR_MULT));
    Update_Position_Based_State_Task<T,enable_constraints,enable_muscles>* task =
        new Update_Position_Based_State_Task<T,enable_constraints,enable_muscles>(                                              
                                                                                  0,
                                                                                  bundles,
                                                                                  u_compact,
                                                                                  p_compact,
                                                                                  d_compact,
                                                                                  muscle_fiber_max_stresses,
                                                                                  muscle_activations,
                                                                                  specialized_data);
    task->Run();
    delete task;
#endif
}


template void PhysBAM::
Update_Position_Based_State_With_Specialized_Kernels<float,true,true>(
    float (*u_compact)[3][27],float (*p_compact)[8],float (*d_compact)[3][27],
    const int,const int,
    const ARRAY<float>&,const ARRAY<float>&,
    SPECIALIZED_KERNEL_DATA<float,3>&);


template void PhysBAM::
Update_Position_Based_State_With_Specialized_Kernels<float,true,false>(
    float (*u_compact)[3][27],float (*p_compact)[8],float (*d_compact)[3][27],
    const int,const int,
    const ARRAY<float>&,const ARRAY<float>&,
    SPECIALIZED_KERNEL_DATA<float,3>&);


template void PhysBAM::
Update_Position_Based_State_With_Specialized_Kernels<float,false,true>(
    float (*u_compact)[3][27],float (*p_compact)[8],float (*d_compact)[3][27],
    const int,const int,
    const ARRAY<float>&,const ARRAY<float>&,
    SPECIALIZED_KERNEL_DATA<float,3>&);


template void PhysBAM::
Update_Position_Based_State_With_Specialized_Kernels<float,false,false>(
    float (*u_compact)[3][27],float (*p_compact)[8],float (*d_compact)[3][27],
    const int,const int,
    const ARRAY<float>&,const ARRAY<float>&,
    SPECIALIZED_KERNEL_DATA<float,3>&);



//********************************************************************************************
//
//                   Add Force Differential: Threaded Task
//
//********************************************************************************************
    template<class T, bool enable_constraints, bool enable_muscles>
        class Set_Force_Differential_Task : public PhysBAM::PTHREAD_QUEUE::TASK
        {
        private:
            typedef VECTOR<T,3> TV;

            int _bundle_start;
            int _bundle_count;

            const ARRAY<int>& _explicit_bundle_list;
            const bool _use_bundle_list;

            T (*_du_compact)[3][27];
            T (*_dp_compact)[8];   
            T (*_df_compact)[3][27];
            T (*_dq_compact)[8];
            const SPECIALIZED_KERNEL_DATA<T,3>& _specialized_data;            

        public:

            long _time;
//#####################################################################
// Constructor
//#####################################################################

            Set_Force_Differential_Task(
              int bundle_start,
              int bundle_count,                              
              T (*du_compact)[3][27],
              T (*dp_compact)[8],   
              T (*df_compact)[3][27],
              T (*dq_compact)[8],
              const SPECIALIZED_KERNEL_DATA<T,3>& specialized_data
                                        )
                : 
                _bundle_start(bundle_start),
                _bundle_count(bundle_count),
                _du_compact(du_compact),                             
                _dp_compact(dp_compact),                             
                _df_compact(df_compact),                             
                _dq_compact(dq_compact),                             
                _specialized_data(specialized_data),
                _explicit_bundle_list(ARRAY<int>()),
                _use_bundle_list(false),
                _time(0)
                
            {  
            }

//#####################################################################
// Constructor
//#####################################################################

            Set_Force_Differential_Task(
              int bundle_start,
              int bundle_count, 
              const ARRAY<int>& bundle_list,
              T (*du_compact)[3][27],
              T (*dp_compact)[8],   
              T (*df_compact)[3][27],
              T (*dq_compact)[8],
              const SPECIALIZED_KERNEL_DATA<T,3>& specialized_data
                                        )
                : 
                _bundle_start(bundle_start),
                _bundle_count(bundle_count),
                _du_compact(du_compact),                             
                _dp_compact(dp_compact),                             
                _df_compact(df_compact),                             
                _dq_compact(dq_compact),                             
                _specialized_data(specialized_data),
                _explicit_bundle_list(bundle_list),
                _use_bundle_list(true),
                _time(0)
                
            {  
            }

//#####################################################################
// Function: Run
//#####################################################################

void Run()
    {
        typedef SELECTED_ARCH T_ARCH;

        _time = 0;

    STATIC_ASSERT(T_ARCH_DESC<T_ARCH>::T_DATA_SIZE == SPECIALIZED_KERNEL_DATA<T,3>::VECTOR_WIDTH);

    for( int _bundle = 0; _bundle < _bundle_count; _bundle++)
        {
            int bundle = _bundle + _bundle_start + 1;
                       
            if(_use_bundle_list)
                bundle = _explicit_bundle_list(bundle);

            typedef T pressure_p[T_ARCH_DESC<T_ARCH>::T_DATA_SIZE];
            
            //__declspec(align(SPECIALIZED_KERNEL_DATA<T,3>::ALIGNMENT)) T due[3][8][T_ARCH_DESC<T_ARCH>::T_DATA_SIZE];
            //__declspec(align(SPECIALIZED_KERNEL_DATA<T,3>::ALIGNMENT)) T dfe[3][8][T_ARCH_DESC<T_ARCH>::T_DATA_SIZE];
            __attribute__ ((aligned( BUNDLE_WIDTH * 4))) T due[3][8][T_ARCH_DESC<T_ARCH>::T_DATA_SIZE];
            __attribute__ ((aligned( BUNDLE_WIDTH * 4))) T dfe[3][8][T_ARCH_DESC<T_ARCH>::T_DATA_SIZE];
            //__declspec(align(SPECIALIZED_KERNEL_DATA<T,3>::ALIGNMENT)) pressure_p dpe;
            //__declspec(align(SPECIALIZED_KERNEL_DATA<T,3>::ALIGNMENT)) pressure_p dqe;
            __attribute__ ((aligned( BUNDLE_WIDTH * 4))) pressure_p dpe;
            __attribute__ ((aligned( BUNDLE_WIDTH * 4))) pressure_p dqe;

            //          ~0.0002
            {
                memset(dfe,0,sizeof(T)*3*8*T_ARCH_DESC<T_ARCH>::T_DATA_SIZE);
                memset(dqe,0,sizeof(T)*T_ARCH_DESC<T_ARCH>::T_DATA_SIZE);
            }
         

            //          ~0.0020 s
            {
                int first_block = (bundle-1)*SPECIALIZED_KERNEL_DATA<T,3>::VECTOR_MULT;
#ifdef ENABLE_MIC
                Block_Duplicate<T_ARCH>(reinterpret_cast<T (&) [2][3][27]>(_du_compact[first_block]), KCAST(fArr_3_8_W, due[0][0][0]));
                for( int i = 0; i < T_ARCH_DESC<T_ARCH>::T_DATA_SIZE; i++){
                    dpe[i] = reinterpret_cast<T (&) [16]>(_dp_compact[first_block])[i];
                }
#else
                for( int i = 0; i < T_ARCH_DESC<T_ARCH>::T_DATA_SIZE; i += 8, first_block++){
                    Block_Duplicate<T_ARCH>(_du_compact[first_block], KCAST(fArr_3_8_W, due[0][0][i]));
                    for(int j=0; j<8; j++)
                        KCAST(fArr_W, dpe[i])[j] = _dp_compact[first_block][j];
                }
#endif
            }

                
            //long start = __rdtsc();

            //          ~0.0025 s
            {
                for( int i = 0; i < T_ARCH_DESC<T_ARCH>::T_DATA_SIZE; i += T_ARCH_DESC<T_ARCH>::WIDTH)
                    ::Add_Force_Differential<
                        T_ARCH,
                        typename  T_ARCH_DESC<T_ARCH>::T_DATA,
                        typename T_ARCH_DESC<T_ARCH>::I_DATA>(
                                                              KCAST( cfArr_3_8_W, due[0][0][i]), 
                                                              KCAST( cfArr_W, dpe[i]),
                                                              KCAST( cfArr_W, (_specialized_data.alpha_sqr_over_kappa_bundled)(bundle).data[i]),
                                                              KCAST( cfArr_W, (_specialized_data.alpha_bundled)(bundle).data[i]),
                                                              KCAST( cfArr_W, (_specialized_data.one_over_h_bundled)(1).data[i]),
                                                              KCAST( cfArr_W, (_specialized_data.cell_volume_bundled)(1).data[i]),
                                                              KCAST( cfArr_3_W, (_specialized_data.Q_hat_bundled)(bundle).data[0][i]),
                                                              KCAST( cfArr_9_W, (_specialized_data.U_bundled)(bundle).data[0][i]),
                                                              KCAST( cfArr_9_W, (_specialized_data.V_bundled)(bundle).data[0][i]),
                                                              KCAST( cfArr_12_W, (_specialized_data.dPdF_bundled)(bundle).data[0][i]),
                                                              KCAST( fArr_3_8_W, dfe[0][0][i]),
                                                              KCAST( fArr_W, dqe[i]));
            }
            //long end = __rdtsc();
            //_time += (end - start);  
            
            {
                for( int i = 0; i < T_ARCH_DESC<T_ARCH>::T_DATA_SIZE; i += T_ARCH_DESC<T_ARCH>::WIDTH)
                    ::Force_Stabilization<
                        T_ARCH,
                        typename T_ARCH_DESC<T_ARCH>::T_DATA,
                        typename T_ARCH_DESC<T_ARCH>::I_DATA>(
                                                              KCAST( cfArr_3_8_W, due[0][0][i]), 
                                                              KCAST( cfArr_W, (_specialized_data.mu_stab_bundled)(bundle).data[i]),
                                                              KCAST( fArr_3_8_W, dfe[0][0][i]));
            }



            {

                if(enable_muscles){
                    int muscle_bundle_offset = (_specialized_data.muscle_base_offsets)(bundle);
                    
                    for( int layer=0; layer < (_specialized_data.muscle_base_lengths)(bundle); layer++)
                        {
                            int muscle_bundle = muscle_bundle_offset + layer;

                            for( int i = 0; i < T_ARCH_DESC<T_ARCH>::T_DATA_SIZE; i += T_ARCH_DESC<T_ARCH>::WIDTH)
                                ::Muscle_Differential_Forces<
                                    T_ARCH,
                                    typename T_ARCH_DESC<T_ARCH>::T_DATA,
                                    typename T_ARCH_DESC<T_ARCH>::I_DATA>(
                                    KCAST( fArr_3_8_W, dfe[0][0][i]),
                                    KCAST( cfArr_3_8_W, due[0][0][i]),
                                    KCAST( cfArr_3_W, (_specialized_data.muscle_fiber)(muscle_bundle).data[0][i]),
                                    KCAST( cfArr_3_W, (_specialized_data.muscle_F_fiber)(muscle_bundle).data[0][i]),
                                    KCAST( cfArr_W, (_specialized_data.muscle_c1)(muscle_bundle).data[i]),
                                    KCAST( cfArr_W, (_specialized_data.muscle_c2)(muscle_bundle).data[i]),
                                    KCAST( cfArr_W, (_specialized_data.one_over_h_bundled)(1).data[i]),
                                    KCAST( cfArr_W, (_specialized_data.cell_volume_bundled)(1).data[i]));
                        }
                } 
                
                if(enable_constraints){
                    int constraint_bundle_offset = (_specialized_data.constraint_base_offsets)(bundle);
                    
                    for( int layer=0; layer < (_specialized_data.constraint_base_lengths)(bundle); layer++)
                        {
                            int constraint_bundle = constraint_bundle_offset + layer;
                            for( int i = 0; i < T_ARCH_DESC<T_ARCH>::T_DATA_SIZE; i += T_ARCH_DESC<T_ARCH>::WIDTH)
                                ::Collision_Force_Differential<
                                    T_ARCH,
                                    typename T_ARCH_DESC<T_ARCH>::T_DATA,
                                    typename T_ARCH_DESC<T_ARCH>::I_DATA>(
                                    KCAST( fArr_3_8_W, dfe[0][0][i]),
                                    KCAST( cfArr_3_8_W, due[0][0][i]),
                                    KCAST( cfArr_3_W, (_specialized_data.constraint_multilinear_coordinates)(constraint_bundle).data[0][i]),
                                    KCAST( ciArr_W, (_specialized_data.spring_id)(constraint_bundle).data[i]),
                                    (_specialized_data.collision_spring_constants)->Get_Array_Pointer());
                        }
                }


            }

            //          ~0.0020 s
            {
                int first_block = (bundle-1)*SPECIALIZED_KERNEL_DATA<T,3>::VECTOR_MULT;
#ifdef ENABLE_MIC
                Block_DeDuplicate<T_ARCH>(KCAST(fArr_3_8_W, dfe[0][0][0]), reinterpret_cast<T (&) [2][3][27]>(_df_compact[first_block]));
                for( int i = 0; i < T_ARCH_DESC<T_ARCH>::T_DATA_SIZE; i++){
                    reinterpret_cast<T (&) [16]>(_dq_compact[first_block])[i] = dqe[i]; 
                }
#else
                for( int i = 0; i < T_ARCH_DESC<T_ARCH>::T_DATA_SIZE; i += 8, first_block++){
                    Block_DeDuplicate<T_ARCH>(KCAST(fArr_3_8_W, dfe[0][0][i]), _df_compact[first_block] );
                    for(int j=0; j<8; j++)
                        _dq_compact[first_block][j] = KCAST(fArr_W, dqe[i])[j];
                }               
#endif
            }





        }



    }

    };

//********************************************************************************************
//
//                   Add_Force_Differential: Master Function
//
//********************************************************************************************

template<typename T, bool enable_constraints, bool enable_muscles>  void PhysBAM::
Add_Force_Differential_With_Specialized_Kernels( T (*du_or_df_compact)[3][27],
                                                 T (*dp_or_dq_compact)[8],   
                                                 const int number_of_blocks,
                                                 const int number_of_partitions,
                                                 const SPECIALIZED_KERNEL_DATA<T,3>& specialized_data
                                                     )
{
    
#ifdef USE_THREADED_KERNELS
    int bundles = (int)(ceil(number_of_blocks / (T)SPECIALIZED_KERNEL_DATA<T,3>::VECTOR_MULT));
    Set_Force_Differential_Task<T,enable_constraints,enable_muscles>** tasks = new Set_Force_Differential_Task<T,enable_constraints,enable_muscles>*[number_of_partitions];

    for(int partition=0;
        partition<number_of_partitions;
        partition++){
        int bundle_begin=(partition)*(bundles/number_of_partitions)+min(partition,bundles%number_of_partitions);
        int bundle_end=(partition+1)*(bundles/number_of_partitions)+min((partition+1),bundles%number_of_partitions);
        tasks[partition] =
            new Set_Force_Differential_Task<T,enable_constraints,enable_muscles>(                                              
                                               bundle_begin,
                                               bundle_end-bundle_begin,
                                               du_or_df_compact,
                                               dp_or_dq_compact,
                                               du_or_df_compact, 
                                               dp_or_dq_compact, 
                                               specialized_data );

	}
    
    //long total_time = 0;
    //long extern_ticks = 0;

 {
//LOG::SCOPE scope("SPECIALIZED_KERNELS::Add_Force_Differential");
     //long start = __rdtsc();
#pragma omp parallel for
 for(int partition=0;partition<number_of_partitions;partition++){
     tasks[partition]->Run();
     delete tasks[partition];
 }
 delete[] tasks;
 //long end = __rdtsc();
 //extern_ticks = ( end - start );
 }

 //for(int partition=0;partition<number_of_partitions;partition++){
 //       total_time += tasks[partition]->_time;
 //delete tasks[partition];
 //   }

 //   delete[] tasks;

//    LOG::cout << "Extern clock ticks: " << extern_ticks << std::endl;
//    LOG::cout << "AFD clock ticks: " << total_time << std::endl;


#else
   int bundles = (int)(ceil(number_of_blocks / (T)SPECIALIZED_KERNEL_DATA<T,3>::VECTOR_MULT));
   Set_Force_Differential_Task<T,enable_constraints,enable_muscles>* task =
       new Set_Force_Differential_Task<T,enable_constraints,enable_muscles>(0,
                                                                            bundles,
                                                                            du_or_df_compact,
                                                                            dp_or_dq_compact,
                                                                            du_or_df_compact, 
                                                                            dp_or_dq_compact, 
                                                                            specialized_data);
   task->Run();
   delete task;
#endif
}


template void PhysBAM::
Add_Force_Differential_With_Specialized_Kernels<float,true,true>(
    
    float (*du_or_df_compact)[3][27],float (*dp_or_dq_compact)[8],
    const int, const int, const SPECIALIZED_KERNEL_DATA<float,3>&);

template void PhysBAM::
Add_Force_Differential_With_Specialized_Kernels<float,true,false>(
    
    float (*du_or_df_compact)[3][27],float (*dp_or_dq_compact)[8],
    const int, const int, const SPECIALIZED_KERNEL_DATA<float,3>&);

template void PhysBAM::
Add_Force_Differential_With_Specialized_Kernels<float,false,true>(
    
    float (*du_or_df_compact)[3][27],float (*dp_or_dq_compact)[8],
    const int, const int, const SPECIALIZED_KERNEL_DATA<float,3>&);

template void PhysBAM::
Add_Force_Differential_With_Specialized_Kernels<float,false,false>(
    
    float (*du_or_df_compact)[3][27],float (*dp_or_dq_compact)[8],
    const int, const int, const SPECIALIZED_KERNEL_DATA<float,3>&);


//********************************************************************************************
//
//                   Add_Force_Differential_Domain: Master Function
//
//********************************************************************************************

template<typename T, bool enable_constraints, bool enable_muscles>  void PhysBAM::
Add_Force_Differential_With_Specialized_Kernels_Domain( T (*du_or_df_compact)[3][27],
                                                        T (*dp_or_dq_compact)[8],   
                                                        const ARRAY<int>& blocks,
                                                        const int number_of_partitions,
                                                        const SPECIALIZED_KERNEL_DATA<T,3>& specialized_data
                                                        )
{
    //LOG::SCOPE scope("SPECIALIZED_KERNELS::Add_Force_Differential");
    ARRAY<int> bundle_list;
    for(int block_idx = 1; block_idx <= blocks.m; block_idx++){
        int block = blocks(block_idx);
        int bundle = ((block-1)/(T)SPECIALIZED_KERNEL_DATA<T,3>::VECTOR_MULT)+1;
        bundle_list.Append(bundle);
    }
    bundle_list.Prune_Duplicates();
    
#ifdef USE_THREADED_KERNELS
    Set_Force_Differential_Task<T,enable_constraints,enable_muscles>** tasks = new Set_Force_Differential_Task<T,enable_constraints,enable_muscles>*[number_of_partitions];
    
    for(int partition=0;
        partition<number_of_partitions;
        partition++){
        int bundle_begin=(partition)*(bundle_list.m/number_of_partitions)+min(partition,bundle_list.m%number_of_partitions);
        int bundle_end=(partition+1)*(bundle_list.m/number_of_partitions)+min((partition+1),bundle_list.m%number_of_partitions);
        tasks[partition] =
            new Set_Force_Differential_Task<T,enable_constraints,enable_muscles>(                                              
                                               bundle_begin,
                                               bundle_end-bundle_begin,
                                               bundle_list,
                                               du_or_df_compact,
                                               dp_or_dq_compact,
                                               du_or_df_compact, 
                                               dp_or_dq_compact, 
                                               specialized_data );

	}

#pragma omp parallel for num_threads(number_of_partitions)
    for(int partition=0;partition<number_of_partitions;partition++){
        tasks[partition]->Run();
        delete tasks[partition];
    }
    delete[] tasks;

#else
   Set_Force_Differential_Task<T,enable_constraints,enable_muscles>* task =
       new Set_Force_Differential_Task<T,enable_constraints,enable_muscles>(0,
                                                                            bundle_list.m,
                                                                            bundle_list,
                                                                            du_or_df_compact,
                                                                            dp_or_dq_compact,
                                                                            du_or_df_compact, 
                                                                            dp_or_dq_compact, 
                                                                            specialized_data);
   task->Run();
   delete task;
#endif
}


template void PhysBAM::
Add_Force_Differential_With_Specialized_Kernels_Domain<float,true,true>(   
    float (*du_or_df_compact)[3][27],float (*dp_or_dq_compact)[8],
    const ARRAY<int>&, const int, const SPECIALIZED_KERNEL_DATA<float,3>&);

template void PhysBAM::
Add_Force_Differential_With_Specialized_Kernels_Domain<float,true,false>(   
    float (*du_or_df_compact)[3][27],float (*dp_or_dq_compact)[8],
    const ARRAY<int>&, const int, const SPECIALIZED_KERNEL_DATA<float,3>&);

template void PhysBAM::
Add_Force_Differential_With_Specialized_Kernels_Domain<float,false,true>(
    float (*du_or_df_compact)[3][27],float (*dp_or_dq_compact)[8],
    const ARRAY<int>&, const int, const SPECIALIZED_KERNEL_DATA<float,3>&);

template void PhysBAM::
Add_Force_Differential_With_Specialized_Kernels_Domain<float,false,false>(   
    float (*du_or_df_compact)[3][27],float (*dp_or_dq_compact)[8],
    const ARRAY<int>&, const int, const SPECIALIZED_KERNEL_DATA<float,3>&);


//********************************************************************************************
//
//                   Add Force First Order: Threaded Task
//
//********************************************************************************************

    template<class T,  bool enable_constraints, bool enable_muscles>
        class Add_Force_First_Order_Task : public PhysBAM::PTHREAD_QUEUE::TASK
        {
        private:
            typedef VECTOR<T,3> TV;

            int _bundle_start;
            int _bundle_count;

            T (*_u_compact)[3][27];
            T (*_p_compact)[8];   
            T (*_f_compact)[3][27];
            T (*_q_compact)[8];
            SPECIALIZED_KERNEL_DATA<T,3>& _specialized_data;


        public:
//#####################################################################
// Constructor
//#####################################################################

            Add_Force_First_Order_Task(
              int bundle_start,
              int bundle_count,                              
              T (*u_compact)[3][27],
              T (*p_compact)[8],   
              T (*f_compact)[3][27],
              T (*q_compact)[8],
              SPECIALIZED_KERNEL_DATA<T,3>& specialized_data
                                       )
                : 
                _bundle_start(bundle_start),
                _bundle_count(bundle_count),
                _u_compact(u_compact),                             
                _p_compact(p_compact),                             
                _f_compact(f_compact),                             
                _q_compact(q_compact),                             
                _specialized_data(specialized_data)                             
            {  
            }


//#####################################################################
// Function: Run
//#####################################################################

void Run()
    {
        typedef SELECTED_ARCH T_ARCH;

    STATIC_ASSERT(T_ARCH_DESC<T_ARCH>::T_DATA_SIZE == SPECIALIZED_KERNEL_DATA<T,3>::VECTOR_WIDTH);

        for( int _bundle = 0; _bundle < _bundle_count; _bundle++)
        {
            int bundle = _bundle + _bundle_start + 1;
                       
            typedef T pressure_p[T_ARCH_DESC<T_ARCH>::T_DATA_SIZE];
            
            //__declspec(align(SPECIALIZED_KERNEL_DATA<T,3>::ALIGNMENT)) T ue[3][8][T_ARCH_DESC<T_ARCH>::T_DATA_SIZE];
            //__declspec(align(SPECIALIZED_KERNEL_DATA<T,3>::ALIGNMENT)) T fe[3][8][T_ARCH_DESC<T_ARCH>::T_DATA_SIZE];
            __attribute__ ((aligned( BUNDLE_WIDTH * 4))) T ue[3][8][T_ARCH_DESC<T_ARCH>::T_DATA_SIZE];
            __attribute__ ((aligned( BUNDLE_WIDTH * 4))) T fe[3][8][T_ARCH_DESC<T_ARCH>::T_DATA_SIZE];
            memset(fe,0,sizeof(T)*3*8*T_ARCH_DESC<T_ARCH>::T_DATA_SIZE);
            //__declspec(align(SPECIALIZED_KERNEL_DATA<T,3>::ALIGNMENT)) pressure_p pe;
            //__declspec(align(SPECIALIZED_KERNEL_DATA<T,3>::ALIGNMENT)) pressure_p qe;
            __attribute__ ((aligned( BUNDLE_WIDTH * 4))) pressure_p pe;
            __attribute__ ((aligned( BUNDLE_WIDTH * 4))) pressure_p qe;
            memset(qe,0,sizeof(T)*T_ARCH_DESC<T_ARCH>::T_DATA_SIZE);

            int first_block = (bundle-1)*SPECIALIZED_KERNEL_DATA<T,3>::VECTOR_MULT;
            for( int i = 0; i < T_ARCH_DESC<T_ARCH>::T_DATA_SIZE; i += 8, first_block++){
                Block_Duplicate<T_ARCH>(_u_compact[first_block], KCAST(fArr_3_8_W, ue[0][0][i]));
                for(int j=0; j<8; j++)
                    KCAST(fArr_W, pe[i])[j] = _p_compact[first_block][j];
            }

                for( int i = 0; i < T_ARCH_DESC<T_ARCH>::T_DATA_SIZE; i += T_ARCH_DESC<T_ARCH>::WIDTH)
                    ::Add_Force_First_Order<
                        T_ARCH_DESC<T_ARCH>::T_MATERIAL,
                        T_ARCH,
                        typename  T_ARCH_DESC<T_ARCH>::T_DATA,
                        typename T_ARCH_DESC<T_ARCH>::I_DATA>::Run(
                        KCAST( cfArr_3_8_W, ue[0][0][i]), 
                        KCAST( cfArr_W, pe[i]),
                        KCAST( cfArr_W, (_specialized_data.mu_bundled)(bundle).data[i]),
                        KCAST( cfArr_W, (_specialized_data.alpha_bundled)(bundle).data[i]),
                        KCAST( cfArr_W, (_specialized_data.alpha_sqr_over_kappa_bundled)(bundle).data[i]),
                        KCAST( cfArr_W, (_specialized_data.kappa_bundled)(bundle).data[i]),
                        KCAST( cfArr_W, (_specialized_data.one_over_h_bundled)(1).data[i]),
                        KCAST( cfArr_W, (_specialized_data.cell_volume_bundled)(1).data[i]),
                        KCAST( cfArr_9_W, (_specialized_data.U_bundled)(bundle).data[0][i]),
                        KCAST( cfArr_9_W, (_specialized_data.V_bundled)(bundle).data[0][i]),
                        KCAST( cfArr_3_W, (_specialized_data.Sigma_bundled)(bundle).data[0][i]),
                        KCAST( cfArr_3_W, (_specialized_data.Q_hat_bundled)(bundle).data[0][i]),
                        KCAST( fArr_3_W, (_specialized_data.P_hat_bundled)(bundle).data[0][i]),
                        KCAST( fArr_3_8_W, fe[0][0][i]),
                        KCAST( fArr_W, qe[i]));
                
                for( int i = 0; i < T_ARCH_DESC<T_ARCH>::T_DATA_SIZE; i += T_ARCH_DESC<T_ARCH>::WIDTH)
                    ::Force_Stabilization<
                        T_ARCH,
                        typename T_ARCH_DESC<T_ARCH>::T_DATA,
                        typename T_ARCH_DESC<T_ARCH>::I_DATA>(
                        KCAST( cfArr_3_8_W, ue[0][0][i]), 
                        KCAST( cfArr_W, (_specialized_data.mu_stab_bundled)(bundle).data[i]),
                        KCAST( fArr_3_8_W, fe[0][0][i]));
                
                if(enable_muscles){
                    int muscle_bundle_offset = (_specialized_data.muscle_base_offsets)(bundle);
                    
                    for( int layer=0; layer < (_specialized_data.muscle_base_lengths)(bundle); layer++)
                        {
                            int muscle_bundle = muscle_bundle_offset + layer;
                            
                            for( int i = 0; i < T_ARCH_DESC<T_ARCH>::T_DATA_SIZE; i += T_ARCH_DESC<T_ARCH>::WIDTH)
                                ::Muscle_Forces<
                                    T_ARCH,
                                    typename T_ARCH_DESC<T_ARCH>::T_DATA,
                                    typename T_ARCH_DESC<T_ARCH>::I_DATA>(
                                    KCAST( fArr_3_8_W, fe[0][0][i]),
                                    KCAST( cfArr_3_W, (_specialized_data.muscle_fiber)(muscle_bundle).data[0][i]),
                                    KCAST( cfArr_3_W, (_specialized_data.muscle_F_fiber)(muscle_bundle).data[0][i]),
                                    KCAST( cfArr_W, (_specialized_data.muscle_c1)(muscle_bundle).data[i]),
                                    KCAST( cfArr_W, (_specialized_data.one_over_h_bundled)(1).data[i]),
                                    KCAST( cfArr_W, (_specialized_data.cell_volume_bundled)(1).data[i]));
                        }
                }
                
                if(enable_constraints){
                    int constraint_bundle_offset = (_specialized_data.constraint_base_offsets)(bundle);
                    
                    for( int layer=0; layer < (_specialized_data.constraint_base_lengths)(bundle); layer++)
                        {
                            int constraint_bundle = constraint_bundle_offset + layer;
                            for( int i = 0; i < T_ARCH_DESC<T_ARCH>::T_DATA_SIZE; i += T_ARCH_DESC<T_ARCH>::WIDTH) {
                                ::Collision_Forces<T_ARCH,typename T_ARCH_DESC<T_ARCH>::T_DATA,typename T_ARCH_DESC<T_ARCH>::I_DATA>(
                                                   KCAST( fArr_3_8_W, fe[0][0][i]),
                                                   KCAST( cfArr_3_8_W, ue[0][0][i]), 
                                                   KCAST( cfArr_3_W, (_specialized_data.constraint_node_positions)(constraint_bundle).data[0][i]),
                                                   KCAST( cfArr_3_W, (_specialized_data.constraint_multilinear_coordinates)(constraint_bundle).data[0][i]),
                                                   KCAST( cfArr_W, (_specialized_data.h_bundled)(1).data[i]),
                                                   KCAST( ciArr_W, (_specialized_data.spring_id)(constraint_bundle).data[i]),
                                                   KCAST( ciArr_W, (_specialized_data.spring_id_X)(constraint_bundle).data[i]),
                                                   KCAST( ciArr_W, (_specialized_data.spring_id_Y)(constraint_bundle).data[i]),
                                                   KCAST( ciArr_W, (_specialized_data.spring_id_Z)(constraint_bundle).data[i]),
                                                   (float*)((_specialized_data.collision_spring_locations)->Get_Array_Pointer()),
                                                   (_specialized_data.collision_spring_constants)->Get_Array_Pointer());                                
                                        }
                             }
                      }

                first_block = (bundle-1)*SPECIALIZED_KERNEL_DATA<T,3>::VECTOR_MULT;
                for( int i = 0; i < T_ARCH_DESC<T_ARCH>::T_DATA_SIZE; i += 8, first_block++){
                    Block_DeDuplicate<T_ARCH>(KCAST(fArr_3_8_W, fe[0][0][i]), _f_compact[first_block] );
                    for(int j=0; j<8; j++)
                        _q_compact[first_block][j] = KCAST(fArr_W, qe[i])[j];
                }   
                
            }
    }


    };



//********************************************************************************************
//
//                   Add_Force_Differential: Master Function
//
//********************************************************************************************


template<typename T, bool enable_constraints, bool enable_muscles>  void PhysBAM::
Add_Force_First_Order_With_Specialized_Kernels( T (*u_or_f_compact)[3][27],
                                                T (*p_or_q_compact)[8], 
                                                const int number_of_blocks,
                                                const int number_of_partitions,
                                                SPECIALIZED_KERNEL_DATA<T,3>& specialized_data
                                                     )

{
    //LOG::SCOPE scope("SPECIALIZED_KERNELS::Add_Force");
#ifdef USE_THREADED_KERNELS
    int bundles = (int)(ceil(number_of_blocks / (T)SPECIALIZED_KERNEL_DATA<T,3>::VECTOR_MULT));  
    Add_Force_First_Order_Task<T,enable_constraints,enable_muscles>** tasks = new Add_Force_First_Order_Task<T,enable_constraints,enable_muscles>*[number_of_partitions];

    for(int partition=0;
        partition<number_of_partitions;
        partition++){
        int bundle_begin=(partition)*(bundles/number_of_partitions)+min(partition,bundles%number_of_partitions);
        int bundle_end=(partition+1)*(bundles/number_of_partitions)+min((partition+1),bundles%number_of_partitions);
        tasks[partition] =
            new Add_Force_First_Order_Task<T,enable_constraints,enable_muscles>(                                              
                                               bundle_begin,
                                               bundle_end-bundle_begin,
                                               u_or_f_compact,
                                               p_or_q_compact,
                                               u_or_f_compact,
                                               p_or_q_compact,
                                               specialized_data);

	}

#pragma omp parallel for num_threads(number_of_partitions)
    for(int partition=0;partition<number_of_partitions;partition++){
        tasks[partition]->Run();
        delete tasks[partition];
    }
    delete[] tasks;

#else
    int bundles = (int)(ceil(number_of_blocks / (T)SPECIALIZED_KERNEL_DATA<T,3>::VECTOR_MULT));
    Add_Force_First_Order_Task<T,enable_constraints,enable_muscles>* task =
        new Add_Force_First_Order_Task<T,enable_constraints,enable_muscles>(0,
                                                                            bundles,
                                                                            u_or_f_compact,
                                                                            p_or_q_compact,
                                                                            u_or_f_compact,
                                                                            p_or_q_compact,
                                                                            specialized_data);
    task->Run();
    delete task;
#endif

}


template void PhysBAM::
Add_Force_First_Order_With_Specialized_Kernels<float,true,true>(
    float (*u_or_f_compact)[3][27],float (*p_or_q_compact)[8],
    const int, const int, SPECIALIZED_KERNEL_DATA<float,3>&);


template void PhysBAM::
Add_Force_First_Order_With_Specialized_Kernels<float,true,false>(
    float (*u_or_f_compact)[3][27],float (*p_or_q_compact)[8],
    const int, const int, SPECIALIZED_KERNEL_DATA<float,3>&);


template void PhysBAM::
Add_Force_First_Order_With_Specialized_Kernels<float,false,true>(
    float (*u_or_f_compact)[3][27],float (*p_or_q_compact)[8],
    const int, const int, SPECIALIZED_KERNEL_DATA<float,3>&);


template void PhysBAM::
Add_Force_First_Order_With_Specialized_Kernels<float,false,false>(
    float (*u_or_f_compact)[3][27],float (*p_or_q_compact)[8],
    const int, const int, SPECIALIZED_KERNEL_DATA<float,3>&);

#endif

#ifndef __SPECIALIZED_KERNELS_H__
#define __SPECIALIZED_KERNELS_H__

#include "SPECIALIZED_KERNELS_DATA.h"
#include "NONLINEAR_ELASTICITY.h"
#include <Common/ALIGNED_ARRAY.h>
#include <Thread_Queueing/PTHREAD_QUEUE.h>

#ifdef USE_SPECIALIZED_KERNELS

namespace PhysBAM
{

    template<typename T, bool enable_constraints, bool enable_muscles>
    void Update_Position_Based_State_With_Specialized_Kernels(
        T (*u_compact)[3][27],
        T (*p_compact)[8], 
        T (*d_compact)[3][27],
        const int number_of_blocks,
        const int number_of_partitions,             
        const ARRAY<T>& muscle_fiber_max_stresses,
        const ARRAY<T>& muscle_activations,
        SPECIALIZED_KERNEL_DATA<T,3>& specialized_data);
    
    template<typename T,  bool enable_constraints, bool enable_muscles>
        void Add_Force_Differential_With_Specialized_Kernels( 
            T (*du_or_df_compact)[3][27],
            T (*dp_or_dq_compact)[8],  
            const int number_of_blocks,
            const int number_of_partitions,
            const SPECIALIZED_KERNEL_DATA<T,3>& specialized_data);  

    template<typename T,  bool enable_constraints, bool enable_muscles>
        void Add_Force_Differential_With_Specialized_Kernels_Domain( 
            T (*du_or_df_compact)[3][27],
            T (*dp_or_dq_compact)[8],  
            const ARRAY<int>& domain_blocks,
            const int number_of_partitions,
            const SPECIALIZED_KERNEL_DATA<T,3>& specialized_data);  
    
    template<typename T,  bool enable_constraints, bool enable_muscles>
        void Add_Force_First_Order_With_Specialized_Kernels(
            T (*u_or_f_compact)[3][27],
            T (*p_or_q_compact)[8], 
            const int number_of_blocks,
            const int number_of_partitions,
            SPECIALIZED_KERNEL_DATA<T,3>& specialized_data);


};


#endif

#endif

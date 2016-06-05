//#####################################################################
// Copyright 2002-2007, Robert Bridson, Ronald Fedkiw, Eran Guendelman, Geoffrey Irving, Nipun Kwatra, Frank Losasso, Duc Nguyen, Craig Schroeder, Andrew Selle, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ARRAY
//#####################################################################
#if 1
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#else
#ifndef __ARRAYS_ND__
#define __ARRAYS_ND__

#include <stdlib.h>
#include <malloc.h>
#include <stdint.h>
#include <errno.h>
#include <PhysBAM_Tools/Arrays_Computations/ARRAY_COPY.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND_BASE.h>
#include <PhysBAM_Tools/Math_Tools/ONE.h>
#include <PhysBAM_Tools/Math_Tools/RANGE.h>
#include <PhysBAM_Tools/Utilities/EXCEPTIONS.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
namespace PhysBAM{

#define ALIGNMENT 64
//#define DEBUG_ALIGNMENT

    template<class T,int d>
        class ARRAY<T,VECTOR<int,d> >: public ARRAYS_ND_BASE<VECTOR<T,d> >
        {
        public:
            typedef VECTOR<int,d> TV_INT;
            enum WORKAROUND1 {dimension=d};
            template<class T2> struct REBIND{typedef ARRAY<T2,TV_INT> TYPE;};
            typedef T ELEMENT;typedef TV_INT INDEX;typedef T& RESULT_TYPE;typedef const T& CONST_RESULT_TYPE;
            typedef ARRAYS_ND_BASE<VECTOR<T,d> > BASE;

            using BASE::array; // one-dimensional data storage
            using BASE::domain;using BASE::counts;
            using BASE::Calculate_Acceleration_Constants;
        private:
            using BASE::base_pointer;
            void* raw_ptr; 
        public:

        ARRAY()
            :BASE(), raw_ptr(NULL)
                {
#ifdef DEBUG_ALIGNMENT
                    LOG::cout << "Created NULL Aligned Array ND at: " << std::hex << (uintptr_t)(raw_ptr)
                              << std::dec <<std::endl;
#endif
                    Calculate_Acceleration_Constants();
                }
    
        ARRAY(const RANGE<TV_INT>& domain_input,const bool initialize_using_default_constructor=true)
            :BASE(domain_input), raw_ptr(NULL)
                {
                    assert(counts.Min()>=0);int size=counts.Product();
                    int error=posix_memalign(&raw_ptr, ALIGNMENT, Value(size)*sizeof(T));
#ifdef DEBUG_ALIGNMENT
                    LOG::cout << "Created Aligned Array ND at: " << std::hex << (uintptr_t)(raw_ptr)
                              << std::dec <<std::endl;
#endif
                    switch(error){
                    case EINVAL:
                        PHYSBAM_FATAL_ERROR("Could not align to a non-power of two or non-multiple of pointer size.");
                    case ENOMEM:
                        PHYSBAM_FATAL_ERROR("Insufficient memory to perform aligned memory allocation.");
                    default:
                        break;
                    }        
                    
                    ARRAY_VIEW<T> new_array(size,new (raw_ptr) T[size]);new_array.Exchange(array);
                    Calculate_Acceleration_Constants();
                    if(initialize_using_default_constructor) ARRAYS_COMPUTATIONS::Fill(array,T()); // initialize array using default constructor
                }

        ARRAY(const int m_start_input,const int m_end_input,const int n_start_input,const int n_end_input,const int mn_start_input,const int mn_end_input,
              const bool initialize_using_default_constructor=true)
            :BASE(RANGE<TV_INT>(TV_INT(m_start_input,n_start_input,mn_start_input),TV_INT(m_end_input,n_end_input,mn_end_input))), raw_ptr(NULL)
                {
                    assert(counts.Min()>=0);int size=counts.Product();
                    int error=posix_memalign(&raw_ptr, ALIGNMENT, Value(size)*sizeof(T));
#ifdef DEBUG_ALIGNMENT
                    LOG::cout << "Created Aligned Array ND at: " << std::hex << (uintptr_t)(raw_ptr)
                              << std::dec <<std::endl;
#endif
                    switch(error){
                    case EINVAL:
                        PHYSBAM_FATAL_ERROR("Could not align to a non-power of two or non-multiple of pointer size.");
                    case ENOMEM:
                        PHYSBAM_FATAL_ERROR("Insufficient memory to perform aligned memory allocation.");
                    default:
                        break;
                    }            
                    ARRAY_VIEW<T> new_array(size,new (raw_ptr) T[size]);new_array.Exchange(array);
                    Calculate_Acceleration_Constants();
                    if(initialize_using_default_constructor) ARRAYS_COMPUTATIONS::Fill(array,T()); // initialize array using default constructor
                }

        ARRAY(const int m_start_input,const int m_end_input,const int n_start_input,const int n_end_input,const bool initialize_using_default_constructor=true)
            :BASE(RANGE<TV_INT>(TV_INT(m_start_input,n_start_input),TV_INT(m_end_input,n_end_input))), raw_ptr(NULL)
                {
                    assert(counts.Min()>=0);int size=counts.Product();
                    int error=posix_memalign(&raw_ptr, ALIGNMENT, Value(size)*sizeof(T));
#ifdef DEBUG_ALIGNMENT
                    LOG::cout << "Created Aligned Array ND at: " << std::hex << (uintptr_t)(raw_ptr)
                              << std::dec <<std::endl;
#endif
                    switch(error){
                    case EINVAL:
                        PHYSBAM_FATAL_ERROR("Could not align to a non-power of two or non-multiple of pointer size.");
                    case ENOMEM:
                        PHYSBAM_FATAL_ERROR("Insufficient memory to perform aligned memory allocation.");
                    default:
                        break;
                    }  
                    ARRAY_VIEW<T> new_array(size,new (raw_ptr) T[size]);new_array.Exchange(array);
                    Calculate_Acceleration_Constants();
                    if(initialize_using_default_constructor) ARRAYS_COMPUTATIONS::Fill(array,T()); // initialize array using default constructor
                }

        ARRAY(const int m_start_input,const int m_end_input,const bool initialize_using_default_constructor=true)
            :BASE(RANGE<TV_INT>(TV_INT(m_start_input),TV_INT(m_end_input))), raw_ptr(NULL)
                {
                    assert(counts.Min()>=0);int size=counts.Product();
                    int error=posix_memalign(&raw_ptr, ALIGNMENT, Value(size)*sizeof(T));
#ifdef DEBUG_ALIGNMENT
                    LOG::cout << "Created Aligned Array ND at: " << std::hex << (uintptr_t)(raw_ptr)
                              << std::dec <<std::endl;
#endif
                    switch(error){
                    case EINVAL:
                        PHYSBAM_FATAL_ERROR("Could not align to a non-power of two or non-multiple of pointer size.");
                    case ENOMEM:
                        PHYSBAM_FATAL_ERROR("Insufficient memory to perform aligned memory allocation.");
                    default:
                        break;
                    }  
                    ARRAY_VIEW<T> new_array(size,new (raw_ptr) T[size]);new_array.Exchange(array);
                    Calculate_Acceleration_Constants();
                    if(initialize_using_default_constructor) ARRAYS_COMPUTATIONS::Fill(array,T()); // initialize array using default constructor
                }

        ARRAY(const ARRAY& old_array,const bool initialize_with_old_array=true)
            :BASE(old_array.domain), raw_ptr(NULL)
                {
                    int size=old_array.array.Size();
                    int error=posix_memalign(&raw_ptr, ALIGNMENT, Value(size)*sizeof(T));
#ifdef DEBUG_ALIGNMENT
                    LOG::cout << "Created Aligned Array ND at: " << std::hex << (uintptr_t)(raw_ptr)
                              << std::dec <<std::endl;
#endif
                    switch(error){
                    case EINVAL:
                        PHYSBAM_FATAL_ERROR("Could not align to a non-power of two or non-multiple of pointer size.");
                    case ENOMEM:
                        PHYSBAM_FATAL_ERROR("Insufficient memory to perform aligned memory allocation.");
                    default:
                        break;
                    }  
                    ARRAY_VIEW<T> new_array(size,new (raw_ptr) T[size]);new_array.Exchange(array);
                    Calculate_Acceleration_Constants();
                    if(initialize_with_old_array) array=old_array.array;
                }

            ~ARRAY()
                {
                    if(raw_ptr)
                        free(raw_ptr);
                    raw_ptr=NULL;
                }

            void Clean_Memory()
            {Resize(RANGE<TV_INT>(TV_INT::All_Ones_Vector(),TV_INT()),false,false);}

            void Delete_Pointers_And_Clean_Memory() // only valid if T is a pointer type
            {for(int i=1;i<=array.Size();i++) delete array(i);Clean_Memory();}

            ARRAY& operator=(const ARRAY& source)
                {
                    if(array.Size()!=source.array.Size()){
                        if(raw_ptr)
                            free(raw_ptr);
                        int error=posix_memalign(&raw_ptr, ALIGNMENT, Value(source.array.Size())*sizeof(T));
#ifdef DEBUG_ALIGNMENT
                        LOG::cout << "Created Aligned Array ND at: " << std::hex << (uintptr_t)(raw_ptr)
                        << std::dec <<std::endl;
#endif
                        switch(error){
                        case EINVAL:
                            PHYSBAM_FATAL_ERROR("Could not align to a non-power of two or non-multiple of pointer size.");
                        case ENOMEM:
                            PHYSBAM_FATAL_ERROR("Insufficient memory to perform aligned memory allocation.");
                        default:
                            break;
                        }  
                        ARRAY_VIEW<T> new_array(source.array.Size(),new (raw_ptr) T[source.array.Size()]);new_array.Exchange(array);
                    }
                    else if(this==&source) return *this;
                    counts=source.counts;domain=source.domain;
                    Calculate_Acceleration_Constants();
                    array=source.array;return *this;
                }
            
            template<class T_ARRAY2>
                ARRAY& operator=(const T_ARRAY2& source)
                {STATIC_ASSERT(IS_SAME<ELEMENT,typename T_ARRAY2::ELEMENT>::value);
                 if(counts!=source.Size()){
                     if(raw_ptr)
                         free(raw_ptr);
                     int source_array_size=source.Size().Product();
                     int error=posix_memalign(&raw_ptr, ALIGNMENT, Value(source_array_size)*sizeof(T));
#ifdef DEBUG_ALIGNMENT
                     LOG::cout << "Created Aligned Array ND at: " << std::hex << (uintptr_t)(raw_ptr)
                               << std::dec <<std::endl;
#endif
                     switch(error){
                     case EINVAL:
                         PHYSBAM_FATAL_ERROR("Could not align to a non-power of two or non-multiple of pointer size.");
                     case ENOMEM:
                         PHYSBAM_FATAL_ERROR("Insufficient memory to perform aligned memory allocation.");
                     default:
                         break;
                     }  
                     ARRAY_VIEW<T> new_array(source_array_size,new (raw_ptr) T[source_array_size]);new_array.Exchange(array);}
                 counts=source.Size();domain=source.Domain_Indices();
                 Calculate_Acceleration_Constants();
                 ARRAY_BASE<T,BASE,TV_INT>::operator=(source);return *this;
                }


            void Resize(int m_start_new,int m_end_new,int n_start_new,int n_end_new,int mn_start_new,int mn_end_new,const bool initialize_new_elements=true,
                        const bool copy_existing_elements=true,const T& initialization_value=T())
            {STATIC_ASSERT(d==3);RANGE<TV_INT> box(m_start_new,m_end_new,n_start_new,n_end_new,mn_start_new,mn_end_new);
                Resize(box,initialize_new_elements,copy_existing_elements,initialization_value);}

            void Resize(int m_start_new,int m_end_new,int n_start_new,int n_end_new,const bool initialize_new_elements=true,
                        const bool copy_existing_elements=true,const T& initialization_value=T())
            {STATIC_ASSERT(d==2);RANGE<TV_INT> box(m_start_new,m_end_new,n_start_new,n_end_new);Resize(box,initialize_new_elements,copy_existing_elements,initialization_value);}

            void Resize(int m_start_new,int m_end_new,const bool initialize_new_elements=true,
                        const bool copy_existing_elements=true,const T& initialization_value=T())
            {STATIC_ASSERT(d==1);RANGE<TV_INT> box(m_start_new,m_end_new);Resize(box,initialize_new_elements,copy_existing_elements,initialization_value);}

        private:
            void Resize_Helper(const RANGE<VECTOR<int,3> >& box,ARRAY_VIEW<T> array_new,const TV_INT& counts_new)
            {TV_INT m1=TV_INT::Componentwise_Max(domain.min_corner,box.min_corner),m2=TV_INT::Componentwise_Min(domain.max_corner,box.max_corner),i;
                for(i.x=m1.x;i.x<=m2.x;i.x++) for(i.y=m1.y;i.y<=m2.y;i.y++) for(i.z=m1.z;i.z<=m2.z;i.z++){
                            TV_INT diff_old(i-domain.min_corner),diff_new(i-box.min_corner);
                            array_new(((diff_new.x*counts_new.y+diff_new.y)*counts_new.z+diff_new.z)+1)=array(((diff_old.x*counts.y+diff_old.y)*counts.z+diff_old.z)+1);}}

            void Resize_Helper(const RANGE<VECTOR<int,2> >& box,ARRAY_VIEW<T> array_new,const TV_INT& counts_new)
            {TV_INT m1=TV_INT::Componentwise_Max(domain.min_corner,box.min_corner),m2=TV_INT::Componentwise_Min(domain.max_corner,box.max_corner),i;
                for(i.x=m1.x;i.x<=m2.x;i.x++) for(i.y=m1.y;i.y<=m2.y;i.y++){
                        TV_INT diff_old(i-domain.min_corner),diff_new(i-box.min_corner);
                        array_new((diff_new.x*counts_new.y+diff_new.y)+1)=array((diff_old.x*counts.y+diff_old.y)+1);}}

            void Resize_Helper(const RANGE<VECTOR<int,1> >& box,ARRAY_VIEW<T> array_new,const TV_INT& counts_new)
            {TV_INT m1=TV_INT::Componentwise_Max(domain.min_corner,box.min_corner),m2=TV_INT::Componentwise_Min(domain.max_corner,box.max_corner),i;
                for(i.x=m1.x;i.x<=m2.x;i.x++){
                    TV_INT diff_old(i-domain.min_corner),diff_new(i-box.min_corner);
                    array_new(diff_new.x+1)=array(diff_old.x+1);}}
        public:

            void Resize(const RANGE<TV_INT>& box,const bool initialize_new_elements=true,const bool copy_existing_elements=true,const T& initialization_value=T())
            {
                if(box==domain) return;
                TV_INT counts_new(box.Edge_Lengths()+1);
                assert(counts_new.Min()>=0);
                int size_new=counts_new.Product();
                void* tmp_ptr;
                int error=posix_memalign(&tmp_ptr, ALIGNMENT, Value(size_new)*sizeof(T));
#ifdef DEBUG_ALIGNMENT
                LOG::cout << "Resized Aligned Array ND at: " << std::hex << (uintptr_t)(tmp_ptr)
                          << std::dec <<std::endl;
#endif
                switch(error){
                case EINVAL:
                    PHYSBAM_FATAL_ERROR("Could not align to a non-power of two or non-multiple of pointer size.");
                case ENOMEM:
                    PHYSBAM_FATAL_ERROR("Insufficient memory to perform aligned memory allocation.");
                default:
                    break;
                }  
                ARRAY_VIEW<T> array_new(size_new,new (tmp_ptr) T[size_new]);
                if(initialize_new_elements) ARRAYS_COMPUTATIONS::Fill(array_new,initialization_value);
                if(copy_existing_elements) Resize_Helper(box,array_new,counts_new);
                domain=box;counts=counts_new;
                if(raw_ptr)
                    free(raw_ptr);
                raw_ptr=tmp_ptr;
                array.Exchange(array_new);Calculate_Acceleration_Constants();
            }
    
            // TODO: this function is extremely broken, since all the functions in the class depends on size being correct
            void Resize_In_Place(const int m_start_new,const int m_end_new,const int n_start_new,const int n_end_new,const int mn_start_new,const int mn_end_new)
            {RANGE<TV_INT> box(m_start_new,m_end_new,n_start_new,n_end_new,mn_start_new,mn_end_new);Resize_In_Place(box);}
    
            void Resize_In_Place(const RANGE<TV_INT>& box)
            {TV_INT counts_new(box.Edge_Lengths()+1);
                if(array.Size()>=counts_new.Product()){domain=box;counts=counts_new;Calculate_Acceleration_Constants();}
                else Resize(box,false,false);}

            static void Exchange_Arrays(ARRAY& a,ARRAY& b)
            {a.array.Exchange(b.array);
                exchange(a.domain,b.domain);exchange(a.counts,b.counts);
                a.Calculate_Acceleration_Constants();b.Calculate_Acceleration_Constants();}
//#####################################################################
        };
}

#undef ALIGNMENT

#endif
#endif

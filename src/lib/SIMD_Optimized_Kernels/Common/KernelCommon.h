//#####################################################################
//  Copyright (c) 2011-2013 Nathan Mitchell, Eftychios Sifakis.
//  This file is covered by the FreeBSD license. Please refer to the 
//  license.txt file for more information.
//#####################################################################


#ifndef __KERNEL_COMMON_H__
#define __KERNEL_COMMON_H__

#include <type_traits>

//
//   Don't define here. Use FIB=1 on the makefile command line.
//
//#define FORCE_IDENTICAL_BEHAVIOR

//
//   Enable use of c++ std io to print Number's
//
#define ENABLE_IO_SUPPORT

//
//   Enable INTEL VECTOR ARCHETECTURES
//
//#define ENABLE_SSE_INSTRUCTION_SET
//#define ENABLE_AVX_INSTRUCTION_SET
//#define ENABLE_MIC_INSTRUCTION_SET

//
//   Enable ARM VECTOR ARCHETECTURES
//
//#define ENABLE_NEON_INSTRUCTION_SET




#if defined(ENABLE_SSE_INSTRUCTION_SET)
#include <immintrin.h>
#ifndef __INTEL_COMPILER
#include <pmmintrin.h>
#endif
#endif

#if defined(ENABLE_AVX_INSTRUCTION_SET)
#include <immintrin.h>
#ifndef __INTEL_COMPILER
#include <pmmintrin.h>
#endif
#endif

#if defined(ENABLE_MIC_INSTRUCTION_SET)
#include <immintrin.h>
#include <zmmintrin.h>
#endif

#if defined(ENABLE_NEON_INSTRUCTION_SET)
#include <arm_neon.h>
#endif


#include "NumberPolicy.h"
#include "Mask.h"
#include "Number.h"
#include "Discrete.h"
#include "Vector3.h"
#include "Constants.h"


// Load Architecture Specific Versions of Number

template <typename Vtype>
struct VTYPE_POLICY{
    const static int V_WIDTH=0;
};

template <>
struct VTYPE_POLICY<float>{
    const static int V_WIDTH=1;
};

#if defined(ENABLE_SSE_INSTRUCTION_SET)
template <>
struct VTYPE_POLICY<__m128>{
    const static int V_WIDTH=4;
};
#include "arch/x86_64/Mask.SSE.h"
#include "arch/x86_64/Number.SSE.h"
#include "arch/x86_64/Discrete.SSE.h"
#endif


#if defined(ENABLE_AVX_INSTRUCTION_SET)
template <>
struct VTYPE_POLICY<__m256>{
    const static int V_WIDTH=8;
};
#include "arch/x86_64/Mask.AVX.h"
#include "arch/x86_64/Number.AVX.h"
#include "arch/x86_64/Discrete.AVX.h"
#endif


#if defined(ENABLE_MIC_INSTRUCTION_SET)
template <>
struct VTYPE_POLICY<__m512>{
    const static int V_WIDTH=16;
};
#include "arch/mic/Mask.MIC.h"
#include "arch/mic/Number.MIC.h"
#include "arch/mic/Discrete.MIC.h"
#endif


#if defined(ENABLE_NEON_INSTRUCTION_SET)
template <>
struct VTYPE_POLICY<float32x4_t>{
    const static int V_WIDTH=4;
};
#include "arch/arm/Mask.NEON.h"
#include "arch/arm/Number.NEON.h"
#include "arch/arm/Discrete.NEON.h"
#endif

// Define Architecture Specific Helper Macros and Miscellanious 

#if defined(ENABLE_SSE_INSTRUCTION_SET) || defined(ENABLE_AVX_INSTRUCTION_SET)
#define KERNEL_MEM_ALIGN __declspec(align(32))
#else
#define KERNEL_MEM_ALIGN
#endif


// Define this to switch to non-augmented material formulas
#define USE_NONMIXED_FORMULAS

struct COROTATED_TAG;
struct NEOHOOKEAN_TAG;
struct BIPHASIC_TAG;

namespace nm_biphasic {
    extern float *biphasic_threshold;
    extern float *biphasic_factor;
    //BUILD_CONSTANT( biphasic_threshold, 1.07f );
    //BUILD_CONSTANT( biphasic_factor, 1e4);
}

template <class T, int width>
struct widetype{
   typedef T type[width];
};

template <class T>
struct widetype<T,1>{
   typedef T type;
};

#define WIDETYPE(TYPE,WIDTH) widetype<TYPE, WIDTH>::type 
#define VWIDTH(TYPE) VTYPE_POLICY<TYPE>::V_WIDTH

#define INSTANCE_KERNEL_SCALAR(KERNEL,DATA_WIDTH) template void KERNEL<float,WIDETYPE(float,DATA_WIDTH),WIDETYPE(int,DATA_WIDTH)>( INSTANCE_KERNEL_ ## KERNEL(DATA_WIDTH) );
#define INSTANCE_KERNEL_SIMD_FLOAT(KERNEL,DATA_WIDTH) template void KERNEL<float,WIDETYPE(float,DATA_WIDTH),WIDETYPE(int,DATA_WIDTH)>( INSTANCE_KERNEL_ ## KERNEL(DATA_WIDTH) );
#ifdef ENABLE_SSE_INSTRUCTION_SET
#define INSTANCE_KERNEL_SIMD_SSE(KERNEL,DATA_WIDTH) template void KERNEL<__m128,WIDETYPE(float,DATA_WIDTH),WIDETYPE(int,DATA_WIDTH)>( INSTANCE_KERNEL_ ## KERNEL(DATA_WIDTH) ); 
#else
#define INSTANCE_KERNEL_SIMD_SSE(KERNEL,DATA_WIDTH)
#endif
#ifdef ENABLE_AVX_INSTRUCTION_SET
#define INSTANCE_KERNEL_SIMD_AVX(KERNEL,DATA_WIDTH) template void KERNEL<__m256,WIDETYPE(float,DATA_WIDTH),WIDETYPE(int,DATA_WIDTH)>( INSTANCE_KERNEL_ ## KERNEL(DATA_WIDTH) );
#else
#define INSTANCE_KERNEL_SIMD_AVX(KERNEL,DATA_WIDTH)
#endif
#ifdef ENABLE_MIC_INSTRUCTION_SET
#define INSTANCE_KERNEL_SIMD_MIC(KERNEL,DATA_WIDTH) template void KERNEL<__m512,WIDETYPE(float,DATA_WIDTH),WIDETYPE(int,DATA_WIDTH)>( INSTANCE_KERNEL_ ## KERNEL(DATA_WIDTH) );
#else
#define INSTANCE_KERNEL_SIMD_MIC(KERNEL,DATA_WIDTH)
#endif
#ifdef ENABLE_NEON_INSTRUCTION_SET
#define INSTANCE_KERNEL_SIMD_NEON(KERNEL,DATA_WIDTH) template void KERNEL<float32x4_t,WIDETYPE(float,DATA_WIDTH),WIDETYPE(int,DATA_WIDTH)>( INSTANCE_KERNEL_ ## KERNEL(DATA_WIDTH) );
#else
#define INSTANCE_KERNEL_SIMD_NEON(KERNEL,DATA_WIDTH)
#endif
#define INSTANCE_KERNEL( KERNEL )               \
    INSTANCE_KERNEL_SCALAR( KERNEL, 1)          \
    INSTANCE_KERNEL_SIMD_FLOAT( KERNEL, 4)      \
    INSTANCE_KERNEL_SIMD_FLOAT( KERNEL, 8)      \
    INSTANCE_KERNEL_SIMD_FLOAT( KERNEL, 16)     \
    INSTANCE_KERNEL_SIMD_SSE( KERNEL, 4)        \
    INSTANCE_KERNEL_SIMD_SSE( KERNEL, 8)        \
    INSTANCE_KERNEL_SIMD_SSE( KERNEL, 16)       \
    INSTANCE_KERNEL_SIMD_AVX( KERNEL, 8)        \
    INSTANCE_KERNEL_SIMD_AVX( KERNEL, 16)       \
    INSTANCE_KERNEL_SIMD_MIC( KERNEL, 16)       \
    INSTANCE_KERNEL_SIMD_NEON( KERNEL, 4)       \
    INSTANCE_KERNEL_SIMD_NEON( KERNEL, 8)       \
    INSTANCE_KERNEL_SIMD_NEON( KERNEL, 16)      


#define INSTANCE_KERNEL_MATERIAL_SCALAR(KERNEL,MATERIAL,DATA_WIDTH) template void KERNEL<MATERIAL,float,WIDETYPE(float,DATA_WIDTH),WIDETYPE(int,DATA_WIDTH)>::Run( INSTANCE_KERNEL_ ## KERNEL(DATA_WIDTH) );
#define INSTANCE_KERNEL_MATERIAL_SIMD_FLOAT(KERNEL,MATERIAL,DATA_WIDTH) template void KERNEL<MATERIAL,float,WIDETYPE(float,DATA_WIDTH),WIDETYPE(int,DATA_WIDTH)>::Run( INSTANCE_KERNEL_ ## KERNEL(DATA_WIDTH) );
#ifdef ENABLE_SSE_INSTRUCTION_SET
#define INSTANCE_KERNEL_MATERIAL_SIMD_SSE(KERNEL,MATERIAL,DATA_WIDTH) template void KERNEL<MATERIAL,__m128,WIDETYPE(float,DATA_WIDTH),WIDETYPE(int,DATA_WIDTH)>::Run( INSTANCE_KERNEL_ ## KERNEL(DATA_WIDTH) );
#else
#define INSTANCE_KERNEL_MATERIAL_SIMD_SSE(KERNEL,MATERIAL,DATA_WIDTH)
#endif
#ifdef ENABLE_AVX_INSTRUCTION_SET
#define INSTANCE_KERNEL_MATERIAL_SIMD_AVX(KERNEL,MATERIAL,DATA_WIDTH) template void KERNEL<MATERIAL,__m256,WIDETYPE(float,DATA_WIDTH),WIDETYPE(int,DATA_WIDTH)>::Run( INSTANCE_KERNEL_ ## KERNEL(DATA_WIDTH) );
#else
#define INSTANCE_KERNEL_MATERIAL_SIMD_AVX(KERNEL,MATERIAL,DATA_WIDTH)
#endif
#ifdef ENABLE_MIC_INSTRUCTION_SET
#define INSTANCE_KERNEL_MATERIAL_SIMD_MIC(KERNEL,MATERIAL,DATA_WIDTH) template void KERNEL<MATERIAL,__m512,WIDETYPE(float,DATA_WIDTH),WIDETYPE(int,DATA_WIDTH)>::Run( INSTANCE_KERNEL_ ## KERNEL(DATA_WIDTH) );
#else
#define INSTANCE_KERNEL_MATERIAL_SIMD_MIC(KERNEL,MATERIAL,DATA_WIDTH)
#endif
#ifdef ENABLE_NEON_INSTRUCTION_SET
#define INSTANCE_KERNEL_MATERIAL_SIMD_NEON(KERNEL,MATERIAL,DATA_WIDTH) template void KERNEL<MATERIAL,float32x4_t,WIDETYPE(float,DATA_WIDTH),WIDETYPE(int,DATA_WIDTH)>::Run( INSTANCE_KERNEL_ ## KERNEL(DATA_WIDTH) );
#else
#define INSTANCE_KERNEL_MATERIAL_SIMD_NEON(KERNEL,MATERIAL,DATA_WIDTH)
#endif
#define INSTANCE_KERNEL_MATERIAL( KERNEL, MATERIAL)                \
    INSTANCE_KERNEL_MATERIAL_SCALAR( KERNEL, MATERIAL, 1)          \
    INSTANCE_KERNEL_MATERIAL_SIMD_FLOAT( KERNEL, MATERIAL, 4)      \
    INSTANCE_KERNEL_MATERIAL_SIMD_FLOAT( KERNEL, MATERIAL, 8)      \
    INSTANCE_KERNEL_MATERIAL_SIMD_FLOAT( KERNEL, MATERIAL, 16)     \
    INSTANCE_KERNEL_MATERIAL_SIMD_SSE( KERNEL, MATERIAL, 4)        \
    INSTANCE_KERNEL_MATERIAL_SIMD_SSE( KERNEL, MATERIAL, 8)        \
    INSTANCE_KERNEL_MATERIAL_SIMD_SSE( KERNEL, MATERIAL, 16)       \
    INSTANCE_KERNEL_MATERIAL_SIMD_AVX( KERNEL, MATERIAL, 8)        \
    INSTANCE_KERNEL_MATERIAL_SIMD_AVX( KERNEL, MATERIAL, 16)       \
    INSTANCE_KERNEL_MATERIAL_SIMD_MIC( KERNEL, MATERIAL, 16)       \
    INSTANCE_KERNEL_MATERIAL_SIMD_NEON( KERNEL, MATERIAL, 4)       \
    INSTANCE_KERNEL_MATERIAL_SIMD_NEON( KERNEL, MATERIAL, 8)       \
    INSTANCE_KERNEL_MATERIAL_SIMD_NEON( KERNEL, MATERIAL, 16)      

#endif

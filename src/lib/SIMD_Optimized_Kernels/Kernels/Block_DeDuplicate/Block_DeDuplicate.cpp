//#####################################################################
//  Copyright (c) 2011-2013 Nathan Mitchell, Eftychios Sifakis.
//  This file is covered by the FreeBSD license. Please refer to the 
//  license.txt file for more information.
//#####################################################################


#ifndef SUBROUTINE_Block_DeDuplicate
#include "KernelCommon.h"
#else
namespace {
#endif

#ifdef __INTEL_COMPILER
#pragma warning( disable : 592 )
#endif

#include "Block_DeDuplicate.h"
#include <string.h>


template<class T_RAW>
#ifdef SUBROUTINE_Block_DeDuplicate
inline
#endif
void Block_DeDuplicate(const float u_duplicated[3][8][8], float u_compact[3][27])
{
    typedef float T3333[3][3][3][3];

    T3333& u_compact_3333 = *(reinterpret_cast<T3333*>(u_compact));

    memset( u_compact_3333, 0, sizeof(float)*3*3*3*3 );

    for( int v=0; v<3; v++)
        for( int i_base=0; i_base<2; i_base++)
            for( int j_base=0; j_base<2; j_base++)
                for( int k_base=0; k_base<2; k_base++)
                    for( int i=0; i<2; i++)
                        for( int j=0; j<2; j++)
                            for( int k=0; k<2; k++)
                                {
                                    u_compact_3333[v][i+i_base][j+j_base][k+k_base] += 
                                        u_duplicated[v][i*4+j*2+k][i_base*4+j_base*2+k_base];
                                }

}


#ifdef ENABLE_AVX_INSTRUCTION_SET
template<>
#ifdef SUBROUTINE_Block_DeDuplicate
inline
#endif
void Block_DeDuplicate<__m256>(const float u_duplicated[3][8][8],float u_compact[3][27])
{
    __m256 rzero;
    rzero=_mm256_xor_ps(rzero,rzero);

    for(int v=0;v<3;v++)
    {
        __m256 rtmp1,rtmp2;

        __m256 rA=_mm256_loadu_ps(u_duplicated[v][0]);
        __m256 rB=_mm256_loadu_ps(u_duplicated[v][1]);
        __m256 rC=_mm256_loadu_ps(u_duplicated[v][2]);
        __m256 rD=_mm256_loadu_ps(u_duplicated[v][3]);
        __m256 rE=_mm256_loadu_ps(u_duplicated[v][4]);
        __m256 rF=_mm256_loadu_ps(u_duplicated[v][5]);
        __m256 rG=_mm256_loadu_ps(u_duplicated[v][6]);
        __m256 rH=_mm256_loadu_ps(u_duplicated[v][7]);

        __m256 r1=_mm256_permute2f128_ps(rA,rC,0x20);
        __m256 r2=_mm256_permute2f128_ps(rB,rD,0x20);
        __m256 r3=_mm256_permute2f128_ps(rE,rG,0x20);
        __m256 r4=_mm256_permute2f128_ps(rF,rH,0x20);
        __m256 r3a=_mm256_permute2f128_ps(rA,rC,0x31);
        __m256 r4a=_mm256_permute2f128_ps(rB,rD,0x31);
        __m256 r5=_mm256_permute2f128_ps(rE,rG,0x31);
        __m256 r6=_mm256_permute2f128_ps(rF,rH,0x31);
        r3=_mm256_add_ps(r3,r3a);
        r4=_mm256_add_ps(r4,r4a);

        // Aggregate values 0-8

        rtmp1=_mm256_permute2f128_ps(r1,r2,0x30);
        rtmp2=_mm256_permute2f128_ps(r1,r2,0x21);
        rtmp2=_mm256_permute_ps(rtmp2,0x4e);
        rtmp1=_mm256_add_ps(rtmp1,rtmp2);
        r1=_mm256_blend_ps(r1,rzero,0x0c);
        r1=_mm256_blend_ps(r1,rtmp1,0x30);
        r2=_mm256_blend_ps(r2,rzero,0x30);
        r2=_mm256_blend_ps(r2,rtmp1,0x0c);

        _mm256_storeu_ps(&u_compact[v][0],r1);
        rtmp1=_mm256_loadu_ps(&u_compact[v][1]);
        rtmp1=_mm256_blend_ps(rtmp1,rzero,0x80);
        r2=_mm256_add_ps(r2,rtmp1);
        _mm256_storeu_ps(&u_compact[v][1],r2);

        // Aggregate values 9-17

        rtmp1=_mm256_permute2f128_ps(r3,r4,0x30);
        rtmp2=_mm256_permute2f128_ps(r3,r4,0x21);
        rtmp2=_mm256_permute_ps(rtmp2,0x4e);
        rtmp1=_mm256_add_ps(rtmp1,rtmp2);
        r3=_mm256_blend_ps(r3,rzero,0x0c);
        r3=_mm256_blend_ps(r3,rtmp1,0x30);
        r4=_mm256_blend_ps(r4,rzero,0x30);
        r4=_mm256_blend_ps(r4,rtmp1,0x0c);

        _mm256_storeu_ps(&u_compact[v][9],r3);
        rtmp1=_mm256_loadu_ps(&u_compact[v][10]);
        rtmp1=_mm256_blend_ps(rtmp1,rzero,0x80);
        r4=_mm256_add_ps(r4,rtmp1);
        _mm256_storeu_ps(&u_compact[v][10],r4);

        // Aggregate values 18-26

        rtmp1=_mm256_permute2f128_ps(r5,r6,0x30);
        rtmp2=_mm256_permute2f128_ps(r5,r6,0x21);
        rtmp2=_mm256_permute_ps(rtmp2,0x4e);
        rtmp1=_mm256_add_ps(rtmp1,rtmp2);
        r5=_mm256_blend_ps(r5,rzero,0x0c);
        r5=_mm256_blend_ps(r5,rtmp1,0x30);
        r6=_mm256_blend_ps(r6,rzero,0x30);
        r6=_mm256_blend_ps(r6,rtmp1,0x0c);

        _mm256_storeu_ps(&u_compact[v][18],r5);
        rtmp1=_mm256_loadu_ps(&u_compact[v][19]);
        rtmp1=_mm256_blend_ps(rtmp1,rzero,0x80);
        r6=_mm256_add_ps(r6,rtmp1);
        _mm256_storeu_ps(&u_compact[v][19],r6);
    }
}
#endif

template<class T_RAW>
#ifdef SUBROUTINE_Block_DeDuplicate
inline
#endif
void Block_DeDuplicate(const float u_duplicated[3][8][16], float u_compact[3][27])
{
    typedef float T3333[3][3][3][3];

    T3333& u_compact_3333 = *(reinterpret_cast<T3333*>(u_compact));

    memset( u_compact_3333, 0, sizeof(float)*3*3*3*3 );

    for( int v=0; v<3; v++)
        for( int i_base=0; i_base<2; i_base++)
            for( int j_base=0; j_base<2; j_base++)
                for( int k_base=0; k_base<2; k_base++)
                    for( int i=0; i<2; i++)
                        for( int j=0; j<2; j++)
                            for( int k=0; k<2; k++)
                                {
                                    u_compact_3333[v][i+i_base][j+j_base][k+k_base] += 
                                        u_duplicated[v][i*4+j*2+k][i_base*4+j_base*2+k_base];
                                }

}


#ifdef ENABLE_AVX_INSTRUCTION_SET
template<>
#ifdef SUBROUTINE_Block_DeDuplicate
inline
#endif
void Block_DeDuplicate<__m256>(const float u_duplicated[3][8][16],float u_compact[3][27])
{
    __m256 rzero;
    rzero=_mm256_xor_ps(rzero,rzero);

    for(int v=0;v<3;v++)
    {
        __m256 rtmp1,rtmp2;

        __m256 rA=_mm256_loadu_ps(u_duplicated[v][0]);
        __m256 rB=_mm256_loadu_ps(u_duplicated[v][1]);
        __m256 rC=_mm256_loadu_ps(u_duplicated[v][2]);
        __m256 rD=_mm256_loadu_ps(u_duplicated[v][3]);
        __m256 rE=_mm256_loadu_ps(u_duplicated[v][4]);
        __m256 rF=_mm256_loadu_ps(u_duplicated[v][5]);
        __m256 rG=_mm256_loadu_ps(u_duplicated[v][6]);
        __m256 rH=_mm256_loadu_ps(u_duplicated[v][7]);

        __m256 r1=_mm256_permute2f128_ps(rA,rC,0x20);
        __m256 r2=_mm256_permute2f128_ps(rB,rD,0x20);
        __m256 r3=_mm256_permute2f128_ps(rE,rG,0x20);
        __m256 r4=_mm256_permute2f128_ps(rF,rH,0x20);
        __m256 r3a=_mm256_permute2f128_ps(rA,rC,0x31);
        __m256 r4a=_mm256_permute2f128_ps(rB,rD,0x31);
        __m256 r5=_mm256_permute2f128_ps(rE,rG,0x31);
        __m256 r6=_mm256_permute2f128_ps(rF,rH,0x31);
        r3=_mm256_add_ps(r3,r3a);
        r4=_mm256_add_ps(r4,r4a);

        // Aggregate values 0-8

        rtmp1=_mm256_permute2f128_ps(r1,r2,0x30);
        rtmp2=_mm256_permute2f128_ps(r1,r2,0x21);
        rtmp2=_mm256_permute_ps(rtmp2,0x4e);
        rtmp1=_mm256_add_ps(rtmp1,rtmp2);
        r1=_mm256_blend_ps(r1,rzero,0x0c);
        r1=_mm256_blend_ps(r1,rtmp1,0x30);
        r2=_mm256_blend_ps(r2,rzero,0x30);
        r2=_mm256_blend_ps(r2,rtmp1,0x0c);

        _mm256_storeu_ps(&u_compact[v][0],r1);
        rtmp1=_mm256_loadu_ps(&u_compact[v][1]);
        rtmp1=_mm256_blend_ps(rtmp1,rzero,0x80);
        r2=_mm256_add_ps(r2,rtmp1);
        _mm256_storeu_ps(&u_compact[v][1],r2);

        // Aggregate values 9-17

        rtmp1=_mm256_permute2f128_ps(r3,r4,0x30);
        rtmp2=_mm256_permute2f128_ps(r3,r4,0x21);
        rtmp2=_mm256_permute_ps(rtmp2,0x4e);
        rtmp1=_mm256_add_ps(rtmp1,rtmp2);
        r3=_mm256_blend_ps(r3,rzero,0x0c);
        r3=_mm256_blend_ps(r3,rtmp1,0x30);
        r4=_mm256_blend_ps(r4,rzero,0x30);
        r4=_mm256_blend_ps(r4,rtmp1,0x0c);

        _mm256_storeu_ps(&u_compact[v][9],r3);
        rtmp1=_mm256_loadu_ps(&u_compact[v][10]);
        rtmp1=_mm256_blend_ps(rtmp1,rzero,0x80);
        r4=_mm256_add_ps(r4,rtmp1);
        _mm256_storeu_ps(&u_compact[v][10],r4);

        // Aggregate values 18-26

        rtmp1=_mm256_permute2f128_ps(r5,r6,0x30);
        rtmp2=_mm256_permute2f128_ps(r5,r6,0x21);
        rtmp2=_mm256_permute_ps(rtmp2,0x4e);
        rtmp1=_mm256_add_ps(rtmp1,rtmp2);
        r5=_mm256_blend_ps(r5,rzero,0x0c);
        r5=_mm256_blend_ps(r5,rtmp1,0x30);
        r6=_mm256_blend_ps(r6,rzero,0x30);
        r6=_mm256_blend_ps(r6,rtmp1,0x0c);

        _mm256_storeu_ps(&u_compact[v][18],r5);
        rtmp1=_mm256_loadu_ps(&u_compact[v][19]);
        rtmp1=_mm256_blend_ps(rtmp1,rzero,0x80);
        r6=_mm256_add_ps(r6,rtmp1);
        _mm256_storeu_ps(&u_compact[v][19],r6);
    }
}
#endif


#ifdef ENABLE_MIC_INSTRUCTION_SET
namespace {
    const int Offsets[16] __attribute__((aligned(64))) {
        0+0+0+0,
            0+0+0+1,
            0+0+3+0,
            0+0+3+1,
            0+9+0+0,
            0+9+0+1,
            0+9+3+0,
            0+9+3+1,
            81+0+0+0,
            81+0+0+1,
            81+0+3+0,
            81+0+3+1,
            81+9+0+0,
            81+9+0+1,
            81+9+3+0,
            81+9+3+1
            };
}
template<>
#ifdef SUBROUTINE_Block_Duplicate
inline
#endif
void Block_DeDuplicate<__m512>(const float u_duplicated[3][8][16], float u_compact[2][3][27])
{
    // Replace this?
    memset( u_compact, 0, sizeof(float)*2*3*3*3*3 );

    __m512 data0, data1, data2, data3, data4, data5, data6, data7;
    __m512i offsets=_mm512_load_epi32(Offsets);

    for(int coord=0;coord<3;coord++){
        data0=_mm512_load_ps(&u_duplicated[coord][0][0]);
        data1=_mm512_load_ps(&u_duplicated[coord][1][0]);
        data2=_mm512_load_ps(&u_duplicated[coord][2][0]);
        data3=_mm512_load_ps(&u_duplicated[coord][3][0]);
        data4=_mm512_load_ps(&u_duplicated[coord][4][0]);
        data5=_mm512_load_ps(&u_duplicated[coord][5][0]);
        data6=_mm512_load_ps(&u_duplicated[coord][6][0]);
        data7=_mm512_load_ps(&u_duplicated[coord][7][0]);

        // First store 0
        _mm512_i32scatter_ps((void*)&u_compact[0][coord][0],offsets,data0,4);

        // Load, Add, and Store 1
        data0 =_mm512_i32gather_ps(offsets,(void*)&u_compact[0][coord][1],4);
        data0 = _mm512_add_ps( data0, data1 );
        _mm512_i32scatter_ps((void*)&u_compact[0][coord][1],offsets,data0,4);


        // Load, Add, and Store 2
        data0 =_mm512_i32gather_ps(offsets,(void*)&u_compact[0][coord][3],4);
        data0 = _mm512_add_ps( data0, data2 );
        _mm512_i32scatter_ps((void*)&u_compact[0][coord][3],offsets,data0,4);


        // Load, Add, and Store 3
        data0 =_mm512_i32gather_ps(offsets,(void*)&u_compact[0][coord][4],4);
        data0 = _mm512_add_ps( data0, data3 );
        _mm512_i32scatter_ps((void*)&u_compact[0][coord][4],offsets,data0,4);


        // Load, Add, and Store 4
        data0 =_mm512_i32gather_ps(offsets,(void*)&u_compact[0][coord][9],4);
        data0 = _mm512_add_ps( data0, data4 );
        _mm512_i32scatter_ps((void*)&u_compact[0][coord][9],offsets,data0,4);


        // Load, Add, and Store 5
        data0 =_mm512_i32gather_ps(offsets,(void*)&u_compact[0][coord][10],4);
        data0 = _mm512_add_ps( data0, data5 );
        _mm512_i32scatter_ps((void*)&u_compact[0][coord][10],offsets,data0,4);


        // Load, Add, and Store 6
        data0 =_mm512_i32gather_ps(offsets,(void*)&u_compact[0][coord][12],4);
        data0 = _mm512_add_ps( data0, data6 );
        _mm512_i32scatter_ps((void*)&u_compact[0][coord][12],offsets,data0,4);


        // Load, Add, and Store 7 
        data0 =_mm512_i32gather_ps(offsets,(void*)&u_compact[0][coord][13],4);
        data0 = _mm512_add_ps( data0, data7 );
        _mm512_i32scatter_ps((void*)&u_compact[0][coord][13],offsets,data0,4);

        __m512 data2=_mm512_i32gather_ps(offsets,(void*)&u_compact[0][coord][3],4);
        __m512 data3=_mm512_i32gather_ps(offsets,(void*)&u_compact[0][coord][4],4);
        __m512 data4=_mm512_i32gather_ps(offsets,(void*)&u_compact[0][coord][9],4);
        __m512 data5=_mm512_i32gather_ps(offsets,(void*)&u_compact[0][coord][10],4);
        __m512 data6=_mm512_i32gather_ps(offsets,(void*)&u_compact[0][coord][12],4);
        __m512 data7=_mm512_i32gather_ps(offsets,(void*)&u_compact[0][coord][13],4);

    }
}
#endif



#ifdef __INTEL_COMPILER
#pragma warning( default : 592 )
#endif

#ifdef SUBROUTINE_Block_DeDuplicate
}
#else
template void Block_DeDuplicate<float>(const float u_duplicated[3][8][16], float u_compact[3][27]);
#ifdef ENABLE_SSE_INSTRUCTION_SET
template void Block_DeDuplicate<__m128>(const float u_duplicated[3][8][16], float u_compact[3][27]);
#endif
#ifdef ENABLE_MIC_INSTRUCTION_SET
template void Block_DeDuplicate<__m512>(const float u_duplicated[3][8][16], float u_compact[3][27]);
template void Block_DeDuplicate<__m512>(const float u_duplicated[3][8][16], float u_compact[2][3][27]);
#endif
template void Block_DeDuplicate<float>(const float u_duplicated[3][8][8], float u_compact[3][27]);
#ifdef ENABLE_SSE_INSTRUCTION_SET
template void Block_DeDuplicate<__m128>(const float u_duplicated[3][8][8], float u_compact[3][27]);
#endif
#endif

//#####################################################################
//  Copyright (c) 2011-2013 Nathan Mitchell, Eftychios Sifakis.
//  This file is covered by the FreeBSD license. Please refer to the 
//  license.txt file for more information.
//#####################################################################


#ifndef SUBROUTINE_Block_Duplicate
#include "KernelCommon.h"
#else
namespace {
#endif

#include "Block_Duplicate.h"

template<class T_RAW>
#ifdef SUBROUTINE_Block_Duplicate
inline
#endif
void Block_Duplicate(const float u_compact[3][27], float u_duplicated[3][8][8])
{
    typedef float T3333[3][3][3][3];

    const T3333& u_compact_3333 = *(reinterpret_cast<const T3333*>(u_compact));

    for( int v=0; v<3; v++)
        for( int i_base=0; i_base<2; i_base++)
            for( int j_base=0; j_base<2; j_base++)
                for( int k_base=0; k_base<2; k_base++)
                    for( int i=0; i<2; i++)
                        for( int j=0; j<2; j++)
                            for( int k=0; k<2; k++)
                                {
                                    u_duplicated[v][i*4+j*2+k][i_base*4+j_base*2+k_base] =
                                        u_compact_3333[v][i+i_base][j+j_base][k+k_base];
                                }

}


#ifdef ENABLE_AVX_INSTRUCTION_SET
template<>
#ifdef SUBROUTINE_Block_Duplicate
inline
#endif
void Block_Duplicate<__m256>(const float u_compact[3][27], float u_duplicated[3][8][8])
{
    for(int v=0;v<3;v++)
    {
        __m256 rtmp;

        __m256 r1=_mm256_loadu_ps(&u_compact[v][0]);
        __m256 r2=_mm256_loadu_ps(&u_compact[v][1]);
        rtmp=_mm256_permute2f128_ps(r1,r2,0x21);
        rtmp=_mm256_permute_ps(rtmp,0x4e);
        r1=_mm256_blend_ps(r1,r2,0x0c);
        r2=_mm256_blend_ps(r2,rtmp,0x0c);
        r2=_mm256_blend_ps(r2,r1,0x30);
        r1=_mm256_blend_ps(r1,rtmp,0x30);

        __m256 r3=_mm256_loadu_ps(&u_compact[v][9]);
        __m256 r4=_mm256_loadu_ps(&u_compact[v][10]);
        rtmp=_mm256_permute2f128_ps(r3,r4,0x21);
        rtmp=_mm256_permute_ps(rtmp,0x4e);
        r3=_mm256_blend_ps(r3,r4,0x0c);
        r4=_mm256_blend_ps(r4,rtmp,0x0c);
        r4=_mm256_blend_ps(r4,r3,0x30);
        r3=_mm256_blend_ps(r3,rtmp,0x30);

        __m256 r5=_mm256_loadu_ps(&u_compact[v][18]);
        __m256 r6=_mm256_loadu_ps(&u_compact[v][19]);
        rtmp=_mm256_permute2f128_ps(r5,r6,0x21);
        rtmp=_mm256_permute_ps(rtmp,0x4e);
        r5=_mm256_blend_ps(r5,r6,0x0c);
        r6=_mm256_blend_ps(r6,rtmp,0x0c);
        r6=_mm256_blend_ps(r6,r5,0x30);
        r5=_mm256_blend_ps(r5,rtmp,0x30);

        __m256 rA=_mm256_permute2f128_ps(r1,r3,0x20);
        __m256 rB=_mm256_permute2f128_ps(r2,r4,0x20);
        __m256 rC=_mm256_permute2f128_ps(r1,r3,0x31);
        __m256 rD=_mm256_permute2f128_ps(r2,r4,0x31);

        __m256 rE=_mm256_permute2f128_ps(r3,r5,0x20);
        __m256 rF=_mm256_permute2f128_ps(r4,r6,0x20);
        __m256 rG=_mm256_permute2f128_ps(r3,r5,0x31);
        __m256 rH=_mm256_permute2f128_ps(r4,r6,0x31);

        _mm256_storeu_ps(u_duplicated[v][0],rA);
        _mm256_storeu_ps(u_duplicated[v][1],rB);
        _mm256_storeu_ps(u_duplicated[v][2],rC);
        _mm256_storeu_ps(u_duplicated[v][3],rD);
        _mm256_storeu_ps(u_duplicated[v][4],rE);
        _mm256_storeu_ps(u_duplicated[v][5],rF);
        _mm256_storeu_ps(u_duplicated[v][6],rG);
        _mm256_storeu_ps(u_duplicated[v][7],rH);
    }
}
#endif



template<class T_RAW>
#ifdef SUBROUTINE_Block_Duplicate
inline
#endif
void Block_Duplicate(const float u_compact[3][27], float u_duplicated[3][8][16])
{
    typedef float T3333[3][3][3][3];

    const T3333& u_compact_3333 = *(reinterpret_cast<const T3333*>(u_compact));

    for( int v=0; v<3; v++)
        for( int i_base=0; i_base<2; i_base++)
            for( int j_base=0; j_base<2; j_base++)
                for( int k_base=0; k_base<2; k_base++)
                    for( int i=0; i<2; i++)
                        for( int j=0; j<2; j++)
                            for( int k=0; k<2; k++)
                                {
                                    u_duplicated[v][i*4+j*2+k][i_base*4+j_base*2+k_base] =
                                        u_compact_3333[v][i+i_base][j+j_base][k+k_base];
                                }

}


#ifdef ENABLE_AVX_INSTRUCTION_SET
template<>
#ifdef SUBROUTINE_Block_Duplicate
inline
#endif
void Block_Duplicate<__m256>(const float u_compact[3][27], float u_duplicated[3][8][16])
{
    for(int v=0;v<3;v++)
    {
        __m256 rtmp;

        __m256 r1=_mm256_loadu_ps(&u_compact[v][0]);
        __m256 r2=_mm256_loadu_ps(&u_compact[v][1]);
        rtmp=_mm256_permute2f128_ps(r1,r2,0x21);
        rtmp=_mm256_permute_ps(rtmp,0x4e);
        r1=_mm256_blend_ps(r1,r2,0x0c);
        r2=_mm256_blend_ps(r2,rtmp,0x0c);
        r2=_mm256_blend_ps(r2,r1,0x30);
        r1=_mm256_blend_ps(r1,rtmp,0x30);

        __m256 r3=_mm256_loadu_ps(&u_compact[v][9]);
        __m256 r4=_mm256_loadu_ps(&u_compact[v][10]);
        rtmp=_mm256_permute2f128_ps(r3,r4,0x21);
        rtmp=_mm256_permute_ps(rtmp,0x4e);
        r3=_mm256_blend_ps(r3,r4,0x0c);
        r4=_mm256_blend_ps(r4,rtmp,0x0c);
        r4=_mm256_blend_ps(r4,r3,0x30);
        r3=_mm256_blend_ps(r3,rtmp,0x30);

        __m256 r5=_mm256_loadu_ps(&u_compact[v][18]);
        __m256 r6=_mm256_loadu_ps(&u_compact[v][19]);
        rtmp=_mm256_permute2f128_ps(r5,r6,0x21);
        rtmp=_mm256_permute_ps(rtmp,0x4e);
        r5=_mm256_blend_ps(r5,r6,0x0c);
        r6=_mm256_blend_ps(r6,rtmp,0x0c);
        r6=_mm256_blend_ps(r6,r5,0x30);
        r5=_mm256_blend_ps(r5,rtmp,0x30);

        __m256 rA=_mm256_permute2f128_ps(r1,r3,0x20);
        __m256 rB=_mm256_permute2f128_ps(r2,r4,0x20);
        __m256 rC=_mm256_permute2f128_ps(r1,r3,0x31);
        __m256 rD=_mm256_permute2f128_ps(r2,r4,0x31);

        __m256 rE=_mm256_permute2f128_ps(r3,r5,0x20);
        __m256 rF=_mm256_permute2f128_ps(r4,r6,0x20);
        __m256 rG=_mm256_permute2f128_ps(r3,r5,0x31);
        __m256 rH=_mm256_permute2f128_ps(r4,r6,0x31);

        _mm256_storeu_ps(u_duplicated[v][0],rA);
        _mm256_storeu_ps(u_duplicated[v][1],rB);
        _mm256_storeu_ps(u_duplicated[v][2],rC);
        _mm256_storeu_ps(u_duplicated[v][3],rD);
        _mm256_storeu_ps(u_duplicated[v][4],rE);
        _mm256_storeu_ps(u_duplicated[v][5],rF);
        _mm256_storeu_ps(u_duplicated[v][6],rG);
        _mm256_storeu_ps(u_duplicated[v][7],rH);
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
void Block_Duplicate<__m512>(const float u_compact[2][3][27], float u_duplicated[3][8][16])
{

    __m512i offsets=_mm512_load_epi32(Offsets);

    for(int coord=0;coord<3;coord++){
        __m512 data0=_mm512_i32gather_ps(offsets,(void*)&u_compact[0][coord][0],4);
        __m512 data1=_mm512_i32gather_ps(offsets,(void*)&u_compact[0][coord][1],4);
        __m512 data2=_mm512_i32gather_ps(offsets,(void*)&u_compact[0][coord][3],4);
        __m512 data3=_mm512_i32gather_ps(offsets,(void*)&u_compact[0][coord][4],4);
        __m512 data4=_mm512_i32gather_ps(offsets,(void*)&u_compact[0][coord][9],4);
        __m512 data5=_mm512_i32gather_ps(offsets,(void*)&u_compact[0][coord][10],4);
        __m512 data6=_mm512_i32gather_ps(offsets,(void*)&u_compact[0][coord][12],4);
        __m512 data7=_mm512_i32gather_ps(offsets,(void*)&u_compact[0][coord][13],4);

        _mm512_store_ps(&u_duplicated[coord][0][0],data0);
        _mm512_store_ps(&u_duplicated[coord][1][0],data1);
        _mm512_store_ps(&u_duplicated[coord][2][0],data2);
        _mm512_store_ps(&u_duplicated[coord][3][0],data3);
        _mm512_store_ps(&u_duplicated[coord][4][0],data4);
        _mm512_store_ps(&u_duplicated[coord][5][0],data5);
        _mm512_store_ps(&u_duplicated[coord][6][0],data6);
        _mm512_store_ps(&u_duplicated[coord][7][0],data7);
    }
}
#endif




#ifdef SUBROUTINE_Block_Duplicate
}
#else
template void Block_Duplicate<float>(const float u_compact[3][27], float u_duplicated[3][8][16]);
#ifdef ENABLE_SSE_INSTRUCTION_SET
template void Block_Duplicate<__m128>(const float u_compact[3][27], float u_duplicated[3][8][16]);
#endif
#ifdef ENABLE_MIC_INSTRUCTION_SET
template void Block_Duplicate<__m512>(const float u_compact[3][27], float u_duplicated[3][8][16]);
template void Block_Duplicate<__m512>(const float u_compact[2][3][27], float u_duplicated[3][8][16]);
#endif
template void Block_Duplicate<float>(const float u_compact[3][27], float u_duplicated[3][8][8]);
#ifdef ENABLE_SSE_INSTRUCTION_SET
template void Block_Duplicate<__m128>(const float u_compact[3][27], float u_duplicated[3][8][8]);
#endif
#endif

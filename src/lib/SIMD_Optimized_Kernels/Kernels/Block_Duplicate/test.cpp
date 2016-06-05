#include <immintrin.h>
#include <iostream>

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


int main(int argc,char *argv[])
{
    float u_compact[2][3][27];
    float u_duplicated[3][8][16] __attribute__((aligned(64)));
    float u_reference[3][8][16];

    __m512i offsets=_mm512_load_epi32(Offsets);

    int count=0;
    for(int i=0;i<2;i++)
        for(int j=0;j<3;j++)
            for(int k=0;k<27;k++)
                u_compact[i][j][k]=++count;

    typedef float (&array_type1)[2][3][3][3][3];
    typedef float (&array_type2)[3][2][2][2][2][2][2][2];

    array_type1 u_compact_as_type1=reinterpret_cast<array_type1>(u_compact[0][0][0]);
    array_type2 u_reference_as_type2=reinterpret_cast<array_type2>(u_reference[0][0][0]);
    array_type2 u_duplicated_as_type2=reinterpret_cast<array_type2>(u_duplicated[0][0][0]);

    for(int coord=0;coord<3;coord++){
        int vindex=0;
        for(int ivertex=0;ivertex<=1;ivertex++)
        for(int jvertex=0;jvertex<=1;jvertex++)
        for(int kvertex=0;kvertex<=1;kvertex++){
            int cindex=0;
            for(int block=0;block<2;block++)
            for(int icell=0;icell<=1;icell++)
            for(int jcell=0;jcell<=1;jcell++)
            for(int kcell=0;kcell<=1;kcell++)
                u_reference[coord][vindex][cindex++]=u_compact[block][coord][(ivertex+icell)*9+(jvertex+jcell)*3+kvertex+kcell];
            vindex++;}}

    unsigned long start,stop;
    start=_rdtsc();

#if 0
    for(int coord=0;coord<3;coord++)
        for(int ivertex=0;ivertex<=1;ivertex++)
            for(int jvertex=0;jvertex<=1;jvertex++)
        for(int kvertex=0;kvertex<=1;kvertex++)
            for(int block=0;block<2;block++)
            for(int icell=0;icell<=1;icell++)
            for(int jcell=0;jcell<=1;jcell++)
            for(int kcell=0;kcell<=1;kcell++)
                u_duplicated_as_type2[coord][ivertex][jvertex][kvertex][block][icell][jcell][kcell]=
                    u_compact_as_type1[block][coord][ivertex+icell][jvertex+jcell][kvertex+kcell];
#else    
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
#endif

    stop=_rdtsc();
    std::cout<<"Ticks elapsed : "<<stop-start<<std::endl;

    bool inconsistencies_found=false;
    for(int i=0;i<3;i++)
        for(int j=0;j<8;j++)
            for(int k=0;k<16;k++)
                if(u_reference[i][j][k]!=u_duplicated[i][j][k]){
                    if(!inconsistencies_found)
                        std::cout<<"Inconsistencies found"<<std::endl;
                    inconsistencies_found=true;
                    std::cout<<"u_reference["<<i<<"]["<<j<<"]["<<k<<"]="<<u_reference[i][j][k]<<std::endl;
                    std::cout<<"u_duplicated["<<i<<"]["<<j<<"]["<<k<<"]="<<u_duplicated[i][j][k]<<std::endl;
                }
    if(!inconsistencies_found)
        std::cout<<"No inconsistencies"<<std::endl;
    

    

    return 0;
}

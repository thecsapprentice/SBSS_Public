
#include "string.h"
#include "time.h"
#include "stdlib.h"
#include "stdio.h" 
#include <malloc.h>
#include <iostream>
//#include <new>
#include <istream>
#include <fstream>
#include <vector>

#include "Dot_Product_Helper.h"


int main( int argc, char** argv)
{
    typedef float T;

    const int BASE = 10000000;
    const int DATASIZE = BASE * 16;

    T* data_a;
    T* data_b;

    posix_memalign((void**)&data_a,64,DATASIZE*sizeof(T));
    if(!data_a) exit(1);
    posix_memalign((void**)&data_b,64,DATASIZE*sizeof(T));
    if(!data_b) exit(1);

    int* offsets;
    posix_memalign((void**)&offsets,64,BASE*sizeof(int));
    if(!offsets) exit(1);
    
    for( int i = 0; i < 16*BASE; i++){
        data_a[i] = 1;
        data_b[i] = 2;
    }

    for( int i = 0; i < BASE; i++)
        offsets[i] = i*16;

    
    int number_of_partitions = 1;
    if(argc > 1 )
        number_of_partitions = atoi( argv[1] );

    unsigned long start_time,end_time;

    while(1){
        start_time=_rdtsc();
        double result = MT_Streaming_Kernels::Vector_Dot_Product_Helper<T, 16>( data_a, data_b, offsets, BASE, number_of_partitions );
        end_time=_rdtsc();
        std::cout << result << std::endl;
        unsigned long ticks=end_time-start_time;
        double elapsed_time=(double)ticks*1e-6l/(double)(CPU_MHZ);
        double bandwidth=(BASE*16.0l)*4/elapsed_time/1073741824;
        std::cout<<"[Saxpy] Clock ticks elapsed="<<ticks<<", Time="<<elapsed_time<<", Bandwidth="<<bandwidth<<"GB/s"<<std::endl;
    }

    
    
}


#include "PTHREAD_QUEUE.h"
//using namespace PhysBAM;

extern PhysBAM::PTHREAD_QUEUE* pthread_queue;

#include <cstdlib>
#include <immintrin.h>
#include <iomanip>
#include <iostream>
#include <sys/time.h>
#include "Grid_Compact.h"



struct timeval starttime, stoptime;
void
start_timer ()
{
  gettimeofday (&starttime, NULL);
}

void
stop_timer ()
{
  gettimeofday (&stoptime, NULL);
}

double
get_time ()
{
  return (double) stoptime.tv_sec - (double) starttime.tv_sec +
    (double) 1e-6 *(double) stoptime.tv_usec -
    (double) 1e-6 *(double) starttime.tv_usec;
}

template < class T > T Get_Random (const T a = (T) - 1., const T b = (T) 1.)
{
  return ((b - a) * (T) rand ()) / (T) RAND_MAX + a;
}


int main(int argc, char ** argv)
{

    pthread_queue=new PhysBAM::PTHREAD_QUEUE(12);
    const int domain_x = 256;
    const int domain_y = 256;
    const int domain_z = 256;

    std::cout<<"domain_x="<<domain_x<<std::endl;
    std::cout<<"domain_y="<<domain_y<<std::endl;
    std::cout<<"domain_z="<<domain_z<<std::endl;

    float* U[3];
    U[0] = new float[domain_x*domain_y*domain_z];
    U[1] = new float[domain_x*domain_y*domain_z];
    U[2] = new float[domain_x*domain_y*domain_z];

    float* P;
    P = new float[domain_x*domain_y*domain_z];

    for( int iter=0; iter<3; iter++)
        for( int index=0; index < domain_x*domain_y*domain_z; index++)
            U[iter][index] = Get_Random < float >();

    for( int index=0; index < domain_x*domain_y*domain_z; index++)
        P[index] = Get_Random < float >();
   

    const int X_Stride = domain_y * domain_z;
    const int Y_Stride = domain_z;

    int number_of_blocks = ((domain_x-1)/2)*((domain_y-1)/2)*((domain_z-1)/2);
    int number_of_partitions = 16;
    int* partition_offsets = new int[number_of_partitions];


    std::cout<<"number_of_blocks="<<number_of_blocks<<std::endl;
    std::cout<<"number_of_partitions="<<number_of_partitions<<std::endl;

    for(int partition=0;partition<number_of_partitions;partition++){
        std::cout<<"Partition #"<<partition<<" :"<<std::endl;
        // Compute optimal final block count


        int optimal_block_begin=(partition)*(number_of_blocks/number_of_partitions)+
            std::min(partition,number_of_blocks%number_of_partitions);


        int optimal_block_end=(partition+1)*(number_of_blocks/number_of_partitions)+
            std::min(partition+1,number_of_blocks%number_of_partitions);


        partition_offsets[partition] = optimal_block_begin;


        std::cout<<"Ideal beginning block = "<<optimal_block_begin<<std::endl;
        std::cout<<"Ideal ending block = "<<optimal_block_end<<std::endl;
    }


    typedef float u_compact_type[3][27];
    typedef float p_compact_type[8];
    u_compact_type *u_compact=new u_compact_type[number_of_blocks];
    p_compact_type *p_compact=new p_compact_type[number_of_blocks];
    int *block_offsets=new int[number_of_blocks];
    
    int counter=0;

    const int block_size=8;

    for( int block_i = 0; block_i < domain_x-2; block_i+=block_size)
    for( int block_j = 0; block_j < domain_y-2; block_j+=block_size)
    for( int block_k = 0; block_k < domain_z-2; block_k+=block_size)

    for( int i = block_i; i < block_i + block_size; i+=2)
    for( int j = block_j; j < block_j + block_size; j+=2)
    for( int k = block_k; k < block_k + block_size; k+=2)

        if((i < domain_x-2) && (j < domain_y-2) && (k < domain_z-2))


                block_offsets[counter++]=X_Stride * i + Y_Stride * j + k ;
    std::cout<<"Intended blocks = "<<number_of_blocks<<std::endl;
    std::cout<<"Initialized blocks = "<<counter<<std::endl;


    while(1){
        start_timer ();
        Grid_Compact<float>(u_compact,p_compact,U[0],U[1],U[2],P,block_offsets,number_of_blocks,
                            partition_offsets, number_of_partitions, X_Stride,Y_Stride, X_Stride,Y_Stride);
        stop_timer ();
        std::cout << get_time () << "s" << std::endl;
    }
    
    delete [] U[0];
    delete [] U[1];
    delete [] U[2];
    delete [] P;

}

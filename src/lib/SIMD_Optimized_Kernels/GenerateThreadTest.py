#!/usr/bin/python

#######################################################################
##  Copyright (c) 2011-2013 Nathan Mitchell, Eftychios Sifakis.
##  This file is covered by the FreeBSD license. Please refer to the 
##  license.txt file for more information.
#######################################################################


import os
import os.path
import re
import sys
import Generators.Generate as Generate
import Generators.GenerateParallel as GenerateParallel
import Generators.Parser as Parser

MAX_SUPPORTED_DIM = 16
SUPPORTED_MATERIALS = ("NEOHOOKEAN", "COROTATED")

StreamTest_Template = """
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <sys/time.h>
#include "KernelCommon.h"
%(MATERIAL_TAGS)s

#include "%(KERNEL_HEADER)s"

#include <Thread_Queueing/PTHREAD_QUEUE.h>
#include <Kernel_Serial_Base_Helper.h>


template<class T>
T Get_Random(const T a=(T)-1.,const T b=(T)1.)
{
    return ((b-a)*(T)rand())/(T)RAND_MAX+a;
}

struct timeval starttime,stoptime;
void start_timer(){gettimeofday(&starttime,NULL);}
void stop_timer(){gettimeofday(&stoptime,NULL);}
double get_time(){return (double)stoptime.tv_sec-(double)starttime.tv_sec+(double)1e-6*(double)stoptime.tv_usec-(double)1e-6*(double)starttime.tv_usec;}

%(SCALAR_INIT)s

#ifdef ENABLE_SSE_INSTRUCTION_SET
%(SSE_INIT)s
#endif

#ifdef ENABLE_AVX_INSTRUCTION_SET
%(AVX_INIT)s
#endif

#ifdef ENABLE_NEON_INSTRUCTION_SET
%(NEON_INIT)s
#endif

#ifdef ENABLE_MIC_INSTRUCTION_SET
%(MIC_INIT)s
#endif

int main(int argc,char* argv[])
{
    typedef float T;

    int seed=1;
    int threads=1;
    int threads_max=1;
    int passes=1;
    const int data_size=1000000;
    if(argc>=2) threads=atoi(argv[1]);
    if(argc>=3) threads_max=atoi(argv[2]);
    if(argc>=4) passes=atoi(argv[3]);
    srand(seed);

    pthread_queue=new PTHREAD_QUEUE(threads_max);

    std::cout << "Preparing to Run " << data_size << " of all kernels with " << threads << " threads." << std::endl;

    %(THREAD_TESTS)s

    return 0;

}

"""


Test_Template = """

{
std::cout << "Running Thread Test for %(KERNEL_NAME)s " << std::endl;

//=======================================================
//
//        DEFINE ALL VARIABLES USED BY KERNEL
//
//=======================================================
std::cout << "\\nAllocating all data: ";
std::cout.flush();

start_timer();
%(ARGUMENT_DEFS)s  

%(ARGUMENT_INIT)s
stop_timer();

std::cout << get_time() << "s\\n\\n"<< std::endl;


//=======================================================
//
//             COMPUTE SCALAR RESULTS
//
//=======================================================

{
std::cout << "\tRunning " << data_size << " of SCALAR :  "<< std::endl;

%(SCALAR_CALL)s

}

//=======================================================
//
//             COMPUTE SSE RESULTS
//
//=======================================================
#ifdef ENABLE_SSE_INSTRUCTION_SET
{
std::cout << "\tRunning " << data_size << " of SSE :  "<< std::endl;

%(SSE_CALL)s

}
#endif

//=======================================================
//
//             COMPUTE AVX RESULTS
//
//=======================================================
#ifdef ENABLE_AVX_INSTRUCTION_SET
{
std::cout << "\tRunning " << data_size << " of AVX :  "<< std::endl;

%(AVX_CALL)s

}
#endif

//=======================================================
//
//             COMPUTE NEON RESULTS
//
//=======================================================
#ifdef ENABLE_NEON_INSTRUCTION_SET
{
std::cout << "\tRunning " << data_size << " of NEON :  "<< std::endl;

%(NEON_CALL)s

}
#endif

//=======================================================
//
//             COMPUTE MIC RESULTS
//
//=======================================================
#ifdef ENABLE_MIC_INSTRUCTION_SET
{
std::cout << "\tRunning " << data_size << " of MIC :  "<< std::endl;

%(MIC_CALL)s

}
#endif

//=======================================================
//
//        FREE MEMORY USED BY ALL VARIABLES
//
//=======================================================
std::cout << "\\nFreeing all data: " << std::endl;
std::cout.flush();

%(ARGUMENT_FREE)s

}
"""


if __name__ == "__main__":
    
    if len(sys.argv) < 2:
        exit(1)
    
    KernelName = sys.argv[1]

    KernelHeader = open( os.path.join( "Kernels", KernelName, KernelName+".h"), 'r')
    HeaderArray = Parser.ParseHeader( KernelHeader )
    KernelHeader.close()

    # Generate Stream Test  

    Tests = ""
    SCALAR_Test_Init = ""
    SSE_Test_Init = ""
    AVX_Test_Init = ""
    NEON_Test_Init = ""
    MIC_Test_Init = ""

    for header in HeaderArray:
        if(header["Mat_Spec"]):
            materials = SUPPORTED_MATERIALS
        else:
            materials = [None,]
        for material in materials:
            Arguments = header["Args"]
            ArgDefs = GenerateParallel.GenArgDefs(Arguments,MAX_SUPPORTED_DIM)
            ArgInit = GenerateParallel.GenArgInit(Arguments,MAX_SUPPORTED_DIM)
            ArgFrees = GenerateParallel.GenArgFrees(Arguments,MAX_SUPPORTED_DIM)

            SCALAR_Test_Init += GenerateParallel.GenParallelOp(header["Name"],"SCALAR",Arguments,MAX_SUPPORTED_DIM,header["Mat_Spec"],material=material)
            SCALAR_Test_Call = GenerateParallel.GenParallelCall(header["Name"],"SCALAR",Arguments,MAX_SUPPORTED_DIM,header["Mat_Spec"],material=material)

            SSE_Test_Init += GenerateParallel.GenParallelOp(header["Name"],"SSE",Arguments,MAX_SUPPORTED_DIM,header["Mat_Spec"],material=material)
            SSE_Test_Call = GenerateParallel.GenParallelCall(header["Name"],"SSE",Arguments,MAX_SUPPORTED_DIM,header["Mat_Spec"],material=material)

            AVX_Test_Init += GenerateParallel.GenParallelOp(header["Name"],"AVX",Arguments,MAX_SUPPORTED_DIM,header["Mat_Spec"],material=material)
            AVX_Test_Call = GenerateParallel.GenParallelCall(header["Name"],"AVX",Arguments,MAX_SUPPORTED_DIM,header["Mat_Spec"],material=material)

            NEON_Test_Init += GenerateParallel.GenParallelOp(header["Name"],"NEON",Arguments,MAX_SUPPORTED_DIM,header["Mat_Spec"],material=material)
            NEON_Test_Call = GenerateParallel.GenParallelCall(header["Name"],"NEON",Arguments,MAX_SUPPORTED_DIM,header["Mat_Spec"],material=material)

            MIC_Test_Init += GenerateParallel.GenParallelOp(header["Name"],"MIC",Arguments,MAX_SUPPORTED_DIM,header["Mat_Spec"],material=material)
            MIC_Test_Call = GenerateParallel.GenParallelCall(header["Name"],"MIC",Arguments,MAX_SUPPORTED_DIM,header["Mat_Spec"],material=material)

            Tests +=  Test_Template % {
                "MAX_SUPPORTED_DIM":MAX_SUPPORTED_DIM,
                "KERNEL_NAME":header["Name"],
                "ARGUMENT_DEFS":ArgDefs,
                "ARGUMENT_INIT":ArgInit,
                "SCALAR_CALL":SCALAR_Test_Call,
                "SSE_CALL":SSE_Test_Call,
                "AVX_CALL":AVX_Test_Call,
                "NEON_CALL":NEON_Test_Call,
                "MIC_CALL":MIC_Test_Call,
                "ARGUMENT_FREE":ArgFrees,
                }
            
        
    print StreamTest_Template % {
        "MATERIAL_TAGS": " ".join(["struct %s_TAG;" % mat for mat in SUPPORTED_MATERIALS]),
        "KERNEL_HEADER" : KernelName+".h",
        "KERNEL_REF_HEADER" : KernelName+"_Reference.h",
        "SCALAR_INIT":SCALAR_Test_Init,
        "SSE_INIT":SSE_Test_Init,
        "AVX_INIT":AVX_Test_Init,
        "NEON_INIT":NEON_Test_Init,
        "MIC_INIT":MIC_Test_Init,
        "THREAD_TESTS" : Tests }

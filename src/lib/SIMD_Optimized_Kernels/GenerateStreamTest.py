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
import Generators.Parser as Parser

MAX_SUPPORTED_DIM = 16
SUPPORTED_MATERIALS = ("NEOHOOKEAN", "COROTATED")

StreamTest_Template = """
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <sys/time.h>
#include "KernelCommon.h"
%(MATERIAL_TAGS)s

#include "%(KERNEL_HEADER)s"

#define NUM_TRIALS 1000000

template<class T>
T Get_Random(const T a=(T)-1.,const T b=(T)1.)
{
    return ((b-a)*(T)rand())/(T)RAND_MAX+a;
}

struct timeval starttime,stoptime;
void start_timer(){gettimeofday(&starttime,NULL);}
void stop_timer(){gettimeofday(&stoptime,NULL);}
double get_time(){return (double)stoptime.tv_sec-(double)starttime.tv_sec+(double)1e-6*(double)stoptime.tv_usec-(double)1e-6*(double)starttime.tv_usec;}

int main(int argc,char* argv[])
{
    typedef float T;

    std::cout << "Preparing to Run " << NUM_TRIALS << " of all kernels." << std::endl;

    int seed=1;
    if(argc==2) seed=atoi(argv[1]);
    srand(seed);

    %(STREAM_TESTS)s

    return 0;

}

"""


Test_Template = """

{
std::cout << "Running Stream Test for %(KERNEL_NAME)s " << std::endl;

//=======================================================
//
//        DEFINE ALL VARIABLES USED BY KERNEL
//
//=======================================================

%(ARGUMENT_DEFS)s  

%(ARGUMENT_INIT)s

//=======================================================
//
//             COMPUTE SCALAR RESULTS
//
//=======================================================

{
std::cout << "\tRunning " << NUM_TRIALS << " of SCALAR :  ";
start_timer();
for(int n=0;n<NUM_TRIALS;n++)
{
%(SCALAR_CALL)s
}
stop_timer();
std::cout << get_time() << "s"<< std::endl;
}

//=======================================================
//
//             COMPUTE SSE RESULTS
//
//=======================================================

#ifdef ENABLE_SSE_INSTRUCTION_SET
{
std::cout << "\tRunning " << NUM_TRIALS << " of SSE :  ";
start_timer();
for(int n=0;n<NUM_TRIALS;n++)
{
%(SSE_CALL)s
}
stop_timer();
std::cout << get_time() << "s"<< std::endl;
}
#endif

//=======================================================
//
//             COMPUTE AVX RESULTS
//
//=======================================================

#ifdef ENABLE_AVX_INSTRUCTION_SET
{
std::cout << "\tRunning " << NUM_TRIALS << " of AVX :  ";
start_timer();
for(int n=0;n<NUM_TRIALS;n++)
{
%(AVX_CALL)s
}
stop_timer();
std::cout << get_time() << "s"<< std::endl;
}
#endif

//=======================================================
//
//             COMPUTE NEON RESULTS
//
//=======================================================

#ifdef ENABLE_NEON_INSTRUCTION_SET
{
std::cout << "\tRunning " << NUM_TRIALS << " of NEON :  ";
start_timer();
for(int n=0;n<NUM_TRIALS;n++)
{
%(NEON_CALL)s
}
stop_timer();
std::cout << get_time() << "s"<< std::endl;
}
#endif

//=======================================================
//
//             COMPUTE MIC RESULTS
//
//=======================================================

#ifdef ENABLE_MIC_INSTRUCTION_SET
{
std::cout << "\tRunning " << NUM_TRIALS << " of MIC :  ";
start_timer();
for(int n=0;n<NUM_TRIALS;n++)
{
%(MIC_CALL)s
}
stop_timer();
std::cout << get_time() << "s"<< std::endl;
}
#endif

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
    for header in HeaderArray:
        if(header["Mat_Spec"]):
            materials = SUPPORTED_MATERIALS
        else:
            materials = [None,]
        for material in materials:
            Arguments = header["Args"]
            ArgDefs = Generate.GenArgDefs(Arguments,MAX_SUPPORTED_DIM)
            ArgInit = Generate.GenArgInit(Arguments,MAX_SUPPORTED_DIM)

            SCALAR_Test_Call = Generate.GenTestCall(header["Name"],"SCALAR",Arguments,MAX_SUPPORTED_DIM,header["Mat_Spec"],material=material,use_compare_loops=False)
            SSE_Test_Call = Generate.GenTestCall(header["Name"],"SSE",Arguments,MAX_SUPPORTED_DIM,header["Mat_Spec"],material=material,use_compare_loops=False)
            AVX_Test_Call = Generate.GenTestCall(header["Name"],"AVX",Arguments,MAX_SUPPORTED_DIM,header["Mat_Spec"],material=material,use_compare_loops=False)
            NEON_Test_Call = Generate.GenTestCall(header["Name"],"NEON",Arguments,MAX_SUPPORTED_DIM,header["Mat_Spec"],material=material,use_compare_loops=False)
            MIC_Test_Call = Generate.GenTestCall(header["Name"],"MIC",Arguments,MAX_SUPPORTED_DIM,header["Mat_Spec"],material=material,use_compare_loops=False)

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
                }
            
        
    print StreamTest_Template % {
        "MATERIAL_TAGS": " ".join(["struct %s_TAG;" % mat for mat in SUPPORTED_MATERIALS]),
        "KERNEL_HEADER" : KernelName+".h",
        "KERNEL_REF_HEADER" : KernelName+"_Reference.h",
        "STREAM_TESTS" : Tests }

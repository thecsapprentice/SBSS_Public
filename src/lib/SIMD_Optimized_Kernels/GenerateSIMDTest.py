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
#include <iostream>
#include "KernelCommon.h"

%(MATERIAL_TAGS)s

#include "%(KERNEL_HEADER)s"
#include "%(KERNEL_REF_HEADER)s"

template<class T>
T Get_Random(const T a=(T)-1.,const T b=(T)1.)
{
    return ((b-a)*(T)rand())/(T)RAND_MAX+a;
}

int main(int argc,char* argv[])
{
    typedef float T;

    int seed=1;
    if(argc==2) seed=atoi(argv[1]);
    srand(seed);

    std::cout.precision(10);
    std::cout.setf(std::ios::fixed,std::ios::floatfield);

    %(SIMD_TESTS)s

    std::cout<<"SIMD check successful!"<<std::endl;

    return 0;

}

"""


Test_Template = """

{
std::cout << "Running SIMD Test for %(KERNEL_NAME)s " << std::endl;


//=======================================================
//
//        DEFINE ALL VARIABLES USED BY KERNEL
//
//=======================================================

%(ARGUMENT_DEFS)s  

%(ARGUMENT_INIT)s

//=======================================================
//
//             COMPUTE REFERENCE RESULTS
//
//=======================================================
 
%(REFERENCE_CALL)s

//=======================================================
//
//               COMPUTE SCALAR RESULTS
//
//=======================================================

{
%(SCALAR_TEST_CALL)s
}

//=======================================================
//
//               COMPUTE SSE RESULTS
//
//=======================================================

#ifdef ENABLE_SSE_INSTRUCTION_SET
{
%(SSE_TEST_CALL)s
}
#endif

//=======================================================
//
//               COMPUTE AVX RESULTS
//
//=======================================================

#ifdef ENABLE_AVX_INSTRUCTION_SET
{
%(AVX_TEST_CALL)s
}
#endif

//=======================================================
//
//               COMPUTE NEON RESULTS
//
//=======================================================

#ifdef ENABLE_NEON_INSTRUCTION_SET
{
%(NEON_TEST_CALL)s
}
#endif

//=======================================================
//
//               COMPUTE MIC RESULTS
//
//=======================================================

#ifdef ENABLE_MIC_INSTRUCTION_SET
{
%(MIC_TEST_CALL)s
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

    # Generate SIMD Test  

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
            Reference_Call = Generate.GenRefCall(header["Name"],Arguments,header["Mat_Spec"],material,array_init=True,max_final_dim=MAX_SUPPORTED_DIM,use_single=True)

            SCALAR_Test_Call = Generate.GenTestCall(header["Name"],"SCALAR",Arguments,MAX_SUPPORTED_DIM,header["Mat_Spec"],material=material,use_compare_loops=True)
            SSE_Test_Call = Generate.GenTestCall(header["Name"],"SSE",Arguments,MAX_SUPPORTED_DIM,header["Mat_Spec"],material=material,use_compare_loops=True)
            AVX_Test_Call = Generate.GenTestCall(header["Name"],"AVX",Arguments,MAX_SUPPORTED_DIM,header["Mat_Spec"],material=material,use_compare_loops=True)
            NEON_Test_Call = Generate.GenTestCall(header["Name"],"NEON",Arguments,MAX_SUPPORTED_DIM,header["Mat_Spec"],material=material,use_compare_loops=True)
            MIC_Test_Call = Generate.GenTestCall(header["Name"],"MIC",Arguments,MAX_SUPPORTED_DIM,header["Mat_Spec"],material=material,use_compare_loops=True)


            Tests +=  Test_Template % {
                "MAX_SUPPORTED_DIM":MAX_SUPPORTED_DIM,
                "KERNEL_NAME":header["Name"],
                "ARGUMENT_DEFS":ArgDefs,
                "ARGUMENT_INIT":ArgInit,
                #           "ARGUMENT_TYPEDEF":ArgTypeDefs,
                "SCALAR_TEST_CALL":SCALAR_Test_Call,
                "SSE_TEST_CALL":SSE_Test_Call,
                "AVX_TEST_CALL":AVX_Test_Call,
                "NEON_TEST_CALL":NEON_Test_Call,
                "MIC_TEST_CALL":MIC_Test_Call,
                "REFERENCE_CALL":Reference_Call,
                }
            
        
    print StreamTest_Template % {
        "MATERIAL_TAGS": " ".join(["struct %s_TAG;" % mat for mat in SUPPORTED_MATERIALS]),
        "KERNEL_HEADER" : KernelName+".h",
        "KERNEL_REF_HEADER" : KernelName+"_Reference.h",
        "SIMD_TESTS" : Tests }



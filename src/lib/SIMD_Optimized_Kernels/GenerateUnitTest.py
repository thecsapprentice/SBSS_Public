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

SUPPORTED_MATERIALS = ("NEOHOOKEAN", "COROTATED")

StreamTest_Template = """
#include <cstdlib>
#include <iostream>

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

    %(UNIT_TESTS)s

    return 0;

}

"""


Test_Template = """

{
%(ARGUMENT_DEFS)s  

%(ARGUMENT_INIT)s

%(TEST_CALL)s
%(REFERENCE_CALL)s
if( !(%(COMPARE_CALL)s) ){
   std::cout << "Failed to confirm unit test for %(KERNEL_NAME)s %(MATERIAL_STR)s" << std::endl;
   return 1;
}

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
            ArgDefs = Generate.GenArgDefs(Arguments,1)
            ArgInit = Generate.GenArgInit(Arguments,1)
            Test_Call = Generate.GenTestCall(header["Name"],"SINGLE",Arguments,1,header["Mat_Spec"],material=material)
            Reference_Call = Generate.GenRefCall(header["Name"],Arguments,header["Mat_Spec"],material)
            Compare_Call = Generate.GenCompareCall(header["Name"],Arguments,header["Mat_Spec"],material)

            Tests +=  Test_Template % {
                "KERNEL_NAME":header["Name"],
                "ARGUMENT_DEFS":ArgDefs,
                "ARGUMENT_INIT":ArgInit,
                "TEST_CALL":Test_Call,
                "REFERENCE_CALL":Reference_Call,
                "COMPARE_CALL":Compare_Call,
                "MATERIAL_STR":(header["Mat_Spec"])*("with material %s" % material),
                }
            
        
    print StreamTest_Template % {
        "MATERIAL_TAGS": " ".join(["struct %s_TAG;" % mat for mat in SUPPORTED_MATERIALS]),
        "KERNEL_HEADER" : KernelName+".h",
        "KERNEL_REF_HEADER" : KernelName+"_Reference.h",
        "UNIT_TESTS" : Tests }
        

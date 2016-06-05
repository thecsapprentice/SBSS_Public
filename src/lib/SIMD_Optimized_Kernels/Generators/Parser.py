#######################################################################
##  Copyright (c) 2011-2013 Nathan Mitchell, Eftychios Sifakis.
##  This file is covered by the FreeBSD license. Please refer to the 
##  license.txt file for more information.
#######################################################################


import re

ARG_MATCHER = re.compile("\s*(?P<HEAD>template<[A-Za-z 0-9, =_]+>)\s*(?P<IS_MATERIAL>(struct)|(void))\s*(?P<FunctionName>[A-Za-z_0-9]+)\s*(?:{\s*static\s*void\s*Run)?\s*\(\s*(?P<ARGS>([A-Za-z0-9 _\[\]&\(\)\*]*,?\s*)+)\)(?:};)?")

ARG_SUB_MATCHER = re.compile("\s*(?P<CONST>const)?\s*(?P<TYPE>[A-Za-z_0-9\*]+)\s*(?P<PTR>\*)?\s*\(?&?(?P<NAME>[A-Za-z_0-9]+)\)?\s*(?P<SIZE>(\[\d+\])*)")

def ParseHeader(file_object):
    HeaderString = ""

    for line in file_object:
        HeaderString += line

    Headers = HeaderString.replace( "\n",  " ").split(";")
    HeaderArray = []

    for header in Headers:
        m = ARG_MATCHER.search( header.strip() )
        if m:
            HeaderArray.append( {"Name":m.group("FunctionName"), "Args":m.group("ARGS").split(","),
                                 "Mat_Spec":(m.group("IS_MATERIAL")=="struct")} )

    for header in HeaderArray:
        header["Args"] = ParseArgs(header["Args"])

    return HeaderArray

def ParseArgs(args_raw):
    args_clean = []
        
    for arg in args_raw:
        arg_bits = ARG_SUB_MATCHER.search( arg ).groupdict()
        name = arg_bits["NAME"]
        size = arg_bits["SIZE"]
        ttype = arg_bits["TYPE"]
        size = size.replace("[","")
        size = size.replace("]"," ")
        size = size.strip()
        dims = size.split(" ")

        if len(dims) == 1 and dims[0] == '':
            dims = []

        if arg_bits["CONST"]:
            is_input = True
        else:
            is_input = False
        
        args_clean.append({"NAME":name, "DIMS": dims, "INPUT":is_input, "TYPE":ttype})

    return args_clean

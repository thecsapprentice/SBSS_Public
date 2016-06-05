#######################################################################
##  Copyright (c) 2011-2013 Nathan Mitchell, Eftychios Sifakis.
##  This file is covered by the FreeBSD license. Please refer to the 
##  license.txt file for more information.
#######################################################################


def GenArgDefs(args, max_final_dim,prefix="",postfix=""):
    defs = ""
    dim_str_bit = "[%s]"
    type_widths = {"SINGLE":"1",
                   "SCALAR":"1",
                   "SSE":"4",
                   "AVX":"8",
                   "NEON":"4",
                   "MIC":"16"}

    for arg in args:
        gen_dim_str = ""
        if (arg["TYPE"] == "T_DATA" or arg["TYPE"] == "I_DATA")  and max_final_dim>1:
            final = [max_final_dim,]
        else:
            final = []
            
        for dim in arg["DIMS"] + final:
            gen_dim_str += dim_str_bit%dim

        if arg["TYPE"] == "T_DATA":
            ttype = "T"
        elif arg["TYPE"] == "I_DATA":
            ttype = "int"
        else:
            ttype = arg["TYPE"]

        defs += "%(type)s %(name)s%(dims)s __attribute__((aligned(%(align)d)));\n" % { "name":prefix+arg["NAME"]+postfix, "dims":gen_dim_str, "type":ttype, "align":int(max_final_dim)*4}

        if not arg["INPUT"]:
            defs += "%(type)s %(name)s_reference%(dims)s __attribute__((aligned(%(align)d)));\n" % { "name":prefix+arg["NAME"]+postfix, "dims":gen_dim_str, "type":ttype, "align":int(max_final_dim)*4}
            defs += "%(type)s %(name)s_original%(dims)s __attribute__((aligned(%(align)d)));\n" % { "name":prefix+arg["NAME"]+postfix, "dims":gen_dim_str, "type":ttype, "align":int(max_final_dim)*4}

    return defs

def GenArgInit(args, max_final_dim):
    inits = ""
    init_iter = ("__a", "__b", "__c", "__d", "__e", "__f")
    dim_str_bit = "for(int %(iter)s=0;%(iter)s<%(dim_max)s;%(iter)s++) "
    dim_str_bit2 = "[%s]"

    INIT_MAP = {"T_DATA":"Get_Random<float>()",
                "I_DATA":"Get_Random<int>(1,99)",
                "float*":"new float[100]; for(int __x__=0;__x__<100;__x__++) %(NAME)s[__x__] = Get_Random<float>();",
                "bool":"true",
                }

    for arg in args:          
        gen_dim_str = ""
        dim_count = 0

        if (arg["TYPE"] == "T_DATA" or arg["TYPE"] == "I_DATA") and max_final_dim>1:
            final = [max_final_dim,]
        else:
            final = []

        for dim in arg["DIMS"]+final:
            gen_dim_str += dim_str_bit % {"iter":init_iter[dim_count], "dim_max":dim }
            dim_count += 1
        
        dim_str = ""
        dim_count = 0
        for dim in arg["DIMS"]+final:
            dim_str += dim_str_bit2 % init_iter[dim_count]
            dim_count += 1
        
        if not arg["INPUT"]:
            gen_dim_str+="{"

            gen_dim_str+="\n%s_original" % arg["NAME"]
            gen_dim_str+=dim_str
            gen_dim_str+="=%s;" % INIT_MAP[arg["TYPE"]]  % {"NAME": arg["NAME"]}

            gen_dim_str+="\n%s" % arg["NAME"]
            gen_dim_str+=dim_str
            gen_dim_str+="=%s_original%s;" %(arg["NAME"],dim_str)

            gen_dim_str+="\n%s_reference" % arg["NAME"]
            gen_dim_str+=dim_str
            gen_dim_str+="=%s_original%s;" %(arg["NAME"],dim_str)

            gen_dim_str+="}\n"

        else:
            gen_dim_str+="\n%s" % arg["NAME"]
            gen_dim_str+=dim_str
            gen_dim_str+="=%s;" % INIT_MAP[arg["TYPE"]]  %  {"NAME": arg["NAME"]}

        inits += gen_dim_str

    return inits






def GenTestCall(name,args):
    
    call = ""
    call += "%(name)s%(type)s(%(args)s);" % { "name":name,
                                              "type":"<T,T,1>",
                                              "args":",".join([arg[0] for arg in args]) }
    return call


def GenArgTypedefs(args,  max_final_dim):
    typedefs = ""
    
    refcount=1        
    for arg in args:
        dim_str = ""

        if (arg["TYPE"] == "T_DATA" or arg["TYPE"] == "I_DATA")  and max_final_dim>1:
            final = [max_final_dim,]
        else:
            final = []

        for dim in arg["DIMS"]+final:
            dim_str+="[%s]"%dim

        if arg["TYPE"] == "T_DATA":
            typedefs += "typedef T (&refArray%(refcount)d)%(SIZE)s;" % {"refcount":refcount, "SIZE":dim_str}
        elif arg["TYPE"] == "I_DATA":
            typedefs += "typedef int (&refArray%(refcount)d)%(SIZE)s;" % {"refcount":refcount, "SIZE":dim_str}
        else:
            typedefs += "typedef %(ttype)s (&refArray%(refcount)d)%(SIZE)s;" % {"refcount":refcount, "SIZE":dim_str, "ttype":arg["TYPE"]}

        refcount+=1

    return typedefs

def GenArgTestInit(args, max_final_dim, use_compare_loops=True):
    inits = ""
    init_iter = ("__a", "__b", "__c", "__d", "__e", "__f")
    dim_str_bit = "for(int %(iter)s=0;%(iter)s<%(dim_max)s;%(iter)s++) "
    dim_str_bit2 = "[%s]"

    for arg in args:
        if arg["INPUT"] or use_compare_loops == False:
            continue

        gen_dim_str = ""
        dim_count = 0
        for dim in arg["DIMS"]+[max_final_dim,]:
            if(dim==None):
                continue
            gen_dim_str += dim_str_bit % {"iter":init_iter[dim_count], "dim_max":dim }
            dim_count += 1

        dim_str = ""
        dim_count = 0
        for dim in arg["DIMS"]+[max_final_dim,]:
            if(dim==None):
                continue
            dim_str += dim_str_bit2 % init_iter[dim_count]
            dim_count += 1

        gen_dim_str+="%s" % arg["NAME"]
        gen_dim_str+=dim_str
        gen_dim_str+="=%s_original" % arg["NAME"]
        gen_dim_str+=dim_str
        gen_dim_str+=";"

        inits += gen_dim_str

    return inits


def GenTestCall(name,t,args,max_final_dim,is_material_spec,material=None,use_compare_loops=False):

    type_widths = {"SINGLE":"1",
                   "SCALAR":"1",
                   "SSE":"4",
                   "AVX":"8",
                   "NEON":"4",
                   "MIC":"16"}

    data_widths = {"SINGLE":"1",
                   "SCALAR":str(max_final_dim),
                   "SSE":str(max_final_dim),
                   "AVX":str(max_final_dim),
                   "NEON":str(max_final_dim),
                   "MIC":str(max_final_dim)}

    if is_material_spec:
        type_args = {"SINGLE":"<%(T_MATERIAL_TAG)s,float,float,int>",
                     "SCALAR":"<%%(T_MATERIAL_TAG)s,float,float[%(size)s],int[%(size)s]>" % {"size" : data_widths["SCALAR"]},
                     "SSE":"<%%(T_MATERIAL_TAG)s,__m128,float[%(size)s],int[%(size)s]>" % {"size" : data_widths["SSE"]},
                     "AVX":"<%%(T_MATERIAL_TAG)s,__m256,float[%(size)s],int[%(size)s]>" % {"size" : data_widths["AVX"]},
                     "NEON":"<%%(T_MATERIAL_TAG)s,float32x4_t,float[%(size)s],int[%(size)s]>" % {"size" : data_widths["NEON"]},
                     "MIC":"<%%(T_MATERIAL_TAG)s,__m512,float[%(size)s],int[%(size)s]>" % {"size" : data_widths["MIC"]}}
    else:
        type_args = {"SINGLE":"<float,float,int>",
                     "SCALAR":"<float,float[%(size)s],int[%(size)s]>" % {"size" : data_widths["SCALAR"]},
                     "SSE":"<__m128,float[%(size)s],int[%(size)s]>" % {"size" : data_widths["SSE"]},
                     "AVX":"<__m256,float[%(size)s],int[%(size)s]>" % {"size" : data_widths["AVX"]},
                     "NEON":"<float32x4_t,float[%(size)s],int[%(size)s]>" % {"size" : data_widths["NEON"]},
                     "MIC":"<__m512,float[%(size)s],int[%(size)s]>" % {"size" : data_widths["MIC"]}}


    call = ""

    if(int(data_widths[t])>1):
        call += GenArgTypedefs(args, max_final_dim)  

    if(int(data_widths[t])>1):
        call += GenArgTestInit(args, data_widths[t], use_compare_loops)
    else:
        call += GenArgTestInit(args, None, use_compare_loops)

    call += """for(int i=0; i<%s; i+=%s) {""" % (data_widths[t], type_widths[t])

    if(int(data_widths[t])>1):
        refcount=1        
        for arg in args:
            dim_str = ""                       
            if (arg["TYPE"] == "T_DATA" or arg["TYPE"] == "I_DATA") :
                for dim in arg["DIMS"]:
                    dim_str+="[0]"      
                dim_str+="[i]"
            else:
                for dim in arg["DIMS"]:
                    dim_str+="[0]"           
                
            call += "refArray%(refcount)d %(name)sk=reinterpret_cast<refArray%(refcount)d>(%(name)s%(SIZE)s);" % {"name":arg["NAME"],
                                                                                                                  "refcount":refcount,
                                                                                                                  "SIZE":dim_str}
            refcount+=1
    
        
    if is_material_spec:
        call += "%(name)s%(type)s::Run(%(args)s);" % { "name":name,
                                                       "type":type_args[t] % {"T_MATERIAL_TAG":material+"_TAG"},
                                                       "args":",".join([arg["NAME"]+((int(data_widths[t])>1)*"k") for arg in args]) }
    else:
        call += "%(name)s%(type)s(%(args)s);" % { "name":name,
                                                  "type":type_args[t],
                                                  "args":",".join([arg["NAME"]+((int(data_widths[t])>1)*"k") for arg in args]) }       

    
    call += """}\n"""

    if use_compare_loops:
        call += GenCompareLoops(t,args, data_widths[t])

    return call


def GenRefCall(name,args,is_material_spec,material=None,array_init=False,max_final_dim=1,use_single=False):
    call = ""

    if array_init:

        call += GenArgDefs(args,1,prefix="__m")

        call += "for( int k=0;k<%d;k++){" % max_final_dim

        inits = ""
        init_iter = ("__a", "__b", "__c", "__d", "__e", "__f")
        dim_str_bit = "for(int %(iter)s=0;%(iter)s<%(dim_max)s;%(iter)s++) "
        dim_str_bit2 = "[%s]"

        for arg in args:          
            gen_dim_str = ""
            dim_count = 0
            
            final = []
            if (arg["TYPE"] == "T_DATA" or arg["TYPE"] == "I_DATA") :
                post = "[k]"
            else:
                post = ""

            for dim in arg["DIMS"]+final:
                gen_dim_str += dim_str_bit % {"iter":init_iter[dim_count], "dim_max":dim }
                dim_count += 1
        
            dim_str = ""
            dim_count = 0
            for dim in arg["DIMS"]+final:
                dim_str += dim_str_bit2 % init_iter[dim_count]
                dim_count += 1
        
            if not arg["INPUT"]:
                gen_dim_str+=""

                gen_dim_str+="\n__m%s_reference" % arg["NAME"]
                gen_dim_str+=dim_str
                gen_dim_str+="=%s_reference%s;" %(arg["NAME"],dim_str+post)
                
                gen_dim_str+="\n"
                
            else:
                gen_dim_str+="\n__m%s" % arg["NAME"]
                gen_dim_str+=dim_str
                gen_dim_str+="=%s%s;" % (arg["NAME"],dim_str+post)
                
            call += gen_dim_str

            
        if use_single:
            if is_material_spec:
                call += "%(name)s%(type)s::Run(%(args)s);" % { "name":name,
                                                               "type":"<%(T_MATERIAL_TAG)s,float,float,int>" % {"T_MATERIAL_TAG":material+"_TAG"},
                                                               "args":",".join(["__m"+arg["NAME"]+(not arg["INPUT"])*"_reference" for arg in args]) }
            else:
                call += "%(name)s%(type)s(%(args)s);" % { "name":name,
                                                          "type":"<float,float,int>",
                                                          "args":",".join(["__m"+arg["NAME"]+(not arg["INPUT"])*"_reference" for arg in args]) }       
                
                
        else:
            if(is_material_spec):
                call += "%(name)s_%(material)s_Reference%(type)s(%(args)s);" % { "name":name,
                                                                                 "material":material.capitalize(),
                                                                                 "type":"<float>",
                                                                                 "args":",".join(["__m"+arg["NAME"]+(not arg["INPUT"])*"_reference" for arg in args]) }
                
            else:
                call += "%(name)s_Reference%(type)s(%(args)s);" % { "name":name,
                                                                    "type":"<float>",
                                                                    "args":",".join(["__m"+arg["NAME"]+(not arg["INPUT"])*"_reference" for arg in args]) }
                
        for arg in args:          
            gen_dim_str = ""
            dim_count = 0
            
            final = []
            if (arg["TYPE"] == "T_DATA" or arg["TYPE"] == "I_DATA") :
                post = "[k]"
            else:
                post = ""

            for dim in arg["DIMS"]+final:
                gen_dim_str += dim_str_bit % {"iter":init_iter[dim_count], "dim_max":dim }
                dim_count += 1
        
            dim_str = ""
            dim_count = 0
            for dim in arg["DIMS"]+final:
                dim_str += dim_str_bit2 % init_iter[dim_count]
                dim_count += 1
        
            if not arg["INPUT"]:
                gen_dim_str+=""

                gen_dim_str+="%s_reference%s=" %(arg["NAME"],dim_str+post)
                gen_dim_str+="\n__m%s_reference" % arg["NAME"]
                gen_dim_str+=dim_str + ";"
                
                
                gen_dim_str+="\n"
                
                call += gen_dim_str
            
        call += "}"

    else:

        if(is_material_spec):
            call += "%(name)s_%(material)s_Reference%(type)s(%(args)s);" % { "name":name,
                                                                             "material":material.capitalize(),
                                                                             "type":"<float>",
                                                                             "args":",".join([arg["NAME"]+(not arg["INPUT"])*"_reference" for arg in args]) }
            
        else:
            call += "%(name)s_Reference%(type)s(%(args)s);" % { "name":name,
                                                                "type":"<float>",
                                                                "args":",".join([arg["NAME"]+(not arg["INPUT"])*"_reference" for arg in args]) }
            



    return call



def GenCompareCall(name,args,is_material_spec,material=None):
    
    call = ""
    if(is_material_spec):
        call += "%(name)s_%(material)s_Compare%(type)s(%(args)s,%(refargs)s)" % { "name":name,
                                                                                  "type":"<float>",
                                                                                  "material":material.capitalize(),
                                                                                  "args":",".join([arg["NAME"] for arg in args if not arg["INPUT"]]),
                                                                                  "refargs":",".join([arg["NAME"]+"_reference" for arg in args if not arg["INPUT"]])}
    else:
        call += "%(name)s_Compare%(type)s(%(args)s,%(refargs)s)" % { "name":name,
                                                                     "type":"<float>",
                                                                     "args":",".join([arg["NAME"] for arg in args if not arg["INPUT"]]),
                                                                     "refargs":",".join([arg["NAME"]+"_reference" for arg in args if not arg["INPUT"]])}

    return call


def GenCompareLoops(t,args, max_final_dim):
    inits = ""
    init_iter = ("__a", "__b", "__c", "__d", "__e", "__f")
    dim_str_bit = "for(int %(iter)s=0;%(iter)s<%(dim_max)s;%(iter)s++) "
    dim_str_bit2 = "[%s]"
    err_bit = """<<", %s="<<%s"""

    for arg in args:
        if arg["INPUT"]:
            continue 
        gen_dim_str = ""
        dim_count = 0
        for dim in arg["DIMS"]+[max_final_dim,]:
            gen_dim_str += dim_str_bit % {"iter":init_iter[dim_count], "dim_max":dim }
            dim_count += 1

        dim_bits = ""

        dim_count = 0
        for dim in arg["DIMS"]+[max_final_dim,]:
            dim_bits += dim_str_bit2 % init_iter[dim_count]
            dim_count += 1

        gen_dim_str+="if(std::abs((%(name)s%(dims)s - %(name)s_reference%(dims)s) / (%(name)s_reference%(dims)s)) > 1 ){" % { "name": arg["NAME"], "dims" : dim_bits }
      
        inits += gen_dim_str

        inits += """std::cerr<<"Mismatch detected in %s implementation"<<std::endl;""" % t
        inits += """std::cerr<<"Variable %s:"<<std::endl;""" % arg["NAME"]
        inits += """std::cerr<<"seed="<<seed"""  
        
        dim_count = 0
        for dim in arg["DIMS"]+[max_final_dim,]:
            inits += err_bit % (init_iter[dim_count],init_iter[dim_count])
            dim_count += 1

        inits +="""<<std::endl;"""

        dims = ""
        dim_count = 0
        for dim in arg["DIMS"]+[max_final_dim,]:
            dims += dim_str_bit2 % init_iter[dim_count]
            dim_count += 1

        inits += """std::cerr<<"%(name)s %(type)s=  "<<%(name)s%(dims)s<<std::endl;"""  % {
            "type":t,
            "name":arg["NAME"],
            "dims":dims
            }

        inits += """std::cerr<<"%(name)s Reference=  "<<%(name)s_reference%(dims)s<<std::endl;""" % {
            "type":t,
            "name":arg["NAME"],
            "dims":dims
            }

        inits += """std::cerr<<"%(name)s Rel Difference=  "<< std::abs((%(name)s%(dims)s - %(name)s_reference%(dims)s) / (%(name)s_reference%(dims)s)) << std::endl;""" % {
            "type":t,
            "name":arg["NAME"],
            "dims":dims
            }

        inits += """std::cerr<<"%(name)s Abs Difference=  "<< std::abs(%(name)s%(dims)s - %(name)s_reference%(dims)s) << std::endl;""" % {
            "type":t,
            "name":arg["NAME"],
            "dims":dims
            }

        inits += """return 1;}\n"""

        

    return inits

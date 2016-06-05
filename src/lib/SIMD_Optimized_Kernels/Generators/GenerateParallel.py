#######################################################################
##  Copyright (c) 2011-2013 Nathan Mitchell, Eftychios Sifakis.
##  This file is covered by the FreeBSD license. Please refer to the 
##  license.txt file for more information.
#######################################################################

def GenOpFullTypedefs(args,  max_final_dim):
    typedefs = ""
    
    refcount=1        
    for arg in args:
        dim_str = ""

        if (arg["TYPE"] == "T_DATA" or arg["TYPE"] == "I_DATA")  and max_final_dim>1:
            final = [max_final_dim,]
        else:
            final = []

        for dim in ["SIZE",]+arg["DIMS"]+final:
            dim_str+="[%s]"%dim

        if arg["TYPE"] == "T_DATA":
            typedefs += "typedef float (&fullArray%(refcount)d)%(SIZE)s;" % {"refcount":refcount, "SIZE":dim_str}
        elif arg["TYPE"] == "I_DATA":
            typedefs += "typedef int (&fullArray%(refcount)d)%(SIZE)s;" % {"refcount":refcount, "SIZE":dim_str}
        else:
            typedefs += "typedef %(ttype)s (&fullArray%(refcount)d)%(SIZE)s;" % {"refcount":refcount, "SIZE":dim_str, "ttype":arg["TYPE"]}

        refcount+=1

    return typedefs


def GenOpFullExtractions(args,  max_final_dim):
    typedefs = ""
    
    refcount=1        
    for arg in args:
        typedefs += "fullArray%(refcount)d _r%(name)s =reinterpret_cast<fullArray%(refcount)d>(*_local_%(name)s);" % {"name":arg["NAME"],"refcount":refcount}
        refcount+=1
    return typedefs


def GenChunkTypedefs(args,  max_final_dim):
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
            typedefs += "typedef float (&refArray%(refcount)d)%(SIZE)s;" % {"refcount":refcount, "SIZE":dim_str}
        elif arg["TYPE"] == "I_DATA":
            typedefs += "typedef int (&refArray%(refcount)d)%(SIZE)s;" % {"refcount":refcount, "SIZE":dim_str}
        else:
            typedefs += "typedef %(ttype)s (&refArray%(refcount)d)%(SIZE)s;" % {"refcount":refcount, "SIZE":dim_str, "ttype":arg["TYPE"]}

        refcount+=1

    return typedefs

def GenChunkExtractions(args,  max_final_dim):
    typedefs = ""
    refcount=1        
    for arg in args:
        dim_str = ""                       
        if (arg["TYPE"] == "T_DATA" or arg["TYPE"] == "I_DATA") :
            for dim in arg["DIMS"]:
                dim_str+="[0]"      
            dim_str+="[chunk_offset]"
        else:
            for dim in arg["DIMS"]:
                dim_str+="[0]"           
                
        typedefs += "refArray%(refcount)d %(name)sk=reinterpret_cast<refArray%(refcount)d>(_r%(name)s[index]%(SIZE)s);" % {"name":arg["NAME"],
                                                                                                                  "refcount":refcount,
                                                                                                                  "SIZE":dim_str}
        refcount+=1

    return typedefs


def type_translate(in_t):
    if in_t == "T_DATA":
        return "float"
    
    if in_t == "I_DATA":
        return "int"
    
    return in_t


def GenParallelOp(name,t,args,max_dim,is_material_spec,material):
    op_template = """
template <int SIZE> 
class %(KERNEL_NAME)s_%(TYPE)s_%(MATERIAL)s
    {
    private:
        // Generate Variables Here
        %(OP_LOCALS)s

    public:
        explicit %(KERNEL_NAME)s_%(TYPE)s_%(MATERIAL)s(%(ARGS)s) : %(CONSTRUCTOR_LIST)s {}
        void Execute(int index)
        {
        // full array typedefs
        //typedef int (&refArray)[SIZE][3][8];
        %(FULL_TYPEDEFS)s

        // full array extractions
        //refArray _rA = reinterpret_cast < refArray >(*_A);
        %(FULL_EXTRACTIONS)s

        const int ChunkSize = %(CHUNK)s;
        // chunk typedef
        //typedef int (&refArray1)[3][8];
        %(CHUNK_TYPEDEFS)s

        for( int chunk_offset=0; chunk_offset<%(MAX_DIM)s; chunk_offset+=ChunkSize)
            {
                //refArray1 Ak = reinterpret_cast < refArray1 > (_rA[entry][0][chunk_offset]);
                %(CHUNK_EXTRACTIONS)s
               
                %(CALL)s
             }

        }
    };
"""

    type_widths = {"SINGLE":"1",
                   "SCALAR":"1",
                   "SSE":"4",
                   "AVX":"8",
                   "NEON":"4",
                   "MIC":"16"}

    data_widths = {"SINGLE":"1",
                   "SCALAR":str(max_dim),
                   "SSE":str(max_dim),
                   "AVX":str(max_dim),
                   "NEON":str(max_dim),
                   "MIC":str(max_dim)}

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

    op_locals = []
    op_args = []
    op_constr = []
    for arg in args: 
        op_locals.append("%(type)s * _local_%(name)s;" % { "type":type_translate(arg["TYPE"]),
                                                       "name":arg["NAME"] })

        op_args.append("%(type)s * %(name)s_in" % { "type":type_translate(arg["TYPE"]),
                                                       "name":arg["NAME"] } )

        op_constr.append("_local_%(name)s(%(name)s_in)" % { "name":arg["NAME"] } )

    full_typedefs = GenOpFullTypedefs(args,  max_dim)
    full_extractions = GenOpFullExtractions(args,  max_dim)
    chunk_typedefs = GenChunkTypedefs(args,  max_dim)
    chunk_extractions = GenChunkExtractions(args,  max_dim)

    call = ""
    if is_material_spec:
        call = "%(name)s%(type)s::Run(%(args)s);" % { "name":name,
                                                      "type":type_args[t] % {"T_MATERIAL_TAG":material+"_TAG"},
                                                      "args":",".join([arg["NAME"]+((int(data_widths[t])>1)*"k") for arg in args]) }
    else:
        call = "%(name)s%(type)s(%(args)s);" % { "name":name,
                                                 "type":type_args[t],
                                                 "args":",".join([arg["NAME"]+((int(data_widths[t])>1)*"k") for arg in args]) }       
        

    return op_template % { "KERNEL_NAME":name,
                           "TYPE":t,
                           "MATERIAL":str(material),
                           "MAX_DIM":max_dim,
                           "OP_LOCALS":" ".join(op_locals),
                           "ARGS":",".join(op_args),
                           "CONSTRUCTOR_LIST":",".join(op_constr),
                           "CHUNK":type_widths[t],
                           "FULL_TYPEDEFS":full_typedefs,
                           "FULL_EXTRACTIONS":full_extractions,
                           "CHUNK_TYPEDEFS":chunk_typedefs,
                           "CHUNK_EXTRACTIONS":chunk_extractions,
                           "CALL":call}

def GenParallelCall(name,t,args,max_dim,is_material_spec,material):
    op_template = """
    %(KERNEL_NAME)s_%(TYPE)s_%(MATERIAL)s<data_size> op(%(INPUT_ARGS)s);

    for( int t = threads; t <= threads_max; t+=std::max<int>(((threads_max - threads)/30), 1)){
         std::cout << "Running Test with " << t << " threads." << std::endl;
         MT_Streaming_Kernels::Kernel_Serial_Base_Helper<%(KERNEL_NAME)s_%(TYPE)s_%(MATERIAL)s<data_size> > helper(op,data_size,t);

         double min_time = 10000000;
         double max_time = -1;
         double avg_time = 0;
                  
         for(int i=0; i<passes; i++){   
           start_timer();
           helper.Run_Parallel();
           stop_timer();
           std::cout << get_time() << "s"<< std::endl;
           min_time = std::min<double>( min_time, get_time() );
           max_time = std::max<double>( max_time, get_time() );
           avg_time += get_time();
         }
         avg_time = avg_time / passes;
         std::cout << "Min pass time: " << min_time << std::endl;
         std::cout << "Max pass time: " << max_time << std::endl;
         std::cout << "Avg pass time: " << avg_time << std::endl;
     }


"""

    op_in_args = []
    for arg in args: 
 
        op_in_args.append("(%(type)s *)&%(name)s" % { "type":type_translate(arg["TYPE"]),
                                                            "name":arg["NAME"] } )
    
    return op_template % { "KERNEL_NAME":name,
                           "TYPE":t,
                           "MATERIAL":str(material),
                           "INPUT_ARGS":",".join(op_in_args)}



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
            final = [str(max_final_dim),]
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



        defs += "typedef %(type)s (&%(name)s_type)[data_size]%(dims)s; %(name)s_type %(name)s = reinterpret_cast<%(name)s_type>(*((T*)(_mm_malloc(%(dims_merged)s*sizeof(%(type)s),%(align)d))));\n" % { "name":prefix+arg["NAME"]+postfix,
                                                                                                                                                                                                         "dims":gen_dim_str,
                                                                                                                                                                                                         "dims_merged":"*".join(["data_size",]+arg["DIMS"]+final),
                                                                                                                                                                                                         "type":ttype,
                                                                                                                                                                                                         "align":int(max_final_dim)*4}

        #defs += "%(type)s *%(name)s = new %(type)s[data_size*%(dims)s];\n" % { "name":prefix+arg["NAME"]+postfix, "dims":gen_dim_str, "type":ttype}

        #if not arg["INPUT"]:
        #    defs += "%(type)s *%(name)s_reference = new %(type)s[data_size]%(dims)s;\n" % { "name":prefix+arg["NAME"]+postfix, "dims":gen_dim_str, "type":ttype}
        #    defs += "%(type)s *%(name)s_original = new %(type)s[data_size]%(dims)s;\n" % { "name":prefix+arg["NAME"]+postfix, "dims":gen_dim_str, "type":ttype}

    return defs

def GenArgFrees(args, max_final_dim,prefix="",postfix=""):
    frees = ""
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
            final = [str(max_final_dim),]
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

        frees += "_mm_free( reinterpret_cast< void* >( %(name)s ));\n" % { "name" : prefix+arg["NAME"]+postfix }

    return frees


def GenArgInit(args, max_final_dim):
    inits = ""
    init_iter = ("__a", "__b", "__c", "__d", "__e", "__f")
    dim_str_bit = "for(int %(iter)s=0;%(iter)s<%(dim_max)s;%(iter)s++) "
    dim_str_bit3 = "for(int %(iter)s=0;%(iter)s<%(dim_max)s;%(iter)s++){ "
    dim_str_bit2 = "[%s]"

    INIT_MAP = {"T_DATA":"Get_Random<float>()",
                "I_DATA":"Get_Random<int>(1,99)",
                "float*":"new float[100]; for(int __x__=0;__x__<100;__x__++) %(NAME)s%(DIM)s[__x__] = Get_Random<float>();",
                "bool":"true",
                }

    for arg in args:          
        gen_dim_str = ""
        dim_count = 0

        if (arg["TYPE"] == "T_DATA" or arg["TYPE"] == "I_DATA") and max_final_dim>1:
            final = [max_final_dim,]
        else:
            final = []

        max_count = len(["data_size",]+arg["DIMS"]+final)
        for dim in ["data_size",]+arg["DIMS"]+final:
            if dim_count+1 == max_count:
                gen_dim_str += dim_str_bit3 % {"iter":init_iter[dim_count], "dim_max":dim }
            else:
                gen_dim_str += dim_str_bit % {"iter":init_iter[dim_count], "dim_max":dim }
            dim_count += 1
        
        dim_str = ""
        dim_count = 0
        for dim in ["data_size",]+arg["DIMS"]+final:
            dim_str += dim_str_bit2 % init_iter[dim_count]
            dim_count += 1
        
        gen_dim_str+="\n%s" % arg["NAME"]
        gen_dim_str+=dim_str
        gen_dim_str+="=%s;}" % INIT_MAP[arg["TYPE"]]  %  {"NAME": arg["NAME"], "DIM":dim_str}

        inits += gen_dim_str

    return inits

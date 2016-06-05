include ../Makefile.Library.PreMatter

CXXFLAGS += -I$(INCLUDE_INSTALL_DIRECTORY)/PhysBAM/Public_Library -I$(INCLUDE_INSTALL_DIRECTORY)
LD_FLAGS += -L$(LIBRARY_INSTALL_DIRECTORY)  
DEPENDENCIES =  -lTQ -lSIMD_Optimized_Kernels


Library_Src = Vector_SAXPY_V9/Vector_SAXPY_V9_Helper.cpp \
                Project/Project_Helper.cpp \
                Vector_Subtract/Vector_Subtract_Helper.cpp \
                Vector_SAXPY_V8/Vector_SAXPY_V8_Helper.cpp \
                Vector_SAXPY_A/Vector_SAXPY_A_Helper.cpp \
                Vector_Times/Vector_Times_Helper.cpp \
                Vector_SAXPY_AB/Vector_SAXPY_AB_Helper.cpp \
                Convergence_Norm/Convergence_Norm_Helper.cpp \
                Vector_SAXPY_V/Vector_SAXPY_V_Helper.cpp \
                Dot_Product/Dot_Product_Helper.cpp \
                Vector_SAXPY_V6/Vector_SAXPY_V6_Helper.cpp \
                Vector_Add/Vector_Add_Helper.cpp \
                Vector_Set/Vector_Set_Helper.cpp

Library_Headers = Vector_SAXPY_V9/Vector_SAXPY_V9_Helper.h \
                  Project/Project_Helper.h \
                  Vector_Subtract/Vector_Subtract_Helper.h \
                  Vector_SAXPY_V8/Vector_SAXPY_V8_Helper.h \
                  Vector_SAXPY_A/Vector_SAXPY_A_Helper.h \
                  Vector_Times/Vector_Times_Helper.h \
                  Vector_SAXPY_AB/Vector_SAXPY_AB_Helper.h \
                  Convergence_Norm/Convergence_Norm_Helper.h \
                  Vector_SAXPY_V/Vector_SAXPY_V_Helper.h \
                  Kernel_Base/Kernel_Base_Reducer_Helper.h \
                  Kernel_Base/Kernel_Base_Helper.h \
                  Dot_Product/Dot_Product_Helper.h \
                  Vector_SAXPY_V6/Vector_SAXPY_V6_Helper.h \
                  Vector_Add/Vector_Add_Helper.h \
                  Kernel_Base/Kernel_Mesh_Base_Reducer_Helper.h \
                  Kernel_Base/Kernel_Mesh_Base_Helper.h \
                  Kernel_Base/Kernel_Serial_Base_Reducer_Helper.h \
                  Kernel_Base/Kernel_Serial_Base_Helper.h \
                  Kernel_Base/Kernel_OpenMP_Base_Reducer_Helper.h \
                  Kernel_Base/Kernel_OpenMP_Base_Helper.h \
                  Vector_Set/Vector_Set_Helper.h \
                  Vector_Element_Times/Vector_Element_Times_Helper.h

LIBRARY_NAME=libCGKernels
INCLUDE_NAME=CG_Optimized_Kernels

include ../Makefile.Library.Template




include ../Makefile.Library.PreMatter

CXXFLAGS += -DFP_ARCH=$(FP_ARCH) -I$(INCLUDE_INSTALL_DIRECTORY)/PhysBAM/Public_Library -I$(INCLUDE_INSTALL_DIRECTORY) -std=gnu++0x 


ifdef ENABLE_MIC
CXXFLAGS += -DENABLE_MIC=1
endif

ifdef ENABLE_PHYSBAM_IO
CXXFLAGS += -DENABLE_PHYSBAM_IO=1
endif

LD_FLAGS += -L$(LIBRARY_INSTALL_DIRECTORY)  
DEPENDENCIES = -lEngineInterface -lEngineBackend
EXT_DEPENDENCIES = -lPhysBAM_NewStyle


Library_Src=$(shell find . -name "*.cpp") 
Library_Headers=$(shell find . -name "*.h") 

LIBRARY_NAME=libEngineInterfaceLocal
INCLUDE_NAME=EngineInterfaceLocal

include ../Makefile.Library.Template

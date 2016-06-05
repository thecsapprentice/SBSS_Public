include ../Makefile.Library.PreMatter

CXXFLAGS += -DFP_ARCH=$(FP_ARCH) -I$(INCLUDE_INSTALL_DIRECTORY)/PhysBAM/Public_Library -I$(INCLUDE_INSTALL_DIRECTORY)/PhysBAM/Common_Libraries -I$(INCLUDE_INSTALL_DIRECTORY) -std=gnu++0x 

ifdef ENABLE_MIC
CXXFLAGS += -DENABLE_MIC=1
endif

ifdef ENABLE_PHYSBAM_IO
CXXFLAGS += -DENABLE_PHYSBAM_IO=1
endif

LD_FLAGS += -L$(LIBRARY_INSTALL_DIRECTORY) 

DEPENDENCIES =  -lTQ -lCutter -lCommon
EXT_DEPENDENCIES = -lPhysBAM_NewStyle -lCommon_Geometry
ifdef USE_MPI
DEPENDENCIES += -lEngineInterfaceMPI
else
DEPENDENCIES += -lEngineInterfaceLocal
endif

Library_Src=$(shell find . -name "*.cpp") 
Library_Headers=$(shell find . -name "*.h") 

LIBRARY_NAME=libEngineFrontend
INCLUDE_NAME=EngineFrontend

include ../Makefile.Library.Template

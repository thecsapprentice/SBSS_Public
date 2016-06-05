include ../Makefile.Library.PreMatter

CXXFLAGS += -DFP_ARCH=$(FP_ARCH) -I$(INCLUDE_INSTALL_DIRECTORY)/PhysBAM/Public_Library -I$(INCLUDE_INSTALL_DIRECTORY) -std=gnu++0x 

LD_FLAGS += -L$(LIBRARY_INSTALL_DIRECTORY)
DEPENDENCIES = 
EXT_DEPENDENCIES = -lPhysBAM_NewStyle

Library_Src=$(shell find . -name "*.cpp") 
Library_Headers=$(shell find . -name "*.h") 

LIBRARY_NAME=libEngineInterface
INCLUDE_NAME=EngineInterface

include ../Makefile.Library.Template

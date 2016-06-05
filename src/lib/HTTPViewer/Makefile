include ../Makefile.Library.PreMatter

CXXFLAGS += -I$(INCLUDE_INSTALL_DIRECTORY) -std=gnu++0x

LD_FLAGS += -L$(LIBRARY_INSTALL_DIRECTORY)
DEPENDENCIES = 
EXT_DEPENDENCIES = -lcurl -lmongoose -lboost_system -lboost_filesystem -lboost_thread-mt

Library_Src=$(shell find . -name "*.cpp") 
Library_Headers=$(shell find . -name "*.h") 

LIBRARY_NAME=libHTTPViewer
INCLUDE_NAME=HTTPViewer

include ../Makefile.Library.Template

include ../Makefile.Library.PreMatter

CXXFLAGS += -I$(INCLUDE_INSTALL_DIRECTORY)/PhysBAM/Public_Library -Ilibtar -Ilibtar_listhash -I./ -Icompat
LD_FLAGS +=

Library_Src=$(shell find . -name "*.cpp")
Library_Headers=ARCHIVE.h

LIBRARY_NAME=libArchive
INCLUDE_NAME=Archive

include ../Makefile.Library.Template

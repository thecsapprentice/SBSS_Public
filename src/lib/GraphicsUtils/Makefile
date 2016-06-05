include ../Makefile.Library.PreMatter

CXXFLAGS+= -I$(INCLUDE_INSTALL_DIRECTORY) -Dlinux -std=gnu++0x 
LD_FLAGS+= -L$(LIBRARY_INSTALL_DIRECTORY) 
DEPENDENCIES = 
EXT_DEPENDENCIES = 

Library_Src = wxGraphics.cpp Bitmap.cpp lightsShaders.cpp lines.cpp GLmatrices.cpp staticTriangle.cpp glslTriangle.cpp shapes.cpp textures.cpp trackball.cpp skinGraphics.cpp trianglesUVW.cpp
Library_Headers = $(shell find . -name "*.h")

LIBRARY_NAME=libGraphicsUtils
INCLUDE_NAME=GraphicsUtils

include ../Makefile.Library.Template



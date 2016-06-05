all: install
.PHONY: install

INCLUDES = -I$(INCLUDE_INSTALL_DIRECTORY) -I$(INCLUDE_INSTALL_DIRECTORY)/GraphicsUtils -I$(INCLUDE_INSTALL_DIRECTORY)/json_spirit

ifeq ($(TYPE), debug)
CXXFLAGS+= -g -O0
else
CXXFLAGS+= -g -O3 
endif

ifdef USE_PHYSICS
CXXFLAGS +=  
else
CXXFLAGS += -DNO_PHYSICS
endif


ifdef USE_RPC_INTERFACE
  ENGINE_LIBS = -lCLERPC
  EXT_LIBS = -lthrift
  INCLUDES += -I$(INCLUDE_INSTALL_DIRECTORY)/CLERPC
  CXXFLAGS += -DUSE_RPC_INTERFACE
else
  ifdef USE_NEW_ENGINE
    ifdef USE_ENGINE_OPTIMIZATIONS
      ENGINE_LIBS = -lEngineFrontend -lCGKernels -lSIMD_Optimized_Kernels -lTQ -lCutter -lCommon 
      EXT_LIBS = -lPhysBAM_NewStyle -lboost_system  -lboost_thread-mt
    else
      ENGINE_LIBS = -lNE_Engine -lTQ
      EXT_LIBS = -lPhysBAM_NewStyle  -lboost_system  -lboost_thread-mt
      INCLUDES += -I$(INCLUDE_INSTALL_DIRECTORY)/NE_Engine
    endif
  else
    ENGINE_LIBS = -lCLE_Engine -lPhysBAM_OldStyle  -lTQ
    EXT_LIBS =  -lboost_system -lboost_thread-mt
    INCLUDES += -I$(INCLUDE_INSTALL_DIRECTORY)/CLE_Engine
  endif
endif

ENGINE_LIBS += -lGraphicsUtils -lHTTPViewer

ifdef USE_MPI
ENGINE_LIBS:= $(addsuffix -mpi, $(ENGINE_LIBS))
CXXFLAGS += -DENABLE_MPI
endif

CXXFLAGS+= -Dlinux -std=gnu++0x `wx-config --cxxflags` $(INCLUDES) -lpthread
LD_FLAGS+=-L$(LIBRARY_INSTALL_DIRECTORY) `wx-config --libs --gl-libs` -lGL -lGLU -lGLEW $(ENGINE_LIBS) $(EXT_LIBS)



LocalFlaps_Src = main.cpp cleScene.cpp surgicalActions.cpp hooks.cpp incision.cpp sutures.cpp fences.cpp config.cpp
LocalFlaps_Headers = cleScene.h surgicalActions.h hooks.h incision.h sutures.h fences.h config.h
LocalFlaps_Objs = $(addprefix $(BUILD_DIRECTORY)/SBSS/,$(patsubst %.cpp,%.o,$(LocalFlaps_Src)))
LocalFlaps_LibToolObjs = $(addprefix $(BUILD_DIRECTORY)/SBSS/,$(patsubst %.cpp,%.lo,$(LocalFlaps_Src)))

$(BUILD_DIRECTORY)/SBSS/%.lo: %.cpp $(LocalFlaps_Headers)
	@mkdir -p $(BUILD_DIRECTORY)/SBSS
	@mkdir -p $(dir $@)
	libtool --tag=$(LIBTOOL_TAG_FLAG) --mode=compile $(REALCXX) $(CXXFLAGS) -c -o $@ $< $(LIBTOOL_TYPE_FLAG)

$(BUILD_DIRECTORY)/SBSS/SBSS: $(LocalFlaps_LibToolObjs)
	libtool --tag=$(LIBTOOL_TAG_FLAG) --mode=link $(REALCXX) -o $@ $(LocalFlaps_LibToolObjs) $(LIBTOOL_TYPE_FLAG) -rpath $(LIBRARY_INSTALL_DIRECTORY) $(LD_FLAGS)

install: $(BUILD_DIRECTORY)/SBSS/SBSS
	@libtool --mode=install install -c $(BUILD_DIRECTORY)/SBSS/SBSS  $(EXEC_INSTALL_DIRECTORY)/SBSS.exe

clean:
	rm -rf $(BUILD_DIRECTORY)/SBSS
	rm -rf $(EXEC_INSTALL_DIRECTORY)/SBSS.exe

cmake_minimum_required(VERSION 2.8)

# Name of the project (will be the name of the plugin)
#project(addon)

# Build a shared library named after the project from the files in `src/`
add_library(CLEjs SHARED
 CLEjs.cpp
 )

# Gives our library file a .node extension without any "lib" prefix
set_target_properties(CLEjs PROPERTIES PREFIX "" SUFFIX ".node")

# Essential include files to build a node addon,
# You should add this line in every CMake.js based project
target_include_directories(CLEjs PRIVATE ${CMAKE_JS_INC})
target_include_directories(CLEjs PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/node_modules/nan/ )
#target_include_directories(CLEjs PRIVATE node_modules/nan/ )
target_include_directories(CLEjs PRIVATE /usr/include/node/ )
target_include_directories(CLEjs PRIVATE /nobackup.1/SOFTWARE/include/node/ )
target_include_directories(CLEjs PRIVATE ../../lib/ )


add_custom_target(CLEjs_node_modules ALL
                  COMMAND npm install nan
                  COMMENT "Checking out NaN for native Nodejs compilation..." )

add_dependencies(CLEjs CLEjs_node_modules)

# Essential library files to link to a node addon
# You should add this line in every CMake.js based project
target_link_libraries(CLEjs 
  HTTPViewer
  GraphicsUtils
  Thread_Queueing
  EngineFrontend
  EngineInterfaceLocal
  EngineInterface
  EngineBackend
  Common
  CG_Optimized_Kernels
  SIMD_Optimized_Kernels
  PhysBAM_Geometry
  PhysBAM_Tools
  PhysBAM_Common_Geometry
  PhysBAM_Nonmanifold_Implicit_Objects
  curl
  GL
  GLEW
  pthread
  boost_system
  boost_thread
  )

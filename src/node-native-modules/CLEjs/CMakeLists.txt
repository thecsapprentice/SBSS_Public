cmake_minimum_required(VERSION 2.8)

# Name of the project (will be the name of the plugin)
#project(addon)

# Build a shared library named after the project from the files in `src/`
add_library(CLEjs SHARED
 CLEjs.cpp
 )

# Gives our library file a .node extension without any "lib" prefix
set_target_properties(CLEjs PROPERTIES PREFIX "" SUFFIX ".node")

target_include_directories(CLEjs PRIVATE ../../lib/ )

add_dependencies(CLEjs native-modules-prebuild)

# Essential library files to link to a node addon
# You should add this line in every CMake.js based project
target_link_libraries(CLEjs 
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
  tet
  GL
  GLEW
  pthread
  )

install( TARGETS CLEjs
  LIBRARY DESTINATION share/sbss/node_modules )

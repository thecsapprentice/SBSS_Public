cmake_minimum_required(VERSION 3.0)

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
  )

install( TARGETS CLEjs
  LIBRARY DESTINATION share/sbss/node_modules )

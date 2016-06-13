cmake_minimum_required(VERSION 2.8)

# Build a shared library named after the project from the files in `src/`
add_library(legacy_cutter SHARED
  incision.cpp
  legacy-cutter.cpp
 )

# Gives our library file a .node extension without any "lib" prefix
set_target_properties(legacy_cutter PROPERTIES PREFIX "" SUFFIX ".node")

target_include_directories(legacy_cutter PRIVATE ../../lib/ )

add_dependencies(CLEjs native-modules-prebuild)
        

# Essential library files to link to a node addon
# You should add this line in every CMake.js based project
target_link_libraries(legacy_cutter 
  GraphicsUtils
  )

install( TARGETS legacy_cutter
  LIBRARY DESTINATION share/sbss/node_modules )
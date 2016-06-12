cmake_minimum_required(VERSION 2.8)

# Name of the project (will be the name of the plugin)
#project(addon)

# Build a shared library named after the project from the files in `src/`
add_library(legacy_cutter SHARED
  incision.cpp
  legacy-cutter.cpp
 )

# Gives our library file a .node extension without any "lib" prefix
set_target_properties(legacy_cutter PROPERTIES PREFIX "" SUFFIX ".node")

# Essential include files to build a node addon,
# You should add this line in every CMake.js based project
target_include_directories(legacy_cutter PRIVATE ${CMAKE_JS_INC})
target_include_directories(legacy_cutter PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/node_modules/nan/ )
#target_include_directories(legacy_cutter PRIVATE node_modules/nan/ )
target_include_directories(legacy_cutter PRIVATE /usr/include/node/ )
target_include_directories(legacy_cutter PRIVATE /nobackup.1/SOFTWARE/include/node/ )
target_include_directories(legacy_cutter PRIVATE ../../lib/ )


add_custom_target(legacy_cutter_node_modules ALL
                  COMMAND npm install nan
                  COMMENT "Checking out NaN for native Nodejs compilation..." )

add_dependencies(legacy_cutter legacy_cutter_node_modules)

# Essential library files to link to a node addon
# You should add this line in every CMake.js based project
target_link_libraries(legacy_cutter 
  GraphicsUtils
  )

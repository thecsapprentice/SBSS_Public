cmake_minimum_required(VERSION 3.2.0)
project (SBSS)

# Store all completed executables and libraries just off the main build tree
SET(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)
SET(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)

add_subdirectory( src )

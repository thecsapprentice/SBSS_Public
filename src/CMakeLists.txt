SET(CPP11_FLAGS "-std=gnu++11")
SET(OPENMP_FLAGS "-fopenmp")
SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${CPP11_FLAGS} ${OPENMP_FLAGS}" )

add_subdirectory( lib )
add_subdirectory( app )

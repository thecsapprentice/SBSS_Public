SET(CPP11_FLAGS "-std=gnu++11 -fPIC -g -O0")
SET(OPENMP_FLAGS "-fopenmp")
SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${CPP11_FLAGS} ${OPENMP_FLAGS}" )

set(electron_version "1.2.2")
set(node_archive_hash "d07ef2282d1e718432bc4467cf2d4861")
set(electron_archive_hash_linux "8d773eff0778659f920f7148e69c6833")
set(electron_archive_hash_apple "e50b0c2d8e56fe15104508861300afa7")
set(electron_archive_hash_win32 "6aae5ff8186c3f2bc05385168123183b")

set(nan_version "2.3.5")
set(nan_archive_hash "eb84043b0c5203dcb5693fdee7d4fb4e")


add_subdirectory( lib )
add_subdirectory( node-native-modules )
add_subdirectory( app )

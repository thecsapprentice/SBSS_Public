

set(electron_version "0.36.12")
set(nan_version "2.3.5")

SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -I${CMAKE_CURRENT_BINARY_DIR}/node-v${electron_version}/deps/uv/include/" )
SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -I${CMAKE_CURRENT_BINARY_DIR}/node-v${electron_version}/deps/v8/include/" )
SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -I${CMAKE_CURRENT_BINARY_DIR}/node-v${electron_version}/src/" )
SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -I${CMAKE_CURRENT_BINARY_DIR}/nan-${nan_version}/" )


add_custom_target(native-modules-prebuild ALL
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  DEPENDS
  ${CMAKE_CURRENT_BINARY_DIR}/nan-${nan_version}
  ${CMAKE_CURRENT_BINARY_DIR}/node-v${electron_version}
    )

# Download NaN headers for compiling native modules

set(nan_url "https://github.com/nodejs/nan/archive/v${nan_version}.tar.gz" )
set(nan_archive_path "${CMAKE_CURRENT_BINARY_DIR}/nan-v${nan_version}.tar.gz")
if (NOT EXISTS ${nan_archive_path} )
    file(DOWNLOAD "${nan_url}" ${nan_archive_path})
endif()

add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/nan-${nan_version}
  COMMAND  ${CMAKE_COMMAND} -E tar xzf ${CMAKE_CURRENT_BINARY_DIR}/nan-v${nan_version}.tar.gz
  DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/nan-v${nan_version}.tar.gz
  )

# Download Node headers for correct electron version

set(node_url "https://atom.io/download/atom-shell/v${electron_version}/node-v${electron_version}.tar.gz")
set(node_archive_path "${CMAKE_CURRENT_BINARY_DIR}/node-v${electron_version}.tar.gz")
if (NOT EXISTS ${node_archive_path} )
    file(DOWNLOAD "${node_url}" ${node_archive_path})
  endif()
  
add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/node-v${electron_version}
  COMMAND  ${CMAKE_COMMAND} -E tar xzf ${CMAKE_CURRENT_BINARY_DIR}/node-v${electron_version}.tar.gz
  DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/node-v${electron_version}.tar.gz
  )


add_subdirectory( CLEjs )
add_subdirectory( legacy-cutter )
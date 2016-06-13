message(STATUS "Using npm to install extra nodejs libraries.")

execute_process(
  COMMAND npm install
  WORKING_DIRECTORY $ENV{DESTDIR}/${CMAKE_INSTALL_PREFIX}/share/sbss/app
  OUTPUT_VARIABLE _output
  RESULT_VARIABLE retcode
  OUTPUT_QUIET
  )

message(STATUS "${_output}")

if(NOT "${retcode}" STREQUAL "0")
  message(FATAL_ERROR "Fatal error when using custom process.")
endif()


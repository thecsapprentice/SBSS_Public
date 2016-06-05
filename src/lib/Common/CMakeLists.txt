add_library( Common 
  KEYFRAMED_PARAMETER_ANIMATION.cpp
  Write_Output.cpp
  )

target_include_directories (Common PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../PhysBAMLegacy/)
target_include_directories (Common PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../ )

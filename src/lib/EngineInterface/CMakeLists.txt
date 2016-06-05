add_library( EngineInterface 
  ENGINE_INTERFACE.cpp
  )

target_include_directories (EngineInterface PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../PhysBAMLegacy/)
target_include_directories (EngineInterface PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../ )

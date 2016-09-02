add_library( EngineInterfaceLocal
  ENGINE_INTERFACE_LOCAL.cpp
  )

target_link_libraries( EngineInterfaceLocal
  PRIVATE
  EngineInterface
  EngineBackend
  )

target_include_directories (EngineInterfaceLocal PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../PhysBAMLegacy/)
target_include_directories (EngineInterfaceLocal PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../ )

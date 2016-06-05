add_library( Archive 
  ARCHIVE.cpp
  )

target_include_directories (Archive PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/libtar/)
target_include_directories (Archive PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/libtar_listhash/)
target_include_directories (Archive PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../PhysBAMLegacy/)
target_include_directories (Archive PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../ )

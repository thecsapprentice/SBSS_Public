add_library( HTTPViewer 
  DATA.cpp
  SERVER.cpp
  )

target_include_directories (HTTPViewer PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../PhysBAMLegacy/)
target_include_directories (HTTPViewer PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../ )

add_library( Thread_Queueing 
  PTHREAD_QUEUE.cpp
  )

target_include_directories (Thread_Queueing PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../PhysBAMLegacy/)
target_include_directories (Thread_Queueing PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../ )

if( ENABLE_MPI_OFFLOAD )

  add_library( EngineInterfaceMPI
    ENGINE_INTERFACE_MPI.cpp
    )
  
  target_compile_definitions( EngineInterfaceMPI
    PRIVATE
    ENABLE_MPI_OFFLOAD=1
    )

  find_package(MPI REQUIRED)

  target_link_libraries( EngineInterfaceMPI
    PRIVATE
    ${MPI_CXX_LIBRARIES}
    )

  if(MPI_CXX_COMPILE_FLAGS)
    set_target_properties(EngineInterfaceMPI PROPERTIES
      COMPILE_FLAGS "${MPI_CXX_COMPILE_FLAGS}")
  endif()
  
  if(MPI_CXX_LINK_FLAGS)
    set_target_properties(EngineInterfaceMPI PROPERTIES
      LINK_FLAGS "${MPI_CXX_LINK_FLAGS}")
  endif()

  target_include_directories (EngineInterfaceMPI PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../PhysBAMLegacy/)
  target_include_directories (EngineInterfaceMPI PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../ )
  target_include_directories (EngineInterfaceMPI PRIVATE ${MPI_CXX_INCLUDE_PATH} )
  
else( ENABLE_MPI_OFFLOAD )
  
endif( ENABLE_MPI_OFFLOAD )


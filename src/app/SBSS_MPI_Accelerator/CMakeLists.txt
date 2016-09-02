if( ENABLE_MPI_OFFLOAD )

  add_executable( SBSS_MPI_Accelerator
    CORE_ENGINE_ACCELERATOR.cpp
    main.cpp
    )
  
  target_link_libraries( SBSS_MPI_Accelerator
    PRIVATE
    EngineBackend
    EngineInterface
    EngineInterfaceLocal
    ${MPI_CXX_LIBRARIES}
    )
  
  if(MPI_CXX_COMPILE_FLAGS)
    set_target_properties(SBSS_MPI_Accelerator PROPERTIES
      COMPILE_FLAGS "${MPI_CXX_COMPILE_FLAGS}")
  endif()
  
  if(MPI_CXX_LINK_FLAGS)
    set_target_properties(SBSS_MPI_Accelerator PROPERTIES
      LINK_FLAGS "${MPI_CXX_LINK_FLAGS}")
  endif()

  target_include_directories (SBSS_MPI_Accelerator PRIVATE ${MPI_CXX_INCLUDE_PATH} )
  target_include_directories (SBSS_MPI_Accelerator PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../../lib/)

endif( ENABLE_MPI_OFFLOAD )

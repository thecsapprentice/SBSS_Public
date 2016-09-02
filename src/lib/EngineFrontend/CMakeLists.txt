
set(ENABLE_SELF_COLLISIONS OFF CACHE BOOL "Enable self-collision support. Requires Tetgen")
set(ENABLE_MPI_OFFLOAD OFF CACHE BOOL "Enable MPI offloading support to local Phi Cards")

add_library( EngineFrontend
ANALYTIC_COLLISION.cpp
CLElib.cpp
ELASTIC_LATTICE_DEFORMER.cpp
EMBEDDING_TOOLS.cpp
OPENVDB_LEVELSET_COLLISION_2D.cpp
OPENVDB_LEVELSET_COLLISION_3D.cpp
PHYSBAM_LEVELSET_COLLISION.cpp
PHYSBAM_RIGID_BODY_LEVELSET_COLLISION.cpp
Write_Output.cpp
)

if( ENABLE_SELF_COLLISIONS )
  target_compile_definitions( EngineFrontend
    PRIVATE
    ENABLE_SELF_COLLISIONS=1
    )

  target_link_libraries( EngineFrontend
    PRIVATE
    tet
    )
endif( ENABLE_SELF_COLLISIONS )

if( ENABLE_MPI_OFFLOAD )
  target_compile_definitions( EngineFrontend
    PRIVATE
    ENABLE_MPI_OFFLOAD=1
    )

  target_link_libraries( EngineFrontend
    PRIVATE
    EngineInterfaceMPI      
    )
else( ENABLE_MPI_OFFLOAD )
  target_link_libraries( EngineFrontend
    PRIVATE
    EngineInterfaceLocal    
    )
endif( ENABLE_MPI_OFFLOAD )

target_link_libraries( EngineFrontend
  PRIVATE
  EngineInterface
  Common
  PhysBAM_Geometry
  PhysBAM_Tools
  PhysBAM_Common_Geometry
  PhysBAM_Nonmanifold_Implicit_Objects
  )

target_include_directories (EngineFrontend PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../PhysBAMLegacy/)
target_include_directories (EngineFrontend PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../ )

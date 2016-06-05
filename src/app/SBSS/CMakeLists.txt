add_executable( SBSS 
cleScene.cpp
config.cpp
fences.cpp
#flaptoolsdlg.cpp
hooks.cpp
incision.cpp
#localFlapsViewer.cpp
main.cpp
surgicalActions.cpp
sutures.cpp)

target_include_directories (SBSS PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../../lib/)
target_link_libraries( SBSS 
  HTTPViewer
  GraphicsUtils
  Thread_Queueing
  EngineFrontend
  EngineInterfaceLocal
  EngineInterface
  EngineBackend
  Common
  CG_Optimized_Kernels
  SIMD_Optimized_Kernels
  PhysBAM_Geometry
  PhysBAM_Tools
  PhysBAM_Common_Geometry
  PhysBAM_Nonmanifold_Implicit_Objects
  curl
  GL
  GLEW
  pthread
  boost_system
  boost_thread
  )

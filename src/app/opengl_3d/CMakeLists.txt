add_executable( opengl_3d
main.cpp)

target_include_directories (opengl_3d PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../../lib/)
target_link_libraries( opengl_3d 
  PhysBAM_Rendering
  PhysBAM_Geometry
  PhysBAM_Tools
  GL
  GLEW
  X11
  glut
  GLU
  )

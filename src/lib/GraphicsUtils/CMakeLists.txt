add_library( GraphicsUtils 
Bitmap.cpp
#bmesh.cpp
#disk_cycle.cpp
#element_cycle.cpp
#elements.cpp
GLmatrices.cpp
glslTriangle.cpp
lightsShaders.cpp
lines.cpp
#loop_cycle.cpp
#radial_cycle.cpp
shapes.cpp
skinGraphics.cpp
staticTriangle.cpp
textures.cpp
trackball.cpp
trianglesUVW.cpp
wxGraphics.cpp

  )

target_link_libraries( GraphicsUtils
  PRIVATE
  GL
  GLEW
  )


target_include_directories (GraphicsUtils PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../PhysBAMLegacy/)
target_include_directories (GraphicsUtils PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../ )

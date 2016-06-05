add_library( CG_Optimized_Kernels 
  Vector_SAXPY_V9/Vector_SAXPY_V9_Helper.cpp
  Project/Project_Helper.cpp
  Vector_Subtract/Vector_Subtract_Helper.cpp
  Vector_SAXPY_V8/Vector_SAXPY_V8_Helper.cpp
  Vector_SAXPY_A/Vector_SAXPY_A_Helper.cpp
  Vector_Times/Vector_Times_Helper.cpp
  Vector_SAXPY_AB/Vector_SAXPY_AB_Helper.cpp
  Convergence_Norm/Convergence_Norm_Helper.cpp
  Vector_SAXPY_V/Vector_SAXPY_V_Helper.cpp
  Dot_Product/Dot_Product_Helper.cpp
  Vector_SAXPY_V6/Vector_SAXPY_V6_Helper.cpp
  Vector_Add/Vector_Add_Helper.cpp
  Vector_Set/Vector_Set_Helper.cpp
  )

target_include_directories (CG_Optimized_Kernels PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../PhysBAMLegacy/)
target_include_directories (CG_Optimized_Kernels PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../ )

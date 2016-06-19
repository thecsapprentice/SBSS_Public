
add_custom_target(electron-prebuild ALL
 WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
 DEPENDS
 ${CMAKE_CURRENT_BINARY_DIR}/electron-v${electron_version}/electron
 )

#Download Electron headers for compiling native modules
if(UNIX AND NOT APPLE)
  set(electron_url "https://github.com/electron/electron/releases/download/v${electron_version}/electron-v${electron_version}-linux-x64.zip" )
  set(electron_archive_hash ${electron_archive_hash_linux})
endif()
if(APPLE)
  set(electron_url "https://github.com/electron/electron/releases/download/v${electron_version}/electron-v${electron_version}-darwin-x64.zip" )
  set(electron_archive_hash ${electron_archive_hash_apple})
endif()
if(WIN32)
  set(electron_url "https://github.com/electron/electron/releases/download/v${electron_version}/electron-v${electron_version}-win32-x64.zip" )
  set(electron_archive_hash ${electron_archive_hash_win32})
endif()

set(electron_archive_path "${CMAKE_CURRENT_BINARY_DIR}/electron-v${electron_version}/electron-v${electron_version}.zip")
if (NOT EXISTS ${electron_archive_path} )
  file(DOWNLOAD "${electron_url}" ${electron_archive_path}
    EXPECT_MD5 ${electron_archive_hash}
    )
endif()

add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/electron-v${electron_version}/electron
  COMMAND  ${CMAKE_COMMAND} -E tar xf ${CMAKE_CURRENT_BINARY_DIR}/electron-v${electron_version}/electron-v${electron_version}.zip
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/electron-v${electron_version}
  DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/electron-v${electron_version}/electron-v${electron_version}.zip
  )

#add_subdirectory( SBSS )

install(DIRECTORY electron-SBSS/  DESTINATION share/sbss/app
  PATTERN "node_modules" EXCLUDE
  REGEX "~$" EXCLUDE
  )

install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/electron-v${electron_version}/  DESTINATION share/sbss/
  USE_SOURCE_PERMISSIONS
  PATTERN "node_modules" EXCLUDE
  PATTERN "electron-v${electron_version}.zip" EXCLUDE
  REGEX "~$" EXCLUDE
  )

install(SCRIPT "./install_node_packages.cmake")

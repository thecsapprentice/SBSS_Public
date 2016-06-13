add_subdirectory( SBSS )

install(DIRECTORY electron-SBSS/  DESTINATION share/sbss
  PATTERN "node_modules" EXCLUDE
  REGEX "~$" EXCLUDE
  )
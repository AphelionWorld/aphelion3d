#----------------------------------------------------------------
# Generated CMake target import file for configuration "Debug".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "Aphelion3D::libmath" for configuration "Debug"
set_property(TARGET Aphelion3D::libmath APPEND PROPERTY IMPORTED_CONFIGURATIONS DEBUG)
set_target_properties(Aphelion3D::libmath PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_DEBUG "CXX"
  IMPORTED_LOCATION_DEBUG "${_IMPORT_PREFIX}/lib/libaphelion_math.a"
  )

list(APPEND _cmake_import_check_targets Aphelion3D::libmath )
list(APPEND _cmake_import_check_files_for_Aphelion3D::libmath "${_IMPORT_PREFIX}/lib/libaphelion_math.a" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)

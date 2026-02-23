#----------------------------------------------------------------
# Generated CMake target import file for configuration "Release".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "CLime::lime" for configuration "Release"
set_property(TARGET CLime::lime APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(CLime::lime PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "C"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/liblime.a"
  )

list(APPEND _cmake_import_check_targets CLime::lime )
list(APPEND _cmake_import_check_files_for_CLime::lime "${_IMPORT_PREFIX}/lib/liblime.a" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)

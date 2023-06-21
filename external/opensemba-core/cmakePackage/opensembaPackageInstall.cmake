
set(LIBNAME opensemba)
message(STATUS "Configuring package installation for ${LIBNAME}")

set("OPENSEMBA_VERSION" "0.16")
message(STATUS "OPENSEMBA_VERSION="${OPENSEMBA_VERSION})

# Relative to CMAKE_INSTALL_PREFIX:
set(ConfigPackageLocation "${LIBNAME}/")

# target_include_directories(libexample PUBLIC ${CMAKE_INSTALL_PREFIX})
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/${LIBNAME}ConfigVersion.cmake"
  VERSION       ${OPENSEMBA_VERSION}
  COMPATIBILITY AnyNewerVersion
)

# Copy the FooConfig.cmake to the build/Foo directory:
configure_file(
  cmakePackage/${LIBNAME}Config.cmake
  "${CMAKE_CURRENT_BINARY_DIR}/${LIBNAME}Config.cmake"
  COPYONLY
)

# --- Installation ---
install(
   TARGETS opensemba 
    opensemba_core 
        opensemba_core_geometry 
        opensemba_core_math 
        opensemba_core_outputrequest
        opensemba_core_physicalmodel 
        opensemba_core_source
    opensemba_parsers opensemba_exporters
   EXPORT ${LIBNAME}Targets
   DESTINATION ${ConfigPackageLocation}/lib/
)

# Export the targets (change namespace appropriately):
export(
  EXPORT    ${LIBNAME}Targets
  FILE      "${CMAKE_CURRENT_BINARY_DIR}/${LIBNAME}Targets.cmake"
)

# This also installs relative to CMAKE_INSTALL_PREFIX:
install(
  FILES       cmakePackage/${LIBNAME}Config.cmake
              "${CMAKE_CURRENT_BINARY_DIR}/${LIBNAME}ConfigVersion.cmake"
  DESTINATION ${ConfigPackageLocation}/share/
)

# This is for the source files, change location appropriately:
install(
  DIRECTORY "src/"
  DESTINATION ${ConfigPackageLocation}/include/opensemba/
  FILES_MATCHING PATTERN "*.h"
)

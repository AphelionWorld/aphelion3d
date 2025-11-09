cmake_minimum_required(VERSION 3.21)
project(Aphelion3D LANGUAGES C CXX VERSION 0.1.0)

# Put all build outputs under build/.../bin and build/.../lib
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

# Also cover multi-config generators just in case (VS, Ninja Multi-Config)
foreach(OUTPUTCONFIG Debug Release RelWithDebInfo MinSizeRel)
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_${OUTPUTCONFIG} ${CMAKE_BINARY_DIR}/bin)
  set(CMAKE_LIBRARY_OUTPUT_DIRECTORY_${OUTPUTCONFIG} ${CMAKE_BINARY_DIR}/lib)
  set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY_${OUTPUTCONFIG} ${CMAKE_BINARY_DIR}/lib)
endforeach()

set(CMAKE_COLOR_DIAGNOSTICS ON)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

option(BUILD_SHARED_LIBS "Build shared instead of static libs" OFF)
option(BUILD_TESTING "Build tests" ON)
option(BUILD_DOCS "Build documentation" ON)

include(GNUInstallDirs)
include(CMakeDependentOption)

set(APHELION3D_EXPORT_NAME "Aphelion3DTargets")
set(APHELION3D_NAMESPACE "Aphelion3D::")

# Find system libraries
find_package(OpenMP REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(pmp REQUIRED)
find_package(GDAL REQUIRED)
find_package(Boost REQUIRED COMPONENTS system filesystem graph)
find_package(PDAL REQUIRED)
find_package(nlohmann_json REQUIRED)

# Use earcut from submodule (header-only)
# Add external libs (vendor)
set(TINYGLTF_HEADER_ONLY OFF CACHE BOOL "Use compiled tinygltf library" FORCE)
set(TINYGLTF_INSTALL OFF CACHE BOOL "Do not install tinygltf from vendored subdir" FORCE)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/src/external/tinygltf)

add_subdirectory(src/libutil)
add_subdirectory(src/libmath)
add_subdirectory(src/core)
add_subdirectory(src/app/hello)
add_subdirectory(src/app/viewer)
add_subdirectory(src/app/mesher)  # Add this line

# Add poly2tri from submodule
set(P2T_BUILD_TESTS OFF CACHE BOOL "Disable poly2tri tests" FORCE)
set(P2T_BUILD_TESTBED OFF CACHE BOOL "Disable poly2tri testbed" FORCE)
add_subdirectory(src/external/poly2tri)

if(BUILD_TESTING)
  include(CTest)
  enable_testing()
  add_subdirectory(tests)
endif()

# Documentation
if(BUILD_DOCS)
  find_package(Doxygen)
  if(DOXYGEN_FOUND)
    configure_file(
      ${CMAKE_CURRENT_SOURCE_DIR}/docs/Doxyfile.in
      ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
      @ONLY
    )
    add_custom_target(docs
      COMMAND ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
      COMMENT "Generating API documentation with Doxygen"
      VERBATIM
    )
    message(STATUS "Doxygen found: documentation can be built with 'make docs' or 'ninja docs'")
  else()
    message(STATUS "Doxygen not found: documentation will not be built")
  endif()
endif()

include(CMakePackageConfigHelpers)

install(DIRECTORY src/libmath/include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/aphelion3d)

write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/Aphelion3DConfigVersion.cmake"
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY AnyNewerVersion
)

configure_package_config_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/cmake/Aphelion3DConfig.cmake.in"
  "${CMAKE_CURRENT_BINARY_DIR}/Aphelion3DConfig.cmake"
  INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/Aphelion3D
)

install(FILES
  "${CMAKE_CURRENT_BINARY_DIR}/Aphelion3DConfig.cmake"
  "${CMAKE_CURRENT_BINARY_DIR}/Aphelion3DConfigVersion.cmake"
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/Aphelion3D
)

install(EXPORT ${APHELION3D_EXPORT_NAME}
  NAMESPACE ${APHELION3D_NAMESPACE}
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/Aphelion3D
)

export(EXPORT ${APHELION3D_EXPORT_NAME}
  NAMESPACE ${APHELION3D_NAMESPACE}
  FILE "${CMAKE_CURRENT_BINARY_DIR}/Aphelion3DTargets.cmake"
)

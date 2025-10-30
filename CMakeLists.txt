cmake_minimum_required(VERSION 3.21)
project(Aphelion3D LANGUAGES C CXX VERSION 0.1.0)

set(CMAKE_COLOR_DIAGNOSTICS ON)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

option(BUILD_SHARED_LIBS "Build shared instead of static libs" OFF)
include(GNUInstallDirs)

include(CMakeDependentOption)
option(BUILD_TESTING "Build tests" ON)

set(APHELION3D_EXPORT_NAME "Aphelion3DTargets")
set(APHELION3D_NAMESPACE "Aphelion3D::")

# Find system libraries
find_package(OpenMP REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(pmp REQUIRED)

add_subdirectory(src/libutil)
add_subdirectory(src/libmath)
add_subdirectory(src/core)
add_subdirectory(src/app/hello)
add_subdirectory(src/app/viewer)

if(BUILD_TESTING)
  include(CTest)
  enable_testing()
  add_subdirectory(tests)
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

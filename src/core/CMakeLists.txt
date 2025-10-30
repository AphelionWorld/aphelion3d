file(GLOB_RECURSE CORE_HEADERS "${CMAKE_CURRENT_SOURCE_DIR}/include/*.hpp")
add_library(core SHARED ${CORE_HEADERS} src/mesh.cpp)
add_library(Aphelion3D::core ALIAS core)

target_include_directories(core PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>
)

target_link_libraries(core PUBLIC 
  Aphelion3D::libutil 
  Aphelion3D::libmath
  pmp
  Eigen3::Eigen
)

target_compile_definitions(core PRIVATE CORE_SHARED)

set_target_properties(core PROPERTIES OUTPUT_NAME "aphelion_core")

install(TARGETS core
  EXPORT ${APHELION3D_EXPORT_NAME}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)

install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/aphelion3d)
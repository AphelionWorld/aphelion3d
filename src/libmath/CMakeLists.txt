file(GLOB_RECURSE LIBMATH_HEADERS "${CMAKE_CURRENT_SOURCE_DIR}/include/*.hpp")
add_library(libmath ${LIBMATH_HEADERS} src/math.cpp)
add_library(Aphelion3D::libmath ALIAS libmath)

target_include_directories(libmath PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>
)

target_link_libraries(libmath PUBLIC Aphelion3D::libutil)

target_compile_definitions(libmath PRIVATE $<$<BOOL:${BUILD_SHARED_LIBS}>:LIBMATH_SHARED>)

set_target_properties(libmath PROPERTIES OUTPUT_NAME "aphelion_math")

install(TARGETS libmath
  EXPORT ${APHELION3D_EXPORT_NAME}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)

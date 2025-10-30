add_library(libutil INTERFACE)
add_library(Aphelion3D::libutil ALIAS libutil)

target_include_directories(libutil INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>
)

install(TARGETS libutil
  EXPORT ${APHELION3D_EXPORT_NAME}
)

install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/aphelion3d)

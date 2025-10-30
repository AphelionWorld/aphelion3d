add_executable(viewer main.cpp)
target_link_libraries(viewer PRIVATE Aphelion3D::libmath Aphelion3D::libutil)
target_compile_definitions(viewer PRIVATE APP_NAME="viewer")

install(TARGETS viewer
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

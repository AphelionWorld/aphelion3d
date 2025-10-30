add_executable(hello main.cpp)
target_link_libraries(hello PRIVATE Aphelion3D::libmath Aphelion3D::libutil)

install(TARGETS hello
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

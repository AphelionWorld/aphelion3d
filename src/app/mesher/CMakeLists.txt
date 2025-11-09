add_executable(mesher ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp)

target_link_libraries(mesher PRIVATE 
    Aphelion3D::core
    Aphelion3D::libutil
    Aphelion3D::libmath
)

target_compile_definitions(mesher PRIVATE APP_NAME="mesher")

install(TARGETS mesher
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)
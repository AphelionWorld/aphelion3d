include(FetchContent)
FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/refs/tags/v1.14.0.zip
)
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

add_executable(math_tests math_tests.cpp)
target_link_libraries(math_tests PRIVATE Aphelion3D::libmath gtest_main)

add_test(NAME math_tests COMMAND math_tests)

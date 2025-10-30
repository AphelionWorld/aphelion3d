# CMake generated Testfile for 
# Source directory: D:/AphelionWorld/Dev/aphelion3d/tests
# Build directory: D:/AphelionWorld/Dev/aphelion3d/build/ucrt64-debug/tests
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test([=[math_tests]=] "D:/AphelionWorld/Dev/aphelion3d/build/ucrt64-debug/tests/math_tests.exe")
set_tests_properties([=[math_tests]=] PROPERTIES  _BACKTRACE_TRIPLES "D:/AphelionWorld/Dev/aphelion3d/tests/CMakeLists.txt;12;add_test;D:/AphelionWorld/Dev/aphelion3d/tests/CMakeLists.txt;0;")
subdirs("../_deps/googletest-build")

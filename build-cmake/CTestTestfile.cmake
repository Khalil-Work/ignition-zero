# CMake generated Testfile for 
# Source directory: C:/GitHub-Codex/original-projects/ignition-zero
# Build directory: C:/GitHub-Codex/original-projects/ignition-zero/build-cmake
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
if(CTEST_CONFIGURATION_TYPE MATCHES "^([Dd][Ee][Bb][Uu][Gg])$")
  add_test(ignition_checks "C:/GitHub-Codex/original-projects/ignition-zero/build-cmake/Debug/test_ignition.exe")
  set_tests_properties(ignition_checks PROPERTIES  _BACKTRACE_TRIPLES "C:/GitHub-Codex/original-projects/ignition-zero/CMakeLists.txt;28;add_test;C:/GitHub-Codex/original-projects/ignition-zero/CMakeLists.txt;0;")
elseif(CTEST_CONFIGURATION_TYPE MATCHES "^([Rr][Ee][Ll][Ee][Aa][Ss][Ee])$")
  add_test(ignition_checks "C:/GitHub-Codex/original-projects/ignition-zero/build-cmake/Release/test_ignition.exe")
  set_tests_properties(ignition_checks PROPERTIES  _BACKTRACE_TRIPLES "C:/GitHub-Codex/original-projects/ignition-zero/CMakeLists.txt;28;add_test;C:/GitHub-Codex/original-projects/ignition-zero/CMakeLists.txt;0;")
elseif(CTEST_CONFIGURATION_TYPE MATCHES "^([Mm][Ii][Nn][Ss][Ii][Zz][Ee][Rr][Ee][Ll])$")
  add_test(ignition_checks "C:/GitHub-Codex/original-projects/ignition-zero/build-cmake/MinSizeRel/test_ignition.exe")
  set_tests_properties(ignition_checks PROPERTIES  _BACKTRACE_TRIPLES "C:/GitHub-Codex/original-projects/ignition-zero/CMakeLists.txt;28;add_test;C:/GitHub-Codex/original-projects/ignition-zero/CMakeLists.txt;0;")
elseif(CTEST_CONFIGURATION_TYPE MATCHES "^([Rr][Ee][Ll][Ww][Ii][Tt][Hh][Dd][Ee][Bb][Ii][Nn][Ff][Oo])$")
  add_test(ignition_checks "C:/GitHub-Codex/original-projects/ignition-zero/build-cmake/RelWithDebInfo/test_ignition.exe")
  set_tests_properties(ignition_checks PROPERTIES  _BACKTRACE_TRIPLES "C:/GitHub-Codex/original-projects/ignition-zero/CMakeLists.txt;28;add_test;C:/GitHub-Codex/original-projects/ignition-zero/CMakeLists.txt;0;")
else()
  add_test(ignition_checks NOT_AVAILABLE)
endif()

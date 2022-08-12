cmake_minimum_required(VERSION 3.14)
project(cMchecker)

# GoogleTest requires at least C++14
set(CMAKE_CXX_STANDARD 14)

# Bring the headers into the project
include_directories(include)
# Add src files individualy 
set(CLI_SOURCES src/cMchecker.cpp)
# Or add all src files
# file(GLOB SOURCES "src/*.cpp")
add_executable(cMchecker ${CLI_SOURCES})

include(FetchContent)
FetchContent_Declare(
  googletest
  GIT_REPOSITORY https://github.com/google/googletest.git
  GIT_TAG release-1.12.1
)
# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)


enable_testing()

add_executable(
    run_tests
    test/run_tests.cpp
)

target_link_libraries(
    run_tests
    GTest::gtest_main
)

include(GoogleTest)
gtest_discover_tests(run_tests)

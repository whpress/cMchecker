cmake_minimum_required(VERSION 3.10)
project(cMchecker)

# Bring the headers into the project
include_directories(include)

# Add src files individualy 
# set(SOURCES src/cMchecker.cpp)
# Or add all src files
file(GLOB SOURCES "src/*.cpp")

add_executable(cMchecker ${SOURCES})
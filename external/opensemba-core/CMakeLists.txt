cmake_minimum_required(VERSION 3.0)

set(CMAKE_CXX_STANDARD 17)

if(NOT DEFINED CMAKE_TOOLCHAIN_FILE AND DEFINED ENV{CMAKE_TOOLCHAIN_FILE})
    set(CMAKE_TOOLCHAIN_FILE $ENV{CMAKE_TOOLCHAIN_FILE})
endif()

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

project(opensemba CXX)

message(STATUS "Creating build system for opensemba with the following options:")
option(OPENSEMBA_ENABLE_TESTS OFF)

message(STATUS "OPEMSEMBA_ENABLE_TESTS: " ${OPENSEMBA_ENABLE_TESTS})

include_directories(src/)

add_subdirectory(src/core)
add_subdirectory(src/parsers)
add_subdirectory(src/exporters)

add_library(opensemba "src/opensemba.cpp")
target_link_libraries(opensemba opensemba_core opensemba_parsers opensemba_exporters)

include(cmakePackage/opensembaPackageInstall.cmake)

if(OPENSEMBA_ENABLE_TESTS)
	enable_testing()
    add_subdirectory(test/)
endif()
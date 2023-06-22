cmake_minimum_required(VERSION 3.0)

message(STATUS "Creating build system for ${PROJECT_NAME}_Tests")

add_executable(${PROJECT_NAME}_tests
 "ElectrostaticSolverTest.cpp"
 "mfemTest.cpp"
)

include_directories(./)

find_package(GTest CONFIG REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS})

target_link_libraries(${PROJECT_NAME}_tests 
	pulmtln-solver
	GTest::gtest
	GTest::gtest_main
)


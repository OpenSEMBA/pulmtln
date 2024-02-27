cmake_minimum_required(VERSION 3.21)

message(STATUS "Creating build system for ${PROJECT_NAME}_tests")

find_package(GTest CONFIG REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS})

add_executable(${PROJECT_NAME}_tests
	"ElectrostaticSolverTest.cpp"
	"mfemTest.cpp"
	"DriverTest.cpp"
    "ParserTest.cpp"
	"ParametersTest.cpp"
	"ModelTest.cpp"
	"DirectedGraphTest.cpp"
	"DomainTest.cpp")

include_directories(./)
 
target_link_libraries(${PROJECT_NAME}_tests 
	pulmtln-core
	GTest::gtest
	GTest::gtest_main
)


cmake_minimum_required(VERSION 3.0)

message(STATUS "Creating build system for ${PROJECT_NAME}_Tests")

add_executable(${PROJECT_NAME}_tests
	"pulmtlnTest.cpp"
)

include_directories(./)

target_link_libraries(${PROJECT_NAME}_tests 
	pulmtln-solver
	GTest::gtest
)


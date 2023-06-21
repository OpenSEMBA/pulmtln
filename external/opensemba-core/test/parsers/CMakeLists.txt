cmake_minimum_required(VERSION 3.10)

message(STATUS "Configuring build system for opensemba_parsers_tests.")

find_package(GTest CONFIG REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS})
include(GoogleTest)

add_executable(opensemba_parsers_tests 
    "json/ParserTest.cpp"
    "stl/ParserTest.cpp"
)

target_link_libraries(opensemba_parsers_tests 
    opensemba_parsers
    GTest::gtest     
    GTest::gtest_main
)


gtest_add_tests(
    TARGET opensemba_parsers_tests
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
)
gtest_discover_tests(
    opensemba_parsers_tests
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
)
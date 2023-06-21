cmake_minimum_required(VERSION 3.0)

project(opensemba_parsers CXX)

add_library(opensemba_parsers
    "Parser.cpp"
    "json/Parser.cpp"
    "stl/Parser.cpp"
)

target_link_libraries(opensemba_parsers opensemba_core)

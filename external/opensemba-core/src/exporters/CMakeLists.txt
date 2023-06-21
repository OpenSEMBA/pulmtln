cmake_minimum_required(VERSION 3.0)

project(opensemba_exporters CXX)

add_library(opensemba_exporters
    "Exporter.cpp"
    "ExporterVTK.cpp"
)
target_link_libraries(opensemba_exporters opensemba_core)
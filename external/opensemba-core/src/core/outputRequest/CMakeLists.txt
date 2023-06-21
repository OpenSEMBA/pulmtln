cmake_minimum_required(VERSION 3.0)

set(CMAKE_CXX_STANDARD 17)

project(opensemba_core_outputrequest CXX)

add_library(opensemba_core_outputrequest STATIC
    "BulkCurrent.cpp"
    "Domain.cpp"
    "FarField.cpp"
    "OnLayer.cpp"
    "OnLine.cpp"
    "OnPoint.cpp"
    "OnSurface.cpp"
    "OutputRequest.cpp"
)
target_link_libraries(opensemba_core_outputrequest opensemba_core_geometry)

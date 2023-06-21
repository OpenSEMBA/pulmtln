cmake_minimum_required(VERSION 3.0)

set(CMAKE_CXX_STANDARD 17)

project(opensemba_core_source CXX)

add_library(opensemba_core_source
    "Source.cpp"    
    "PlaneWave.cpp"
    "OnLine.cpp"
    "Generator.cpp"
    "Dipole.cpp"
    
    "port/Port.cpp"
    "port/TEM.cpp"
    "port/TEMCoaxial.cpp"
    "port/Waveguide.cpp"
    "port/WaveguideRectangular.cpp"

    "magnitude/Magnitude.cpp"
    "magnitude/Numerical.cpp"
)
target_link_libraries(opensemba_core_source opensemba_core_geometry)

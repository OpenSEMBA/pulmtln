cmake_minimum_required(VERSION 3.0)

set(CMAKE_CXX_STANDARD 17)

project(opensemba_core_geometry CXX)

add_library(opensemba_core_geometry STATIC
    "coordinate/Conformal.cpp"
    "coordinate/Coordinate.cpp"
    "coordinate/Relative.cpp"
    "element/Element.cpp"
    "element/LineConformal.cpp"
    "element/Tetrahedron.cpp"
    "element/Tetrahedron4.cpp"
    "element/Triangle3.cpp"
    "Layer.cpp"
    "mesh/Geometric.cpp"
    "mesh/Structured.cpp"
    "mesh/Unstructured.cpp"
)


target_link_libraries(opensemba_core_geometry opensemba_core_math)

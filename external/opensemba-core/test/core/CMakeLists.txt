cmake_minimum_required(VERSION 3.10)

message(STATUS "Configuring build system for opensemba_core_tests.")

find_package(GTest CONFIG REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS})
include(GoogleTest)

add_executable(opensemba_core_tests
    "util/GroupIdentifiableUniqueTest.cpp"
    "util/ProjectTest.cpp"
    "geometry/coordinate/CoordinateGroupTest.cpp"
    "geometry/element/Line2Test.cpp"
    "geometry/layer/LayerTest.cpp"
    "geometry/mesh/GeometricTest.cpp"
    "geometry/mesh/StructuredTest.cpp"
    "geometry/mesh/UnstructuredTest.cpp"
    "geometry/BoxTest.cpp"
    "geometry/GridTest.cpp"
    "math/axis/LocalTest.cpp"
    "math/function/GaussianTest.cpp"
    "math/function/LinearInterpolationTest.cpp"
    "math/matrix/DynamicTest.cpp"
    "math/matrix/StaticTest.cpp"
    "math/simplex/LineTest.cpp"
    "math/simplex/TetrahedronTest.cpp"
    "math/simplex/TriangleTest.cpp"
    "math/util/RealTest.cpp"
    "math/vector/CartesianVectorTest.cpp"
    "model/ModelTest.cpp"
    "physicalModel/multiport/RLCTest.cpp"
    "physicalModel/volume/ClassicTest.cpp"
    "physicalModel/wire/ExtremesTest.cpp"
    "physicalModel/wire/WireTest.cpp"
    "physicalModel/BoundTest.cpp"
    "physicalModel/GroupTest.cpp"
    "source/magnitude/NumericalTest.cpp"
    "source/port/TEMCoaxialTest.cpp"
    "source/port/WaveguideRectangularTest.cpp"
    "ProblemDescriptionTest.cpp"
)

target_link_libraries(opensemba_core_tests 
    opensemba_core
    GTest::gtest 
    GTest::gtest_main
)

gtest_add_tests(
    TARGET opensemba_core_tests
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
)
gtest_discover_tests(
    opensemba_core_tests
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
)
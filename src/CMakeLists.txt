cmake_minimum_required(VERSION 3.21)

set(LIBRARY_NAME "pulmtln")
message(STATUS "Creating build system for ${LIBRARY_NAME}")

find_package(nlohmann_json CONFIG REQUIRED)

find_package(MFEM CONFIG REQUIRED)
include_directories(${MFEM_INCLUDE_DIRS})

find_package(Boost REQUIRED graph)
include_directories(${Boost_INCLUDE_DIRS})

add_library(${LIBRARY_NAME} STATIC 
	"FES.cpp"
	"ElectrostaticSolver.cpp"
	"Driver.cpp" 
	"Parser.cpp"
	"Model.cpp"
	"Parameters.cpp" 
	"DirectedGraph.cpp"
	"Domain.cpp")

target_link_libraries(${LIBRARY_NAME} mfem nlohmann_json Boost::graph)
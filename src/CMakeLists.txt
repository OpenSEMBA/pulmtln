cmake_minimum_required(VERSION 3.0)

set(LIBRARY_NAME "pulmtln-solver")
message(STATUS "Creating build system for ${LIBRARY_NAME}")

find_package(MFEM CONFIG REQUIRED)
include_directories(${MFEM_INCLUDE_DIRS})

find_package(nlohmann_json CONFIG REQUIRED)

add_library(${LIBRARY_NAME} STATIC 
	"FES.cpp"
	"ElectrostaticSolver.cpp"
	"Driver.cpp" 
 "Materials.h")

target_link_libraries(${LIBRARY_NAME} mfem nlohmann_json)
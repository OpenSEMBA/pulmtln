cmake_minimum_required(VERSION 3.0)

set(LIBRARY_NAME "pulmtln")
message(STATUS "Creating build system for ${LIBRARY_NAME}")

find_package(nlohmann_json CONFIG REQUIRED)

find_package(MFEM CONFIG REQUIRED)
include_directories(${MFEM_INCLUDE_DIRS})

add_library(${LIBRARY_NAME} STATIC 
	"FES.cpp"
	"ElectrostaticSolver.cpp"
	"Driver.cpp" 
	"Parser.cpp"
	
)

target_link_libraries(${LIBRARY_NAME} mfem nlohmann_json)
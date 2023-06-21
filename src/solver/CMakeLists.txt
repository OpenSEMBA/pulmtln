cmake_minimum_required(VERSION 3.0)

set(LIBRARY_NAME "pulmtln-solver")
message(STATUS "Creating build system for ${LIBRARY_NAME}")

find_package(MFEM CONFIG REQUIRED)
include_directories(${MFEM_INCLUDE_DIRS})

add_library(${LIBRARY_NAME} STATIC 
	"Material.cpp"
	"Model.cpp"
	"FES.cpp"
	"Solver.cpp"
)

get_filename_component(PARENT_DIR ../ ABSOLUTE)
include_directories(${PARENT_DIR})

target_link_libraries(${LIBRARY_NAME} mfem)

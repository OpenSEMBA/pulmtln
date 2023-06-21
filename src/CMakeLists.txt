cmake_minimum_required(VERSION 3.0)


find_package(MFEM CONFIG REQUIRED)
include_directories(${MFEM_INCLUDE_DIRS})

#add_subdirectory(adapter)
add_subdirectory(solver)

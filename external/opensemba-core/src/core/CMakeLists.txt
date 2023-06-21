cmake_minimum_required(VERSION 3.0)

project(opensemba_core CXX)

include_directories(${CMAKE_CURRENT_LIST_DIR})

add_subdirectory(model/ model/)
add_subdirectory(math/ math/)
add_subdirectory(geometry/ geometry/)
add_subdirectory(outputRequest/ outputRequest/)
add_subdirectory(physicalModel/ physicalModel/)
add_subdirectory(source/ source/)

find_package(nlohmann_json CONFIG REQUIRED)

add_library(opensemba_core STATIC
    "ProblemDescription.cpp"
)
    
target_link_libraries(opensemba_core
    opensemba_core_geometry
    opensemba_core_math
    opensemba_core_physicalmodel
    opensemba_core_source
    opensemba_core_outputrequest
    nlohmann_json
)

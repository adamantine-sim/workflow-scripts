cmake_minimum_required(VERSION 3.13.4)
project(snap_tessellated_mesh)

find_package(deal.II 9.5.0 REQUIRED HINTS ${DEAL_II_DIR})

if(NOT DEAL_II_WITH_OPENCASCADE)
  message(FATAL_ERROR "This tool requires deal.II to be built with OpenCASCADE support!")
endif()

add_executable(snap_to_tessellated_mesh.exe snap_to_tessellated_mesh.cpp)

target_link_libraries(snap_to_tessellated_mesh.exe dealii::dealii)

enable_testing()
add_subdirectory(tests)

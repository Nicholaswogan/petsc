cmake_minimum_required(VERSION 3.14)
project(TSPSEUDO LANGUAGES Fortran)

find_package(LAPACK REQUIRED)

add_library(tspseudo pseudo_transient.f90)
target_link_libraries(tspseudo PRIVATE ${LAPACK_LIBRARIES})

add_executable(robertson_ptc_demo robertson_ptc_demo.f90)
target_link_libraries(robertson_ptc_demo PRIVATE tspseudo)

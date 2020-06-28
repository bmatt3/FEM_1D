#Minimum Version
cmake_minimum_required(VERSION 2.8)

#Project Name
project(FEM)

set(CMAKE_CXX_COMPILER gcc)

set(CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR}/bin)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})
set(LIBARY_OUTPUT_PATH ${CMAKE_BINARY_DIR})

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/Modules")
include("${CMAKE_MODULE_PATH}/FindLAPACKE.cmake")
#include("${CMAKE_MODULE_PATH}/FindCBLAS.cmake")
include("${CMAKE_MODULE_PATH}/FindLAPACKEXT.cmake")

find_package(BLAS REQUIRED)
#if(BLAS_FOUND)
#  message(STATUS "BLAS library was found.")
#else()
#  message(STATUS "BLAS library was not found.")
#endif()

find_package(LAPACK REQUIRED)
if(LAPACK_FOUND)
  message(STATUS "LAPACK library was found.")
else()
  message(STATUS "LAPACK library was not found.")
endif()

find_package(LAPACKE REQUIRED)
if(LAPACKE_FOUND)
   message(STATUS "LAPACKE library was found.")
else()
   message(STATUS "LAPACKE library was not found.")
endif()
  
set(CMAKE_EXE_LINKER_FLAGS "-O3 -fopt-info-vec-all")
set(SOURCES src/matrix_utils.c inc/matrix_utils.h src/fem.c)
add_executable(FEM ${SOURCES})
target_link_libraries(FEM -g -O3 -llapacke -llapack -lblas -lm -Wall)

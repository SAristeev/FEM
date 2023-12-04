cmake_minimum_required(VERSION 3.17.0 FATAL_ERROR)
project(FEM VERSION 0.3.0 LANGUAGES CXX C)
set_property(GLOBAL PROPERTY USE_FOLDERS ON)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/CMake")

find_package(MKL REQUIRED)

add_executable(FEM
  FEM.cpp
  FEM.hpp
  fp.hpp
  main.cpp
  external/json/json.hpp
  external/base64/base64.h 
  external/base64/base64.cpp
)

target_include_directories(FEM PUBLIC ${MKL_INCLUDE_DIR})
target_link_libraries(FEM LINK_PUBLIC ${MKL_LIBRARIES})
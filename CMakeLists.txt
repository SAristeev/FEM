cmake_minimum_required(VERSION 3.17.0 FATAL_ERROR)
project(FEM VERSION 0.3.0 LANGUAGES CXX C)
set_property(GLOBAL PROPERTY USE_FOLDERS ON)

find_package(MKL CONFIG REQUIRED PATHS $ENV{MKLROOT})

add_executable(FEM
  FEM.cpp
  FEM.hpp
  fp.hpp
  main.cpp
  external/json/json.hpp
  external/base64/base64.h 
  external/base64/base64.cpp
)

target_link_libraries(FEM PUBLIC MKL::MKL)
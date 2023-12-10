cmake_minimum_required(VERSION 3.17 FATAL_ERROR)
project(vtu_writer CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_library(vtu_writer 
  src/vtu_writer.cpp
)

target_include_directories(vtu_writer 
  PUBLIC include 
  PRIVATE src
  PRIVATE external
)

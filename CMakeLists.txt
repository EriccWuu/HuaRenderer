cmake_minimum_required(VERSION 3.20)

set(ROOT_DIR "${CMAKE_CURRENT_SOURCE_DIR}")
set(RENDERER_ROOT_DIR "${ROOT_DIR}/renderer")
set(THIRD_PARTY_ROOT "${RENDERER_ROOT_DIR}/3rdparty")
set(SDL2_ROOT "${RENDERER_ROOT_DIR}/3rdparty/SDL2-2.30.4")

set(renderer_name HuaRenderer)

project(${renderer_name})

set(EXPORT_COMPILE_COMMANDS ON)

# 设置 C++ 标准
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")

include(cmake/CopyFiles.cmake)

add_subdirectory(renderer)

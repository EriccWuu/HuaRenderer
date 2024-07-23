set(THIRD_PARITY_DIRS "3rdparty")
set(THIRD_PARITY_INCLUDE_DIRS "include")

add_subdirectory(3rdparty)

# 主程序入口
add_executable(${renderer_name} src/main.cpp)

# 添加 include 目录
include_directories(${THIRD_PARITY_INCLUDE_DIRS})
include_directories(src)

# sdl2
include(3rdparty/SDL2-2.30.4/cmake/sdl2-config.cmake)
target_include_directories(${renderer_name} PUBLIC ${SDL2_INCLUDE_DIRS})
target_link_libraries(${renderer_name} ${SDL2_LIBRARIES})

# 链接库
target_link_libraries(${renderer_name} CoreLayer RenderingLayer PresentationLayer)

add_subdirectory(src/CoreLayer)
add_subdirectory(src/RenderingLayer)
add_subdirectory(src/PresentationLayer)
add_subdirectory(test)

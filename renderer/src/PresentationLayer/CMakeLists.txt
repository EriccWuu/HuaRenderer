# PresentationLayer CMakeLists.txt
file(GLOB PRESENTATION_LAYER_SOURCES "*.cpp" "*.h")
add_library(PresentationLayer ${PRESENTATION_LAYER_SOURCES})

# # 链接 SDL2 库
# find_package(SDL2 REQUIRED)
# target_include_directories(PresentationLayer PUBLIC ${SDL2_INCLUDE_DIRS})
# target_link_libraries(PresentationLayer PUBLIC ${SDL2_LIBRARIES})

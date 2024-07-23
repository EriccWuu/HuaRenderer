# src/RenderingLayer/CMakeLists.txt
file(GLOB_RECURSE SOFTWARE_RASTERIZER_SOURCES "SoftwareRasterizer/*.cpp" "SoftwareRasterizer/*.h")
# file(GLOB SOFTWARE_RAYTRACER_SOURCES "SoftwareRayTracer/*.cpp" "SoftwareRayTracer/*.h")
file(GLOB VULKAN_RASTERIZER_SOURCES "VulkanRasterizer/*.cpp" "VulkanRasterizer/*.h")

# Create library
add_library(RenderingLayer ${SOFTWARE_RASTERIZER_SOURCES} ${SOFTWARE_RAYTRACER_SOURCES} ${VULKAN_RASTERIZER_SOURCES})

# set include directiories
target_include_directories(RenderingLayer PRIVATE ${renderer_name})

# link libraries
target_link_libraries(RenderingLayer PRIVATE CoreLayer)
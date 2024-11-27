cmake_minimum_required(VERSION 3.29)
project(Raytracer)

set(CMAKE_CXX_STANDARD 20)

find_package(OpenGL REQUIRED)
find_package(OpenMP REQUIRED)

set(GLFW_BUILD_DOCS OFF CACHE BOOL "" FORCE)
set(GLFW_BUILD_TESTS OFF CACHE BOOL "" FORCE)
set(GLFW_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)

set(JSON_BuildTests OFF CACHE INTERNAL "")
set(JSON_Install OFF CACHE INTERNAL "")

add_subdirectory(lib/glfw)
add_subdirectory(lib/json)

# setup ImGui as a static library
add_library(imgui STATIC
        lib/imgui/imconfig.h
        lib/imgui/imgui.cpp
        lib/imgui/imgui.h
        lib/imgui/imgui_demo.cpp
        lib/imgui/imgui_draw.cpp
        lib/imgui/imgui_internal.h
        lib/imgui/imgui_tables.cpp
        lib/imgui/imgui_widgets.cpp
        lib/imgui/imstb_rectpack.h
        lib/imgui/imstb_textedit.h
        lib/imgui/imstb_truetype.h
        lib/imgui/backends/imgui_impl_glfw.cpp
        lib/imgui/backends/imgui_impl_glfw.h
        lib/imgui/backends/imgui_impl_opengl3.cpp
        lib/imgui/backends/imgui_impl_opengl3.h
        lib/imgui/backends/imgui_impl_opengl3_loader.h)
target_include_directories(imgui PUBLIC
        lib/imgui
        lib/imgui/backends
)
target_link_libraries(imgui PRIVATE glfw)

add_executable(Raytracer
    raytracer.cpp
    raytracer.h
    viewport.cpp
    viewport.h
    camera.cpp
    camera.h
    ray.cpp
    ray.h
        geometry.cpp
    geometry.h
        primitive.cpp
        primitive.h
        scene.cpp
        scene.h
        light.cpp
        light.h
        parser.cpp
        parser.h
        material.cpp
        material.h
        image.cpp
        image.h
)
target_link_libraries(Raytracer PRIVATE OpenGL::GL glfw)
target_link_libraries(Raytracer PRIVATE nlohmann_json::nlohmann_json)
target_link_libraries(Raytracer PRIVATE imgui)

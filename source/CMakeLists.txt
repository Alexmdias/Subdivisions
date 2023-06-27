cmake_minimum_required(VERSION 3.16)
project(example)

list(PREPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# libigl
option(LIBIGL_GLFW "Build target igl::glfw" ON)
include(libigl)

set(
        main_src
        main.cpp
        trimesh.h
        trimesh.cpp
        trimesh_types.h
        loop_scheme.h
        loop_scheme.cpp
        butterfly_scheme.h
        butterfly_scheme.cpp
        sqrt3_scheme.h
        sqrt3_scheme.cpp)

# Add your project files
add_executable(${PROJECT_NAME} ${main_src} )
target_link_libraries(${PROJECT_NAME} PUBLIC igl::glfw)
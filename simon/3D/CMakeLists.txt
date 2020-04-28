cmake_minimum_required(VERSION 3.7)
set(CMAKE_CXX_STANDARD 11)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake/sdl2)
set(SOURCE_FILES libs/sdw/Colour.cpp libs/sdw/Utils.cpp libs/sdw/DrawingWindow.cpp libs/sdw/ModelTriangle.cpp Model.cpp Material.h Light.h Transformable.cpp Camera.cpp Animatable.cpp Movement.cpp Times.cpp VectorUtil.cpp Rigidbody.cpp Updatable.hpp Temporal.cpp Vertex.cpp Triangle.cpp Clipping.cpp Lines.cpp SSAA.cpp RedNoise.cpp)

IF(WIN32)
    set(CMAKE_EXE_LINKER_FLAGS -mconsole)
ENDIF(WIN32)

project(RedNoise)

find_package(SDL2 REQUIRED)
find_package(OpenMP)
include_directories(${SDL2_INCLUDE_DIRS} libs/sdw libs/glm)

add_executable(diagnostic ${SOURCE_FILES})
set_target_properties(diagnostic PROPERTIES OUTPUT_NAME ${CMAKE_PROJECT_NAME})
target_compile_options(diagnostic PUBLIC -Wall -Werror -pedantic -O1 -fsanitize=undefined -fsanitize=address -fno-omit-frame-pointer)
target_link_libraries(diagnostic SDL2::Main -fsanitize=undefined -fsanitize=address -fno-omit-frame-pointer OpenMP::OpenMP_CXX)

add_executable(production EXCLUDE_FROM_ALL ${SOURCE_FILES})
set_target_properties(production PROPERTIES OUTPUT_NAME ${CMAKE_PROJECT_NAME})
target_link_libraries(production SDL2::Main OpenMP::OpenMP_CXX)

add_executable(debug EXCLUDE_FROM_ALL ${SOURCE_FILES})
set_target_properties(debug PROPERTIES OUTPUT_NAME ${CMAKE_PROJECT_NAME})
target_compile_options(debug PUBLIC -Wall -ggdb -g3)
target_link_libraries(debug SDL2::Main OpenMP::OpenMP_CXX)

add_executable(speedy EXCLUDE_FROM_ALL ${SOURCE_FILES})
set_target_properties(speedy PROPERTIES OUTPUT_NAME ${CMAKE_PROJECT_NAME})
target_compile_options(speedy PUBLIC -Wall -Ofast -funsafe-math-optimizations -fno-finite-math-only -march=native)
target_link_libraries(speedy SDL2::Main OpenMP::OpenMP_CXX)
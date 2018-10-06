cmake_minimum_required(VERSION 3.1)
project(photo_fit VERSION 1.0.0 LANGUAGES CXX)

file(GLOB sources ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${CMAKE_CURRENT_SOURCE_DIR}/src/*.h)

add_executable(${PROJECT_NAME} ${sources} ${headers})

target_include_directories(${PROJECT_NAME} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src )

target_compile_options(${PROJECT_NAME} PRIVATE -Wall -pipe -O3 -std=c++11)

find_package(Eigen3 3.3 REQUIRED NO_MODULE)
target_link_libraries(${PROJECT_NAME} Eigen3::Eigen ecf)

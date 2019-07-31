cmake_minimum_required(VERSION 3.8)
project(photo_fit VERSION 1.1.0 LANGUAGES CXX)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release")
endif()

if (CMAKE_BUILD_TYPE EQUAL "DEBUG")
     message("debug mode")
endif()


add_executable(${PROJECT_NAME}
                src/ecf_extension/EvalOp.h
                src/ecf_extension/EvalOp.cpp
                src/ecf_extension/GaussMutOp.h
                src/ecf_extension/GaussMutOp.cpp
                src/ecf_extension/MyFloatingPoint.h
                src/ecf_extension/MyFloatingPoint.cpp
                src/gaussian_functor.h
                src/gaussian_functor.cpp
                src/basis.h
                src/basis.cpp
                src/control_data.h
                src/control_data.cpp
                src/utils.h
                src/utils.cpp
                src/main.cpp
                src/procedures.h
                src/procedures.cpp
                )

target_compile_features(${PROJECT_NAME} PUBLIC
                        cxx_std_11)

target_compile_definitions(${PROJECT_NAME} PUBLIC "$<$<CONFIG:DEBUG>:FIT_DEBUG>")

target_include_directories(${PROJECT_NAME} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src )

target_compile_options(${PROJECT_NAME} PRIVATE -Wall)

find_package(OpenMP REQUIRED)
find_package (Eigen3 3.3 REQUIRED NO_MODULE)

find_library(ECF NAMES libecf.so)
 
target_link_libraries (${PROJECT_NAME} PRIVATE OpenMP::OpenMP_CXX Eigen3::Eigen ${ECF})
 
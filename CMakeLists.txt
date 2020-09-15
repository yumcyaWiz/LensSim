cmake_minimum_required(VERSION 3.12)
project(LensSim LANGUAGES C CXX)

# build type
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to Release as none was specified.")
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build." FORCE)
endif()


# pthread
find_package(Threads REQUIRED)

# LensSim
add_library(LensSim)
# needed for pybind11 binding
set_target_properties(LensSim PROPERTIES POSITION_INDEPENDENT_CODE ON)

# compile settings
target_compile_features(LensSim PUBLIC cxx_std_17)
set_target_properties(LensSim PROPERTIES CXX_STANDARD_REQUIRED ON)
set_target_properties(LensSim PROPERTIES CXX_EXTENSIONS OFF)

# src
target_sources(LensSim PRIVATE 
  src/core/spectrum.cpp 
  src/parallel/parallel.cpp 
  src/samplers/rng.cpp 
  src/lens-system/lens-system.cpp 
  src/film.cpp
)

# ext
add_subdirectory(ext)

# LensSim - include and link
target_include_directories(LensSim PUBLIC src)
target_link_libraries(LensSim PUBLIC Threads::Threads)
target_link_libraries(LensSim PUBLIC pool)
target_link_libraries(LensSim PUBLIC nlohmann_json::nlohmann_json)
target_link_libraries(LensSim PUBLIC stb)
target_link_libraries(LensSim PUBLIC tinyexr)

# IBL Example
add_executable(main src/main.cpp)
target_include_directories(main PRIVATE LensSim)
target_link_libraries(main PRIVATE LensSim)

# pybind11
pybind11_add_module(binding python/binding.cpp)
set_target_properties(binding PROPERTIES OUTPUT_NAME "LensSim")
target_include_directories(binding PRIVATE LensSim)
target_link_libraries(binding PRIVATE LensSim)

# copy build output to python/
add_custom_command(TARGET binding POST_BUILD 
  COMMAND "${CMAKE_COMMAND}" -E copy "$<TARGET_FILE:binding>" "${PROJECT_SOURCE_DIR}/python/$<TARGET_FILE_NAME:binding>" 
  COMMENT "Copying build output to python/"
)
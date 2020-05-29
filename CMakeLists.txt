cmake_minimum_required(VERSION 3.13..3.15)
project(LensSym LANGUAGES C CXX)

# build type
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to Release as none was specified.")
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build." FORCE)
endif()


# pthread
find_package(Threads REQUIRED)

# executable
add_executable(main src/sphere.cpp)

# compile settings
target_compile_features(main PRIVATE cxx_std_17)
set_target_properties(main PROPERTIES CXX_EXTENSIONS OFF)

# src
target_sources(main PRIVATE src/core/spectrum.cpp src/parallel/parallel.cpp src/samplers/rng.cpp)

# ext
add_subdirectory(ext)

# include and link
target_include_directories(main PRIVATE src)
target_link_libraries(main PRIVATE Threads::Threads)
target_link_libraries(main PRIVATE pool)
target_link_libraries(main PRIVATE nlohmann_json::nlohmann_json)
target_link_libraries(main PRIVATE stb)
target_link_libraries(main PRIVATE tinyexr)
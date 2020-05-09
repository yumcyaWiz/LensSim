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
add_executable(main src/main.cpp)

# compile settings
target_compile_features(main PUBLIC cxx_std_17)
set_target_properties(main PROPERTIES CXX_EXTENSIONS OFF)

# src
target_sources(main PUBLIC src/parallel.cpp)

# ext
add_subdirectory(ext)

# include and link
target_include_directories(main PUBLIC src)
target_link_libraries(main PUBLIC Threads::Threads)
target_link_libraries(main PUBLIC pool)

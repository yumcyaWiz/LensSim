cmake_minimum_required(VERSION 3.13..3.15)
project(LensSym LANGUAGES C CXX)

# build type
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to Release as none was specified.")
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build." FORCE)
endif()

# pthread
find_package(Threads REQUIRED)

add_executable(main src/main.cpp)
target_sources(main PUBLIC src/parallel.cpp)

add_subdirectory(ext)

target_include_directories(main PUBLIC src)
target_link_libraries(main PUBLIC Threads::Threads)
target_link_libraries(main PUBLIC pool)

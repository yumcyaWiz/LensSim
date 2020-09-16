add_executable(test test.cpp)
target_include_directories(test PRIVATE LensSim)
target_link_libraries(test PRIVATE LensSim)
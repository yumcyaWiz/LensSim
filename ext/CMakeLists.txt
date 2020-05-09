# ThreadPool
add_library(pool INTERFACE)
target_include_directories(pool SYSTEM INTERFACE
  ThreadPool
)

# json
set(JSON_BuildTests OFF CACHE INTERNAL "")
add_subdirectory(json)

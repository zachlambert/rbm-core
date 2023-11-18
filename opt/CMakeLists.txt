
add_library(mathbox-opt INTERFACE)
target_include_directories(mathbox-opt
    INTERFACE include
)
target_link_libraries(mathbox-opt
    INTERFACE mathbox-types mathbox-diff
)

add_executable(example_opt example/main.cpp)
target_link_libraries(example_opt mathbox-opt)

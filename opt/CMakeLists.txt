
add_library(owl_opt INTERFACE)
target_include_directories(owl_opt
    INTERFACE include
)
target_link_libraries(owl_opt
    INTERFACE owl_types owl_diff
)

add_executable(example_opt example/main.cpp)
target_link_libraries(example_opt owl_opt)


add_library(owl_probability INTERFACE)
target_include_directories(owl_probability
    INTERFACE include
)
target_link_libraries(owl_probability
    INTERFACE owl_types cpp-utils owl_diff owl_transform
)

add_executable(example_probability example/main.cpp)
target_link_libraries(example_probability owl_probability)

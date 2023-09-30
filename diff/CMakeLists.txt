
add_library(owl_diff INTERFACE)
target_include_directories(owl_diff
    INTERFACE include
)
target_link_libraries(owl_diff
    INTERFACE owl_types owl_transform
)

add_executable(example_diff_dfunction example/dfunction.cpp)
target_link_libraries(example_diff_dfunction owl_diff)

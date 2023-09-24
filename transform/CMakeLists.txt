
add_library(math_transform INTERFACE)
target_include_directories(math_transform
    INTERFACE include
)
target_link_libraries(math_transform
    INTERFACE math_types core_serialize
)

add_executable(math_transform_test test/main.cpp)
target_link_libraries(math_transform_test PRIVATE math_transform)

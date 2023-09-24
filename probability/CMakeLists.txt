
add_library(math_probability INTERFACE)
target_include_directories(math_probability
    INTERFACE include
)
target_link_libraries(math_probability
    INTERFACE math_types core_std math_algebra math_transform
)

add_executable(math_probability_test test/main.cpp)
target_link_libraries(math_probability_test math_probability)

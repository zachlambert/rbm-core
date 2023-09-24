
add_library(math_algebra INTERFACE)
target_include_directories(math_algebra
    INTERFACE include
)
target_link_libraries(math_algebra
    INTERFACE math_types math_transform
)

add_executable(math_algebra_test_dfunction test/dfunction.cpp)
target_link_libraries(math_algebra_test_dfunction math_algebra)

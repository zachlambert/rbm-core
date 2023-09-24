
add_library(math_opt INTERFACE)
target_include_directories(math_opt
    INTERFACE include
)
target_link_libraries(math_opt
    INTERFACE math_types math_algebra
)

add_executable(math_opt_test test/main.cpp)
target_link_libraries(math_opt_test math_opt gui)

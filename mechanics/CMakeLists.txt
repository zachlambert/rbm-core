
add_library(math_mechanics INTERFACE)
target_include_directories(math_mechanics
    INTERFACE include
)
target_link_libraries(math_mechanics
    INTERFACE math_types math_transform math_geometry
)

add_executable(math_mechanics_test test/main.cpp)
target_link_libraries(math_mechanics_test math_mechanics)

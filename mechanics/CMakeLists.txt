
add_library(owl_mechanics INTERFACE)
target_include_directories(owl_mechanics
    INTERFACE include
)
target_link_libraries(owl_mechanics
    INTERFACE owl_types owl_transform owl_geometry
)

add_executable(example_mechanics example/main.cpp)
target_link_libraries(example_mechanics owl_mechanics)


add_library(owl_geometry SHARED
    src/height_map_grid.cpp
    src/height_map_conversions.cpp
    src/height_map_functional.cpp
)
target_include_directories(owl_geometry
    PUBLIC include
)
target_link_libraries(owl_geometry
    PUBLIC cpp-utils owl_types owl_transform
)

add_executable(example_geometry_collision_primitives example/collision_primitives.cpp)
target_link_libraries(example_geometry_collision_primitives PRIVATE owl_geometry)

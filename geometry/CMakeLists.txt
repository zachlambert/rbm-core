
add_library(mathbox-geometry SHARED
    src/height_map_grid.cpp
    src/height_map_conversions.cpp
    src/height_map_functional.cpp
)
target_include_directories(mathbox-geometry
    PUBLIC include
)
target_link_libraries(mathbox-geometry
    PUBLIC cpp-utils mathbox-types mathbox-transform
)

add_executable(example_geometry_collision_primitives example/collision_primitives.cpp)
target_link_libraries(example_geometry_collision_primitives PRIVATE mathbox-geometry)

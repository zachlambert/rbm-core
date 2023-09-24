
add_library(math_geometry SHARED
    src/height_map_grid.cpp
    src/height_map_conversions.cpp
    src/height_map_functional.cpp
)
target_include_directories(math_geometry
    PUBLIC include
)
target_link_libraries(math_geometry
    PUBLIC core_std math
)

add_executable(math_geometry_test_collision test/collision.cpp)
target_link_libraries(math_geometry_test_collision PRIVATE math_geometry)

add_executable(math_geometry_test_height_map test/height_map.cpp)
target_link_libraries(math_geometry_test_height_map PRIVATE math_geometry)

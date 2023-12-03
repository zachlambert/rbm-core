add_library(rbm-geometry SHARED
    src/height_map_grid.cpp
    src/height_map_conversions.cpp
    src/height_map_functional.cpp
)
target_include_directories(rbm-geometry PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)
target_link_libraries(rbm-geometry
    PUBLIC cpp-utils rbm-types rbm-transform
)

if (${BUILD_EXAMPLES})
    add_executable(example_geometry_collision_primitives example/collision_primitives.cpp)
    target_link_libraries(example_geometry_collision_primitives PRIVATE rbm-geometry)
endif()

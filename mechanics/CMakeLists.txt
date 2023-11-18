
add_library(mathbox-mechanics INTERFACE)
target_include_directories(mathbox-mechanics
    INTERFACE include
)
target_link_libraries(mathbox-mechanics
    INTERFACE mathbox-types mathbox-transform mathbox-geometry
)

add_executable(example_mechanics example/main.cpp)
target_link_libraries(example_mechanics mathbox-mechanics)

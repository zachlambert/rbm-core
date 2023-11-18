
add_library(mathbox-probability INTERFACE)
target_include_directories(mathbox-probability
    INTERFACE include
)
target_link_libraries(mathbox-probability
    INTERFACE mathbox-types cpp-utils mathbox-diff mathbox-transform
)

add_executable(example_probability example/main.cpp)
target_link_libraries(example_probability mathbox-probability)

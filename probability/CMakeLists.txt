
add_library(mbox-probability INTERFACE)
target_include_directories(mbox-probability
    INTERFACE include
)
target_link_libraries(mbox-probability
    INTERFACE mbox-types cpp-utils mbox-diff mbox-transform
)

add_executable(example_probability example/main.cpp)
target_link_libraries(example_probability mbox-probability)

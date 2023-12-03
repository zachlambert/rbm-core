
add_library(rbm-probability INTERFACE)
target_include_directories(rbm-probability INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)
target_link_libraries(rbm-probability
    INTERFACE rbm-types cpp-utils rbm-diff rbm-transform
)

add_executable(example_probability example/main.cpp)
target_link_libraries(example_probability rbm-probability)

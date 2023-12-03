
add_library(rbm-mechanics INTERFACE)
target_include_directories(rbm-mechanics INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)
target_link_libraries(rbm-mechanics
    INTERFACE rbm-types rbm-transform rbm-geometry
)

if (${BUILD_EXAMPLES})
    add_executable(example_mechanics example/main.cpp)
    target_link_libraries(example_mechanics rbm-mechanics)
endif()

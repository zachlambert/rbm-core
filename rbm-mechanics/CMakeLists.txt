
add_library(mbox-mechanics INTERFACE)
target_include_directories(mbox-mechanics INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)
target_link_libraries(mbox-mechanics
    INTERFACE mbox-types mbox-transform mbox-geometry
)

add_executable(example_mechanics example/main.cpp)
target_link_libraries(example_mechanics mbox-mechanics)

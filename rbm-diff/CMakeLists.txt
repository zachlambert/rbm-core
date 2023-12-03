
add_library(rbm-diff INTERFACE)
target_include_directories(rbm-diff INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)
target_link_libraries(rbm-diff
    INTERFACE rbm-types rbm-transform
)

add_executable(example_diff_dfunction example/dfunction.cpp)
target_link_libraries(example_diff_dfunction rbm-diff)

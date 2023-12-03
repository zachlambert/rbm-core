
add_library(rbm-opt INTERFACE)
target_include_directories(rbm-opt INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)
target_link_libraries(rbm-opt
    INTERFACE rbm-types rbm-diff
)

if(${USE_SVIZ})
    add_executable(example_opt example/main.cpp)
    target_link_libraries(example_opt rbm-opt sviz)
endif()

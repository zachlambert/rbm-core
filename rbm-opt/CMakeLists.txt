
add_library(mbox-opt INTERFACE)
target_include_directories(mbox-opt INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)
target_link_libraries(mbox-opt
    INTERFACE mbox-types mbox-diff
)

if(${USE_SVIZ})
    add_executable(example_opt example/main.cpp)
    target_link_libraries(example_opt mbox-opt sviz)
endif()

add_library(rbm-gui SHARED
    src/camera.cpp
    src/camera_controller.cpp
    src/font.cpp
    src/grid_renderer.cpp
    src/marker_renderer.cpp
    src/mesh_renderer.cpp
    src/shader.cpp
    src/viewport.cpp
    src/window.cpp
)
target_include_directories(rbm-gui PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)
target_link_libraries(rbm-gui
    PUBLIC rbm-cpp rbm-types rbm-transform rbm-geometry
)

if (${BUILD_EXAMPLES})
    add_executable(example_gui_window example/window.cpp)
    target_link_libraries(example_gui_window PRIVATE rbm-gui)
endif()

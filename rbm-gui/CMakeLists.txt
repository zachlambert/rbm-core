set(shader_list)
function(rbm_gui_add_shader shader_list parent_include file label)
    add_custom_command(
        OUTPUT ${CMAKE_CURRENT_LIST_DIR}/shaders/include/shader/${label}.h
        COMMAND ./gen_shader_h.sh ${file} ${label} ${parent_include}
        DEPENDS ${CMAKE_CURRENT_LIST_DIR}/shaders/${file} ${CMAKE_CURRENT_LIST_DIR}/gen_shader_h.sh
        WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}
    )
    set(${shader_list} ${${shader_list}} shaders/include/shader/${label}.h PARENT_SCOPE)
endfunction()

set(rbm-gui-shader-files "")
rbm_gui_add_shader(rbm-gui-shader-files rbm/gui grid.fs grid_fs)
rbm_gui_add_shader(rbm-gui-shader-files rbm/gui grid.vs grid_vs)
rbm_gui_add_shader(rbm-gui-shader-files rbm/gui marker.fs marker_fs)
rbm_gui_add_shader(rbm-gui-shader-files rbm/gui marker.vs marker_vs)
rbm_gui_add_shader(rbm-gui-shader-files rbm/gui mesh.fs mesh_fs)
rbm_gui_add_shader(rbm-gui-shader-files rbm/gui mesh.vs mesh_vs)
add_custom_target(rbm-gui-shaders ALL DEPENDS ${rbm-gui-shader-files})

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
    lib/imgui/imgui.cpp
    lib/imgui/imgui_demo.cpp
    lib/imgui/imgui_draw.cpp
    lib/imgui/imgui_tables.cpp
    lib/imgui/imgui_widgets.cpp
    lib/imgui/backends/imgui_impl_opengl3.cpp
    lib/imgui/backends/imgui_impl_glfw.cpp
)
target_include_directories(rbm-gui
    PRIVATE
        shaders/include
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/lib/imgui>
        $<INSTALL_INTERFACE:include>
)
target_link_libraries(rbm-gui
    PUBLIC rbm-cpp rbm-types rbm-transform rbm-geometry GL glfw GLEW ${CMAKE_DL_LIBS}
)
add_dependencies(rbm-gui rbm-gui-shaders)

if (${BUILD_EXAMPLES})
    add_executable(example_gui_window example/window.cpp)
    target_link_libraries(example_gui_window PRIVATE rbm-gui)
endif()

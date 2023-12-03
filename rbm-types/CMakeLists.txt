add_library(rbm-types INTERFACE)
target_include_directories(rbm-types INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)
target_link_libraries(rbm-types
    INTERFACE Eigen3::Eigen cpp-utils datapack
)

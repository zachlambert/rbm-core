add_library(mbox-types INTERFACE)
target_include_directories(mbox-types INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)
target_link_libraries(mbox-types
    INTERFACE Eigen3::Eigen cpp-utils datapack
)

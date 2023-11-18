add_library(mbox-types INTERFACE)
target_include_directories(mbox-types
    INTERFACE include
)
target_link_libraries(mbox-types
    INTERFACE Eigen3::Eigen cpp-utils datapack
)

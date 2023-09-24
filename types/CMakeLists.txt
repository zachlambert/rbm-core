add_library(owl_types INTERFACE)
target_include_directories(owl_types
    INTERFACE include
)
target_link_libraries(owl_types
    INTERFACE Eigen3::Eigen cpp-utils parrot
)

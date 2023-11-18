add_library(mathbox-types INTERFACE)
target_include_directories(mathbox-types
    INTERFACE include
)
target_link_libraries(mathbox-types
    INTERFACE Eigen3::Eigen cpp-utils datapack
)

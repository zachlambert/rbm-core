
find_package(Eigen3 3.4 REQUIRED NO_MODULE)

add_library(math_types INTERFACE)
target_include_directories(math_types
    INTERFACE include
)
target_link_libraries(math_types
    INTERFACE Eigen3::Eigen core_std core_serialize
)

if(BUILD_ADDITIONAL_TARGETS)
    add_executable(test_types_json test/json.cpp)
    target_link_libraries(test_json parrot GTest::gtest_main)
    gtest_discover_tests(test_json)
endif()

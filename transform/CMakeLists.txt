
add_library(owl_transform INTERFACE)
target_include_directories(owl_transform
    INTERFACE include
)
target_link_libraries(owl_transform
    INTERFACE owl_types parrot
)

if(BUILD_ADDITIONAL_TARGETS)
    add_executable(test_transform_angle test/angle.cpp)
    target_link_libraries(test_transform_angle owl_transform GTest::gtest_main)
    gtest_discover_tests(test_transform_angle)

    add_executable(test_transform_serialize test/serialize.cpp)
    target_link_libraries(test_transform_serialize owl_transform GTest::gtest_main)
    gtest_discover_tests(test_transform_serialize)

    add_executable(test_transform_transform2 test/transform2.cpp)
    target_link_libraries(test_transform_transform2 owl_transform GTest::gtest_main)
    gtest_discover_tests(test_transform_transform2)

    add_executable(test_transform_transform3 test/transform3.cpp)
    target_link_libraries(test_transform_transform3 owl_transform GTest::gtest_main)
    gtest_discover_tests(test_transform_transform3)

    add_executable(benchmark_transform_transform3 benchmark/transform3.cpp)
    target_link_libraries(benchmark_transform_transform3 owl_transform)
endif()

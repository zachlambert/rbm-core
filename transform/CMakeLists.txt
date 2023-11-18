
add_library(mathbox-transform INTERFACE)
target_include_directories(mathbox-transform
    INTERFACE include
)
target_link_libraries(mathbox-transform
    INTERFACE mathbox-types datapack
)

if(BUILD_ADDITIONAL_TARGETS)
    add_executable(test_transform_angle test/angle.cpp)
    target_link_libraries(test_transform_angle mathbox-transform GTest::gtest_main)
    gtest_discover_tests(test_transform_angle)

    add_executable(test_transform_datapack test/datapack.cpp)
    target_link_libraries(test_transform_datapack mathbox-transform GTest::gtest_main)
    gtest_discover_tests(test_transform_datapack)

    add_executable(test_transform_transform2 test/transform2.cpp)
    target_link_libraries(test_transform_transform2 mathbox-transform GTest::gtest_main)
    gtest_discover_tests(test_transform_transform2)

    add_executable(test_transform_transform3 test/transform3.cpp)
    target_link_libraries(test_transform_transform3 mathbox-transform GTest::gtest_main)
    gtest_discover_tests(test_transform_transform3)

    add_executable(benchmark_transform_transform3 benchmark/transform3.cpp)
    target_link_libraries(benchmark_transform_transform3 mathbox-transform)
endif()

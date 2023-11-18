
add_library(mbox-transform INTERFACE)
target_include_directories(mbox-transform INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)
target_link_libraries(mbox-transform
    INTERFACE mbox-types datapack
)

if(BUILD_ADDITIONAL_TARGETS)
    add_executable(test_transform_angle test/angle.cpp)
    target_link_libraries(test_transform_angle mbox-transform GTest::gtest_main)
    gtest_discover_tests(test_transform_angle)

    add_executable(test_transform_datapack test/datapack.cpp)
    target_link_libraries(test_transform_datapack mbox-transform GTest::gtest_main)
    gtest_discover_tests(test_transform_datapack)

    add_executable(test_transform_transform2 test/transform2.cpp)
    target_link_libraries(test_transform_transform2 mbox-transform GTest::gtest_main)
    gtest_discover_tests(test_transform_transform2)

    add_executable(test_transform_transform3 test/transform3.cpp)
    target_link_libraries(test_transform_transform3 mbox-transform GTest::gtest_main)
    gtest_discover_tests(test_transform_transform3)

    add_executable(benchmark_transform_transform3 benchmark/transform3.cpp)
    target_link_libraries(benchmark_transform_transform3 mbox-transform)
endif()


add_library(rbm-transform INTERFACE)
target_include_directories(rbm-transform INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)
target_link_libraries(rbm-transform
    INTERFACE rbm-types datapack
)

if(${BUILD_TESTS})
    add_executable(test_transform_angle test/angle.cpp)
    target_link_libraries(test_transform_angle rbm-transform GTest::gtest_main)
    gtest_discover_tests(test_transform_angle)

    add_executable(test_transform_datapack test/datapack.cpp)
    target_link_libraries(test_transform_datapack rbm-transform GTest::gtest_main)
    gtest_discover_tests(test_transform_datapack)

    add_executable(test_transform_transform2 test/transform2.cpp)
    target_link_libraries(test_transform_transform2 rbm-transform GTest::gtest_main)
    gtest_discover_tests(test_transform_transform2)

    add_executable(test_transform_transform3 test/transform3.cpp)
    target_link_libraries(test_transform_transform3 rbm-transform GTest::gtest_main)
    gtest_discover_tests(test_transform_transform3)
endif()

if(${BUILD_EXAMPLES})
    add_executable(example_transform_transform3 example/transform3.cpp)
    target_link_libraries(example_transform_transform3 rbm-transform)
endif()

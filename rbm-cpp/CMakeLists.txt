add_library(rbm-cpp INTERFACE)
target_include_directories(rbm-cpp INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)


if (${BUILD_TESTS})
    add_executable(test_cpp_optional_tuple test/optional_tuple.cpp)
    target_link_libraries(test_cpp_optional_tuple rbm-cpp GTest::gtest_main)
    gtest_discover_tests(test_cpp_optional_tuple)
endif()

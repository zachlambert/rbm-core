cmake_minimum_required(VERSION 3.14)
project(mathbox)

# Setup

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
include(FetchContent)

set(BUILD_ADDITIONAL_TARGETS OFF)
if(CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
    set(BUILD_ADDITIONAL_TARGETS ON)
endif()

if(BUILD_ADDITIONAL_TARGETS)
    FetchContent_Declare(googletest URL https://github.com/google/googletest/archive/03597a01ee50ed33e9dfd640b249b4be3799d395.zip)
    include(GoogleTest)
    set(INSTALL_GTEST OFF)
    set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
    FetchContent_MakeAvailable(googletest)
    enable_testing()
    # subdirectories will add their own tests
endif()

# Dependencies

find_package(Eigen3 3.4 REQUIRED NO_MODULE)
find_package(datapack REQUIRED)
find_package(cpp-utils REQUIRED)
if(${USE_SVIZ})
    find_package(sviz REQUIRED)
endif()

# Library

add_subdirectory(rbm-types)
add_subdirectory(rbm-transform)
add_subdirectory(rbm-diff)
add_subdirectory(rbm-geometry)
add_subdirectory(rbm-opt)
add_subdirectory(rbm-probability)
add_subdirectory(rbm-mechanics)

# Create install target
# https://stackoverflow.com/questions/47718485/install-and-export-interface-only-library-cmake

# include(CMakePackageConfigHelpers)
# write_basic_package_version_file(
#     "${PROJECT_BINARY_DIR}/rbmConfigVersion.cmake"
#     VERSION 0.1
#     COMPATIBILITY AnyNewerVersion
# )
#
# install(
#     TARGETS
#         rbm-types
#         rbm-transform
#         rbm-diff
#         rbm-geometry
#         rbm-opt
#         rbm-probability
#         rbm-mechanics
#     EXPORT rbmTargets
#     LIBRARY DESTINATION lib COMPONENT Runtime
#     ARCHIVE DESTINATION lib COMPONENT Development
#     RUNTIME DESTINATION bin COMPONENT Runtime
#     PUBLIC_HEADER DESTINATION include COMPONENT Development
#     BUNDLE DESTINATION bin COMPONENT Runtime
# )
#
# include(CMakePackageConfigHelpers)
# configure_package_config_file(
#     "${PROJECT_SOURCE_DIR}/cmake/rbmConfig.cmake.in"
#     "${PROJECT_BINARY_DIR}/rbmConfig.cmake"
#     INSTALL_DESTINATION lib/cmake/rbm
# )
#
# install(EXPORT rbmTargets DESTINATION lib/cmake/rbm)
# install(FILES "${PROJECT_BINARY_DIR}/rbmConfigVersion.cmake"
#               "${PROJECT_BINARY_DIR}/rbmConfig.cmake"
#         DESTINATION lib/cmake/rbm)
#
# install(DIRECTORY ${PROJECT_SOURCE_DIR}/diff/include/ DESTINATION include)
# install(DIRECTORY ${PROJECT_SOURCE_DIR}/geometry/include/ DESTINATION include)
# install(DIRECTORY ${PROJECT_SOURCE_DIR}/mechanics/include/ DESTINATION include)
# install(DIRECTORY ${PROJECT_SOURCE_DIR}/opt/include/ DESTINATION include)
# install(DIRECTORY ${PROJECT_SOURCE_DIR}/probability/include/ DESTINATION include)
# install(DIRECTORY ${PROJECT_SOURCE_DIR}/transform/include/ DESTINATION include)
# install(DIRECTORY ${PROJECT_SOURCE_DIR}/types/include/ DESTINATION include)

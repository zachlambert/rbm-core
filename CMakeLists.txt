cmake_minimum_required(VERSION 3.14)
project(owl)

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

# Library

add_subdirectory(types)
add_subdirectory(transform)
add_subdirectory(diff)
add_subdirectory(geometry)
add_subdirectory(opt)
add_subdirectory(probability)
add_subdirectory(mechanics)

add_library(mbox INTERFACE)
target_link_libraries(mbox INTERFACE
    mbox-types
    mbox-transform
    mbox-diff
    mbox-geometry
    mbox-opt
    mbox-probability
    mbox-mechanics
)

# Create install target
# https://stackoverflow.com/questions/47718485/install-and-export-interface-only-library-cmake

include(CMakePackageConfigHelpers)
write_basic_package_version_file(
    "${PROJECT_BINARY_DIR}/mboxConfigVersion.cmake"
    VERSION 0.1
    COMPATIBILITY AnyNewerVersion
)

install(
    TARGETS
        mbox-types
        mbox-transform
    EXPORT mboxTargets
    LIBRARY DESTINATION lib COMPONENT Runtime
    ARCHIVE DESTINATION lib COMPONENT Development
    RUNTIME DESTINATION bin COMPONENT Runtime
    PUBLIC_HEADER DESTINATION include COMPONENT Development
    BUNDLE DESTINATION bin COMPONENT Runtime
)

include(CMakePackageConfigHelpers)
configure_package_config_file(
    "${PROJECT_SOURCE_DIR}/cmake/mboxConfig.cmake.in"
    "${PROJECT_BINARY_DIR}/mboxConfig.cmake"
    INSTALL_DESTINATION lib/cmake/mbox
)

install(EXPORT mboxTargets DESTINATION lib/cmake/mbox)
install(FILES "${PROJECT_BINARY_DIR}/mboxConfigVersion.cmake"
              "${PROJECT_BINARY_DIR}/mboxConfig.cmake"
        DESTINATION lib/cmake/mbox)

install(DIRECTORY ${PROJECT_SOURCE_DIR}/types/include/ DESTINATION include)
install(DIRECTORY ${PROJECT_SOURCE_DIR}/transform/include/ DESTINATION include)

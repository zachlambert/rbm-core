
add_library(mbox-diff INTERFACE)
target_include_directories(mbox-diff
    INTERFACE include
)
target_link_libraries(mbox-diff
    INTERFACE mbox-types mbox-transform
)

add_executable(example_diff_dfunction example/dfunction.cpp)
target_link_libraries(example_diff_dfunction mbox-diff)

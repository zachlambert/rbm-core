
add_library(mathbox-diff INTERFACE)
target_include_directories(mathbox-diff
    INTERFACE include
)
target_link_libraries(mathbox-diff
    INTERFACE mathbox-types mathbox-transform
)

add_executable(example_diff_dfunction example/dfunction.cpp)
target_link_libraries(example_diff_dfunction mathbox-diff)

cmake_minimum_required(VERSION 3.0.0)

project(STG VERSION 0.1.0)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

add_compile_options(
    -Wall
    -Wold-style-cast
    -Warray-bounds
)

add_subdirectory(src)
file(INSTALL ./test_resources DESTINATION ${CMAKE_BINARY_DIR}/bin)
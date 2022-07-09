set(HEADERS TestsHeader.hpp)
set(SRC main.cpp 
    FEM_FinElemMesh_tests.cpp)

set(STG_tests STG_unit_tests)

add_executable(${STG_tests} ${SRC} ${HEADERS})

target_link_libraries(${STG_tests} PRIVATE ${STG_lib} ${CATCH2_FRAMEWORK})
target_include_directories(${STG_tests} PRIVATE ${CMAKE_SOURCE_DIR}/src/lib/include)


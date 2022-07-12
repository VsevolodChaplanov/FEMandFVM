set(HEADERS TestsHeader.hpp)
set(SRC main.cpp 
    FEM_FinElemMesh_tests.cpp)


add_executable(${PROJECT_NAME}_tests ${SRC} ${HEADERS})

target_link_libraries(${PROJECT_NAME}_tests PRIVATE ${PROJECT_NAME}_lib ${CATCH2_FRAMEWORK})
target_include_directories(${PROJECT_NAME}_tests PRIVATE ${CMAKE_SOURCE_DIR}/src/lib/include)


set(HEADERS 
    include/FEM_FiniteElement.hpp
    include/FEM_Mesh.hpp
    include/FEM_Parser.hpp
    include/FVM_FiniteVolumeElement.hpp
    include/FEM_MeshBuilder.hpp
)
set(SRC
    stg/FEM_FiniteElement.cpp
    stg/FEM_Mesh.cpp
    stg/FEM_Parser.cpp
    stg/FVM_FiniteVolumeElement.cpp
    stg/FEM_MeshBuilder.cpp
)

find_package(Boost REQUIRED)

add_library(${PROJECT_NAME}_lib SHARED ${HEADERS} ${SRC} )
target_include_directories(${PROJECT_NAME}_lib PRIVATE include)
target_include_directories(${PROJECT_NAME}_lib PRIVATE ${Boost_INCLUDE_DIR})
target_link_libraries(${PROJECT_NAME}_lib PRIVATE ${Boost_LIBARARIES})

add_subdirectory(tests)
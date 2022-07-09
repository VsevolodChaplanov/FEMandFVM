set(HEADERS 
    include/FEM_Mesh.hpp
    include/FEM_Parser.hpp
    include/FiniteElement.hpp
    include/FVM_FiniteVolumeElement.hpp
)
set(SRC
    stg/FEM_Mesh.cpp
    stg/FEM_Parser.cpp
    stg/FiniteElement.cpp
    stg/FVM_FiniteVolumeElement.cpp
)

find_package(Boost REQUIRED)

add_library(${PROJECT_NAME}_lib SHARED ${SRC} ${HEADERS})
target_include_directories(${PROJECT_NAME}_lib PRIVATE include)
target_include_directories(${PROJECT_NAME}_lib PRIVATE ${Boost_INCLUDE_DIR})
target_link_libraries(${PROJECT_NAME}_lib PRIVATE ${Boost_LIBARARIES})

add_subdirectory(tests)
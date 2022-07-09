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

set(STG_lib STG_ex)

find_package(Boost REQUIRED)

add_library(${STG_lib} SHARED ${SRC} ${HEADERS})
target_include_directories(${STG_lib} PRIVATE include)
target_include_directories(${STG_lib} PRIVATE ${Boost_INCLUDE_DIR})
target_link_libraries(${STG_lib} PRIVATE ${Boost_LIBARARIES})

add_subdirectory(tests)
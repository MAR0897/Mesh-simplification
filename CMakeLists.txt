cmake_minimum_required(VERSION 3.0)
project(my_project)

set(CMAKE_CXX_STANDARD 11)
set(PROGRAM_NAME main)
set(CPP_FILES main.cpp meshwrap.cpp constraints.cpp)

add_executable(${PROGRAM_NAME} ${CPP_FILES})

add_compile_options(-fsanitize=address)
add_link_options(-fsanitize=address)

find_package(Eigen3 REQUIRED)
target_include_directories(${PROGRAM_NAME} PRIVATE ${Eigen_INCLUDE_DIRS})
target_link_libraries(${PROGRAM_NAME} PRIVATE ${Eigen_LIBRARIES})

find_package(OpenMesh REQUIRED)
target_include_directories(${PROGRAM_NAME} PRIVATE ${OPENMESH_INCLUDE_DIRS})
target_link_libraries(${PROGRAM_NAME} PRIVATE ${OPENMESH_LIBRARIES})
link_directories(${OPENMESH_LIBRARY_DIRS})
target_link_libraries(${PROGRAM_NAME} PRIVATE OpenMeshCore OpenMeshTools)

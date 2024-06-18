cmake_minimum_required(VERSION 3.16)
project(my_project)

set(CMAKE_CXX_STANDARD 11)
set(PROGRAM_NAME main)
set(CPP_FILES main.cpp meshwrap.cpp constraints.cpp)
set(HEADER_FILES pch.h)
#Optimization (O0=nejmene az O3=nejvice, pak taky Ofast...)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Ofast") 


add_executable(${PROGRAM_NAME} ${CPP_FILES})


#Precompiled headers
#target_precompile_headers(${PROGRAM_NAME} PRIVATE ${HEADER_FILES})

# add_compile_options(-fsanitize=address)
# add_link_options(-fsanitize=address)

# pro profilovani:  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pg")

#Include Eigen library
find_package(Eigen3 REQUIRED)
target_include_directories(${PROGRAM_NAME} PRIVATE ${Eigen_INCLUDE_DIRS})
target_link_libraries(${PROGRAM_NAME} PRIVATE ${Eigen_LIBRARIES})

#Include OpenMeshLibrary
find_package(OpenMesh REQUIRED)
target_include_directories(${PROGRAM_NAME} PRIVATE ${OPENMESH_INCLUDE_DIRS})
target_link_libraries(${PROGRAM_NAME} PRIVATE ${OPENMESH_LIBRARIES})
link_directories(${OPENMESH_LIBRARY_DIRS})
target_link_libraries(${PROGRAM_NAME} PRIVATE OpenMeshCore OpenMeshTools)

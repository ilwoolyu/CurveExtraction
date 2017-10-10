cmake_minimum_required(VERSION 2.8)

project(CurveExtraction)
set(CMAKE_BUILD_TYPE Release)

# CMAKE
include_directories(${CMAKE_SOURCE_DIR} ${CMAKE_BINARY_DIR})

# GenerateCLP
find_package(SlicerExecutionModel REQUIRED)
include(${SlicerExecutionModel_USE_FILE})
include(${GenerateCLP_USE_FILE})
set(WRAPPER_MAIN wrapper/main.cpp)
generateCLP(WRAPPER_MAIN ${CMAKE_SOURCE_DIR}/wrapper/wrapper.xml)

# OpenMP
find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

# MESH library
find_package(MeshLib REQUIRED)
include_directories(${Mesh_INCLUDE_DIRS} ${Mesh_INCLUDE_DIRS}/Util)
link_directories(${Mesh_LIBRARY_DIRS})

# sources
include_directories(${CMAKE_SOURCE_DIR}/src)
set(SRC_FILES ${CMAKE_SOURCE_DIR}/src/GyralCurve.cpp ${CMAKE_SOURCE_DIR}/src/GyralPoint.cpp ${CMAKE_SOURCE_DIR}/src/SulcalCurve.cpp ${CMAKE_SOURCE_DIR}/src/SulcalPoint.cpp ${CMAKE_SOURCE_DIR}/src/Extraction.cpp)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ~/bin)

# executable file
add_executable(CurveExtraction ${SRC_FILES} ${WRAPPER_MAIN})
target_link_libraries(CurveExtraction ${Mesh_LIBRARIES})


cmake_minimum_required(VERSION 2.8)

project(CurveExtraction)

# CMAKE
include_directories(${CMAKE_SOURCE_DIR} ${CMAKE_BINARY_DIR})

# GenerateCLP
find_package(SlicerExecutionModel REQUIRED)
include(${SlicerExecutionModel_USE_FILE})
include(${GenerateCLP_USE_FILE})
generateCLP(WRAPPER_XML wrapper/wrapper.xml)
set(WRAPPER_MAIN wrapper/main.cpp)

# OpenMP
find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

# MESH library
set(MESH_LIB_DIR /home/local/VANDERBILT/lyui/Tools/MeshLib)
include_directories(${MESH_LIB_DIR}/include ${MESH_LIB_DIR}/include/Util)
#add_library(Mesh SHARED IMPORTED)
#set_target_properties(Mesh PROPERTIES IMPORTED_LOCATION ${MESH_LIB_DIR}/res/libMesh.so.1.0)
link_directories(${MESH_LIB_DIR}/lib)

# sources
include_directories(src)
#set(INCLUDE_FILES src/GyralCurve.h src/GyralPoint.h src/SulcalCurve.h src/SulcalPoint.h src/Extraction.h)
set(SRC_FILES src/GyralCurve.cpp src/GyralPoint.cpp src/SulcalCurve.cpp src/SulcalPoint.cpp src/Extraction.cpp)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ~/bin)

# executable file
add_executable(CurveExtraction ${SRC_FILES} ${WRAPPER_XML} ${WRAPPER_MAIN})
target_link_libraries(CurveExtraction Mesh)

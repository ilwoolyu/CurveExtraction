cmake_minimum_required(VERSION 3.5)

project(CurveExtraction)

# Set a default build type if none was specified
if (NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
        message(STATUS "Setting build type to 'Release' as none was specified.")
        set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build." FORCE)
        mark_as_advanced(CMAKE_BUILD_TYPE)
        # Set the possible values of build type for cmake-gui
        set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

# GenerateCLP
find_package(SlicerExecutionModel REQUIRED)
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

# sources
include_directories(${CMAKE_SOURCE_DIR}/src)
set(
	SRC_FILES
	${CMAKE_SOURCE_DIR}/src/GyralCurve.cpp
	${CMAKE_SOURCE_DIR}/src/GyralPoint.cpp
	${CMAKE_SOURCE_DIR}/src/SulcalCurve.cpp
	${CMAKE_SOURCE_DIR}/src/SulcalPoint.cpp
	${CMAKE_SOURCE_DIR}/src/Extraction.cpp
	)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ~/bin)

# executable file
add_executable(CurveExtraction ${SRC_FILES} ${WRAPPER_MAIN})

target_link_libraries(CurveExtraction Mesh)

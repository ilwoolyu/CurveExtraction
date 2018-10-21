include(${GenerateCLP_USE_FILE})

include_directories(${CMAKE_SOURCE_DIR}/src)

generateCLP(CurveExtraction_WRAPPER ${CMAKE_SOURCE_DIR}/wrapper/wrapper.xml)

set(
	CurveExtraction_SOURCE
	${CMAKE_SOURCE_DIR}/src/GyralCurve.cpp
	${CMAKE_SOURCE_DIR}/src/GyralPoint.cpp
	${CMAKE_SOURCE_DIR}/src/SulcalCurve.cpp
	${CMAKE_SOURCE_DIR}/src/SulcalPoint.cpp
	${CMAKE_SOURCE_DIR}/src/Extraction.cpp
	)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ~/bin)

add_executable(
	CurveExtraction
	${CurveExtraction_SOURCE}
	${CMAKE_SOURCE_DIR}/wrapper/main.cpp ${CurveExtraction_WRAPPER}
)
	
if (OPENMP_FOUND)
    target_compile_definitions(CurveExtraction PRIVATE _USE_OPENMP)
endif()

target_link_libraries(CurveExtraction Mesh)

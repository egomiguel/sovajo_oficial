cmake_minimum_required(VERSION 3.0.0)
project(UKA_Segmentation VERSION 1.0.1 LANGUAGES CXX)
include(GenerateExportHeader)
include(CMakePackageConfigHelpers)
#find_package(ITK COMPONENTS ITKTransform REQUIRED)
find_package(ITK REQUIRED)
include(${ITK_USE_FILE})
find_package(OpenCV REQUIRED)
include_directories(${OpenCV_INCLUDE_DIRS})
include_directories(eigen3 eigen3/unsupported)
include_directories(surface)
find_package(VTK REQUIRED)
include(${VTK_USE_FILE})
file(GLOB_RECURSE _HDRS "*.hpp")
file(GLOB_RECURSE _SRCS "*.cpp")
set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)
add_library(UKA_Segmentation SHARED ${_HDRS} ${_SRCS})
target_link_libraries(UKA_Segmentation PRIVATE ${ITK_LIBRARIES} ${OpenCV_LIBS} ${VTK_LIBRARIES})
generate_export_header(UKA_Segmentation)
target_include_directories(${PROJECT_NAME} PUBLIC 
							"$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>" 
							"$<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>" )

							install(TARGETS UKA_Segmentation
	EXPORT UKA_SegmentationTargets
	RUNTIME DESTINATION bin
	LIBRARY DESTINATION lib
	ARCHIVE DESTINATION lib
	INCLUDES DESTINATION include)
install(FILES ${_HDRS} DESTINATION include)
install(FILES ${PROJECT_BINARY_DIR}/uka_segmentation_export.h DESTINATION include)
install(EXPORT UKA_SegmentationTargets
        FILE UKA_SegmentationTargets.cmake
        DESTINATION lib/cmake/UKA_Segmentation)
configure_package_config_file(${CMAKE_CURRENT_SOURCE_DIR}/UKA_SegmentationConfig.cmake.in
	"${CMAKE_CURRENT_BINARY_DIR}/UKA_SegmentationConfig.cmake"
	INSTALL_DESTINATION lib/cmake/UKA_Segmentation)
install(FILES
	"${CMAKE_CURRENT_BINARY_DIR}/UKA_SegmentationConfig.cmake"
	DESTINATION lib/cmake/UKA_Segmentation)
set(CPACK_GENERATOR NSIS)
set(CPACK_PACKAGE_FILE_NAME "${CMAKE_PROJECT_NAME}")
set(CPACK_PACKAGE_VENDOR "CMake.org")
set(CPACK_PACKAGE_VERSION "1.0.1")
set(CPACK_PACKAGE_VERSION_MAJOR "1")
set(CPACK_PACKAGE_VERSION_MINOR "0")
set(CPACK_PACKAGE_VERSION_PATCH "1")
set(CPACK_PACKAGE_INSTALL_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}")
set(CPACK_GENERATOR "ZIP")
SET(CPACK_NSIS_MODIFY_PATH ON)
INCLUDE(CPack)
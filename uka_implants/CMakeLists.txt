cmake_minimum_required(VERSION 3.0.0)
project(UKA_Implants VERSION 1.0.1 LANGUAGES CXX)
include(GenerateExportHeader)
include(CMakePackageConfigHelpers)
#find_package(ITK COMPONENTS ITKTransform REQUIRED)
find_package(ITK REQUIRED)
include(${ITK_USE_FILE})
find_package(OpenCV REQUIRED)
include_directories(${OpenCV_INCLUDE_DIRS})
include_directories(eigen3 eigen3/unsupported)
find_package(VTK REQUIRED)
include(${VTK_USE_FILE})
set(Boost_COMPILER "-vc142")
#find_package(PCL 1.10 REQUIRED COMPONENTS io filters)
#include_directories(${PCL_INCLUDE_DIRS})
#link_directories(${PCL_LIBRARY_DIRS})
#add_definitions(${PCL_DEFINITIONS})
file(GLOB_RECURSE _HDRS "*.hpp")
file(GLOB_RECURSE _SRCS "*.cpp")
set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)
add_library(UKA_Implants SHARED ${_HDRS} ${_SRCS})
#target_link_libraries(UKA_Implants PRIVATE ${ITK_LIBRARIES} ${OpenCV_LIBS} ${VTK_LIBRARIES} ${PCL_LIBRARIES})
target_link_libraries(UKA_Implants PRIVATE ${ITK_LIBRARIES} ${OpenCV_LIBS} ${VTK_LIBRARIES})

generate_export_header(UKA_Implants)
target_include_directories(${PROJECT_NAME} PUBLIC 
							"$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>" 
							"$<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>" )
install(TARGETS UKA_Implants
	EXPORT UKA_ImplantsTargets
	RUNTIME DESTINATION bin
	LIBRARY DESTINATION lib
	ARCHIVE DESTINATION lib
	INCLUDES DESTINATION include)
install(FILES ${_HDRS} DESTINATION include)
install(FILES ${PROJECT_BINARY_DIR}/uka_implants_export.h DESTINATION include)
install(EXPORT UKA_ImplantsTargets
        FILE UKA_ImplantsTargets.cmake
        DESTINATION lib/cmake/UKA_Implants)
configure_package_config_file(${CMAKE_CURRENT_SOURCE_DIR}/UKA_ImplantsConfig.cmake.in
	"${CMAKE_CURRENT_BINARY_DIR}/UKA_ImplantsConfig.cmake"
	INSTALL_DESTINATION lib/cmake/UKA_Implants)
install(FILES
	"${CMAKE_CURRENT_BINARY_DIR}/UKA_ImplantsConfig.cmake"
	DESTINATION lib/cmake/UKA_Implants)
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

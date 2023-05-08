cmake_minimum_required(VERSION 3.0.0)
project(THA_Registration VERSION 1.0.1 LANGUAGES CXX)
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
#find_package(PCL 1.10 REQUIRED COMPONENTS commom io visualization filters)
#find_package(PCL 1.10 REQUIRED COMPONENTS commom io filters)
#include_directories(${PCL_INCLUDE_DIRS})
#link_directories(${PCL_LIBRARY_DIRS})
#add_definitions(${PCL_DEFINITIONS})
file(GLOB_RECURSE _HDRS "*.hpp")
file(GLOB_RECURSE _SRCS "*.cpp")
set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)
add_library(THA_Registration SHARED ${_HDRS} ${_SRCS})
#target_link_libraries(THA_Registration PRIVATE ${ITK_LIBRARIES} ${OpenCV_LIBS} ${VTK_LIBRARIES} ${PCL_LIBRARIES})
target_link_libraries(THA_Registration PRIVATE ${ITK_LIBRARIES} ${OpenCV_LIBS} ${VTK_LIBRARIES})

generate_export_header(THA_Registration)
target_include_directories(${PROJECT_NAME} PUBLIC 
							"$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>" 
							"$<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>" )
install(TARGETS THA_Registration
	EXPORT THA_RegistrationTargets
	RUNTIME DESTINATION bin
	LIBRARY DESTINATION lib
	ARCHIVE DESTINATION lib
	INCLUDES DESTINATION include)
install(FILES ${_HDRS} DESTINATION include)
install(FILES ${PROJECT_BINARY_DIR}/tha_registration_export.h DESTINATION include)
install(EXPORT THA_RegistrationTargets
        FILE THA_RegistrationTargets.cmake
        DESTINATION lib/cmake/THA_Registration)
configure_package_config_file(${CMAKE_CURRENT_SOURCE_DIR}/THA_RegistrationConfig.cmake.in
	"${CMAKE_CURRENT_BINARY_DIR}/THA_RegistrationConfig.cmake"
	INSTALL_DESTINATION lib/cmake/THA_Registration)
install(FILES
	"${CMAKE_CURRENT_BINARY_DIR}/THA_RegistrationConfig.cmake"
	DESTINATION lib/cmake/THA_Registration)
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
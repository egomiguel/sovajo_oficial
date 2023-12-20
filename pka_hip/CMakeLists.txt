cmake_minimum_required(VERSION 3.0.0)
project(PKA_Hip VERSION 1.0.1 LANGUAGES CXX)
include(GenerateExportHeader)
include(CMakePackageConfigHelpers)
find_package(OpenCV REQUIRED)
include_directories(${OpenCV_INCLUDE_DIRS})
file(GLOB_RECURSE _HDRS "*.hpp")
file(GLOB_RECURSE _SRCS "*.cpp")
set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)
add_library(PKA_Hip SHARED ${_HDRS} ${_SRCS})
target_link_libraries(PKA_Hip ${OpenCV_LIBS})
generate_export_header(PKA_Hip)
target_include_directories(${PROJECT_NAME} PUBLIC 
							"$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>" 
							"$<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>" )
install(TARGETS PKA_Hip
	EXPORT PKA_HipTargets
	RUNTIME DESTINATION bin
	LIBRARY DESTINATION lib
	ARCHIVE DESTINATION lib
	INCLUDES DESTINATION include)
install(FILES ${_HDRS} DESTINATION include)
install(FILES ${PROJECT_BINARY_DIR}/pka_hip_export.h DESTINATION include)
install(EXPORT PKA_HipTargets
        FILE PKA_HipTargets.cmake
        DESTINATION lib/cmake/PKA_Hip)
configure_package_config_file(${CMAKE_CURRENT_SOURCE_DIR}/PKA_HipConfig.cmake.in
	"${CMAKE_CURRENT_BINARY_DIR}/PKA_HipConfig.cmake"
	INSTALL_DESTINATION lib/cmake/PKA_Hip)
install(FILES
	"${CMAKE_CURRENT_BINARY_DIR}/PKA_HipConfig.cmake"
	DESTINATION lib/cmake/PKA_Hip)
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

cmake_minimum_required(VERSION 3.0.0)
project(UKA_Hip VERSION 1.0.1 LANGUAGES CXX)
include(GenerateExportHeader)
include(CMakePackageConfigHelpers)
find_package(OpenCV REQUIRED)
include_directories(${OpenCV_INCLUDE_DIRS})
file(GLOB_RECURSE _HDRS "*.hpp")
file(GLOB_RECURSE _SRCS "*.cpp")
set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)
add_library(UKA_Hip SHARED ${_HDRS} ${_SRCS})
target_link_libraries(UKA_Hip ${OpenCV_LIBS})
generate_export_header(UKA_Hip)
target_include_directories(${PROJECT_NAME} PUBLIC 
							"$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>" 
							"$<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>" )
install(TARGETS UKA_Hip
	EXPORT UKA_HipTargets
	RUNTIME DESTINATION bin
	LIBRARY DESTINATION lib
	ARCHIVE DESTINATION lib
	INCLUDES DESTINATION include)
install(FILES ${_HDRS} DESTINATION include)
install(FILES ${PROJECT_BINARY_DIR}/uka_hip_export.h DESTINATION include)
install(EXPORT UKA_HipTargets
        FILE UKA_HipTargets.cmake
        DESTINATION lib/cmake/UKA_Hip)
configure_package_config_file(${CMAKE_CURRENT_SOURCE_DIR}/UKA_HipConfig.cmake.in
	"${CMAKE_CURRENT_BINARY_DIR}/UKA_HipConfig.cmake"
	INSTALL_DESTINATION lib/cmake/UKA_Hip)
install(FILES
	"${CMAKE_CURRENT_BINARY_DIR}/UKA_HipConfig.cmake"
	DESTINATION lib/cmake/UKA_Hip)
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

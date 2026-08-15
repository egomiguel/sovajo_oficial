cmake_minimum_required(VERSION 3.0.0)
project(TKA_Hip_Nav VERSION 1.0.1 LANGUAGES CXX)
include(GenerateExportHeader)
include(CMakePackageConfigHelpers)
find_package(OpenCV REQUIRED)
include_directories(${OpenCV_INCLUDE_DIRS})
file(GLOB_RECURSE _HDRS "*.hpp")
file(GLOB_RECURSE _SRCS "*.cpp")
set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)
add_library(TKA_Hip_Nav SHARED ${_HDRS} ${_SRCS})
target_link_libraries(TKA_Hip_Nav ${OpenCV_LIBS})
generate_export_header(TKA_Hip_Nav)
target_include_directories(${PROJECT_NAME} PUBLIC 
							"$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>" 
							"$<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>" )
install(TARGETS TKA_Hip_Nav
	EXPORT TKA_Hip_NavTargets
	RUNTIME DESTINATION bin
	LIBRARY DESTINATION lib
	ARCHIVE DESTINATION lib
	INCLUDES DESTINATION include)
install(FILES ${_HDRS} DESTINATION include)
install(FILES ${PROJECT_BINARY_DIR}/tka_hip_nav_export.h DESTINATION include)
install(EXPORT TKA_Hip_NavTargets
        FILE TKA_Hip_NavTargets.cmake
        DESTINATION lib/cmake/TKA_Hip_Nav)
configure_package_config_file(${CMAKE_CURRENT_SOURCE_DIR}/TKA_Hip_NavConfig.cmake.in
	"${CMAKE_CURRENT_BINARY_DIR}/TKA_Hip_NavConfig.cmake"
	INSTALL_DESTINATION lib/cmake/TKA_Hip_Nav)
install(FILES
	"${CMAKE_CURRENT_BINARY_DIR}/TKA_Hip_NavConfig.cmake"
	DESTINATION lib/cmake/TKA_Hip_Nav)
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
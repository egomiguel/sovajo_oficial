cmake_minimum_required(VERSION 3.0.0)
project(LoadFiles VERSION 1.0.1 LANGUAGES CXX)
include(GenerateExportHeader)
include(CMakePackageConfigHelpers)
find_package(Qt5 REQUIRED Core)
file(GLOB_RECURSE _HDRS "*.hpp")
file(GLOB_RECURSE _SRCS "*.cpp")
set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)
add_library(LoadFiles SHARED ${_HDRS} ${_SRCS})
target_link_libraries(LoadFiles Qt5::Core)
generate_export_header(LoadFiles)
target_include_directories(${PROJECT_NAME} PUBLIC 
							"$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>" 
							"$<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>" )
install(TARGETS LoadFiles
	EXPORT LoadFilesTargets
	RUNTIME DESTINATION bin
	LIBRARY DESTINATION lib
	ARCHIVE DESTINATION lib
	INCLUDES DESTINATION include)
install(FILES ${_HDRS} DESTINATION include)
install(FILES ${PROJECT_BINARY_DIR}/LoadFiles_export.h DESTINATION include)
install(EXPORT LoadFilesTargets
        FILE LoadFilesTargets.cmake
        DESTINATION lib/cmake/LoadFiles)
configure_package_config_file(${CMAKE_CURRENT_SOURCE_DIR}/LoadFilesConfig.cmake.in
	"${CMAKE_CURRENT_BINARY_DIR}/LoadFilesConfig.cmake"
	INSTALL_DESTINATION lib/cmake/LoadFiles)
install(FILES
	"${CMAKE_CURRENT_BINARY_DIR}/LoadFilesConfig.cmake"
	DESTINATION lib/cmake/LoadFiles)
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
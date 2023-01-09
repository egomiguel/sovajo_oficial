cmake_minimum_required(VERSION 3.0.0)
project(CheckImageRod VERSION 1.0.1 LANGUAGES CXX)
include(GenerateExportHeader)
include(CMakePackageConfigHelpers)
find_package(ITK REQUIRED)
include(${ITK_USE_FILE})
file(GLOB_RECURSE _HDRS "*.hpp")
file(GLOB_RECURSE _SRCS "*.cpp")
set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)
add_library(CheckImageRod SHARED ${_HDRS} ${_SRCS})
target_link_libraries(CheckImageRod PRIVATE ${ITK_LIBRARIES})
generate_export_header(CheckImageRod)
target_include_directories(${PROJECT_NAME} PUBLIC 
							"$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>" 
							"$<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>" )
install(TARGETS CheckImageRod
	EXPORT CheckImageRodTargets
	RUNTIME DESTINATION bin
	LIBRARY DESTINATION lib
	ARCHIVE DESTINATION lib
	INCLUDES DESTINATION include)
install(FILES ${_HDRS} DESTINATION include)
install(FILES ${PROJECT_BINARY_DIR}/CheckImageRod_export.h DESTINATION include)
install(EXPORT CheckImageRodTargets
        FILE CheckImageRodTargets.cmake
        DESTINATION lib/cmake/CheckImageRod)
configure_package_config_file(${CMAKE_CURRENT_SOURCE_DIR}/CheckImageRodConfig.cmake.in
	"${CMAKE_CURRENT_BINARY_DIR}/CheckImageRodConfig.cmake"
	INSTALL_DESTINATION lib/cmake/CheckImageRod)
install(FILES
	"${CMAKE_CURRENT_BINARY_DIR}/CheckImageRodConfig.cmake"
	DESTINATION lib/cmake/CheckImageRod)
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
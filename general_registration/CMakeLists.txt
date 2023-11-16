cmake_minimum_required(VERSION 3.0.0)
project(General_Registration VERSION 1.0.1 LANGUAGES CXX)
include(GenerateExportHeader)
include(CMakePackageConfigHelpers)
find_package(ITK REQUIRED)
include(${ITK_USE_FILE})
find_package(OpenCV REQUIRED)
include_directories(${OpenCV_INCLUDE_DIRS})
include_directories(eigen3 eigen3/unsupported)
find_package(VTK REQUIRED)
include(${VTK_USE_FILE})
set(Boost_COMPILER "-vc142")
file(GLOB_RECURSE _HDRS "*.hpp")
file(GLOB_RECURSE _SRCS "*.cpp")
set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)
add_library(General_Registration SHARED ${_HDRS} ${_SRCS})
target_link_libraries(General_Registration PRIVATE ${OpenCV_LIBS} ${ITK_LIBRARIES} ${VTK_LIBRARIES})

generate_export_header(General_Registration)
target_include_directories(${PROJECT_NAME} PUBLIC 
							"$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>" 
							"$<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>" )
install(TARGETS General_Registration
	EXPORT General_RegistrationTargets
	RUNTIME DESTINATION bin
	LIBRARY DESTINATION lib
	ARCHIVE DESTINATION lib
	INCLUDES DESTINATION include)
install(FILES ${_HDRS} DESTINATION include)
install(FILES ${PROJECT_BINARY_DIR}/general_registration_export.h DESTINATION include)
install(EXPORT General_RegistrationTargets
        FILE General_RegistrationTargets.cmake
        DESTINATION lib/cmake/General_Registration)
configure_package_config_file(${CMAKE_CURRENT_SOURCE_DIR}/General_RegistrationConfig.cmake.in
	"${CMAKE_CURRENT_BINARY_DIR}/General_RegistrationConfig.cmake"
	INSTALL_DESTINATION lib/cmake/General_Registration)
install(FILES
	"${CMAKE_CURRENT_BINARY_DIR}/General_RegistrationConfig.cmake"
	DESTINATION lib/cmake/General_Registration)
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

cmake_minimum_required(VERSION 3.7.2)

project(WatershedNative)
set(PROJECT_VERSION_MAJOR 0)
set(PROJECT_VERSION_MINOR 1)
set(PROJECT_VERSION_PATCH 0)

# Manual build type selection (for debugging purposes)
#set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_BUILD_TYPE Release)

find_package(JNI REQUIRED)

string(TOUPPER ${CMAKE_BUILD_TYPE} CMAKE_BUILD_TYPE_UPPER)

# Include directories
include_directories(${JNI_INCLUDE_DIRS})

# Sources
file(GLOB LIB_SRC "*.cpp")

set(LIB_NAME ${PROJECT_NAME})
add_library(${LIB_NAME} SHARED ${LIB_SRC})

install(TARGETS ${LIB_NAME} LIBRARY DESTINATION .)

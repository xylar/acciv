#The CMake Minimum version that is required.
CMAKE_MINIMUM_REQUIRED(VERSION 2.8)



#Detect Build Type if Debug appears in CMAKE_BINARY_DIR path
STRING(REGEX MATCH ".*(Debug|debug).*" DEBUG_MATCH ${CMAKE_BINARY_DIR} )
if(DEBUG_MATCH)
    set(CMAKE_BUILD_TYPE Debug CACHE STRING "The build type for makefile based generators")
    SET(PROJECTNAMEPREFIX Debug)
ELSE()
    set(CMAKE_BUILD_TYPE Release CACHE STRING "The build type for makefile based generators")
    SET(PROJECTNAMEPREFIX Release)
ENDIF()



set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake/;${CMAKE_MODULE_PATH}")
message(STATUS "Module path is now set to: " ${CMAKE_MODULE_PATH} )



find_package(FFTW REQUIRED)
find_package(HDF5 COMPONENTS CXX HL REQUIRED)


SET(CMAKE_DEBUG_POSTFIX "_d")


ADD_SUBDIRECTORY(acciv)

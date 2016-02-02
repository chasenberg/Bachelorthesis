
cmake_minimum_required(VERSION 2.8)

project(bachelor)
set(PROJECT_VERSION 0.0.0.0)
set(PROJECT_BRIEF "\"C++ Template\"")

# Enable/disable testing
option(TESTS_ENABLED "Enable automatic tests" ON)
# Enable/disable test coverage
option(COVERAGE_ENABLED "Enable test coverage" OFF)


set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/")
set(LIBRARY_OUTPUT_PATH "${CMAKE_BINARY_DIR}/lib")
set(EXECUTABLE_OUTPUT_PATH "${CMAKE_BINARY_DIR}/bin")

include_directories(${CMAKE_SOURCE_DIR}/src)

if(COVERAGE_ENABLED)
    set(CMAKE_CXX_FLAGS_DEBUG "-Wall -fPIC -O0 -ggdb -fprofile-arcs -ftest-coverage")
else()
    set(CMAKE_CXX_FLAGS_DEBUG "-Wall -O0 -ggdb")
endif()

set(CMAKE_CXX_FLAGS_RELEASE "-Wall -O2")

include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
if(COMPILER_SUPPORTS_CXX11)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
elseif(COMPILER_SUPPORTS_CXX0X)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
else()
    message(FATAL_ERROR "Compiler does not support C++11")
endif()

#find_package(Doxygen)
#if(DOXYGEN_FOUND)
#   	configure_file(${CMAKE_CURRENT_SOURCE_DIR}/docs/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)
#  	add_custom_target(docs
#     	${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
#    	WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
#   	COMMENT "Generating API documentation with Doxygen" VERBATIM
#)
#endif(DOXYGEN_FOUND)

set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if (NOT CMAKE_BUILD_TYPE)
    message(STATUS "No build type selected, defaulting to Debug")
    set(CMAKE_BUILD_TYPE "Debug")
endif()

# Boost, a collection of useful C++ utilities and idioms
find_package(Boost 1.48 COMPONENTS
    program_options
    filesystem
    thread
    system
    random
    REQUIRED)
include_directories(${Boost_INCLUDE_DIR})
# As the older 'ld' linker version shows a different behaviour,
# boost is not included on PhiDo trough this rather weak system info comparison.
if(${CMAKE_SYSTEM} STREQUAL "Linux-2.6.32-504.8.1.el6.x86_64")
    message(STATUS "We are on PhiDo: do not include Boost_LIBRARIES")
else()
    message(STATUS "We are not on PhiDo: hurry, include Boost_LIBRARIES")
    set(LIBRARIES ${LIBRARIES} ${Boost_LIBRARIES})
endif()

# ROOT, a library for HEP data analysis (contains RooFit)
find_package(ROOT REQUIRED)
include_directories(${ROOT_INCLUDE_DIR})
set(LIBRARIES ${LIBRARIES} ${ROOT_LIBRARIES})
set(LIBRARIES ${LIBRARIES} ${ROOFIT_LIBRARIES})
 
# Eigen3, a template library for high performance linear algebra
#find_package(Eigen3 REQUIRED)
#include_directories(EIGEN3_INCLUDE_DIR)

# glog, a powerful logging framework
#find_package(Glog REQUIRED)
#include_directories(${GLOG_INCLUDE_DIRS})
#set(LIBRARIES ${LIBRARIES} ${GLOG_LIBRARIES})

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "$ENV{DOOFITSYS}/cmake/Modules/")
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "$ENV{DOOCORESYS}/cmake/Modules/")
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "$ENV{DOOSELECTIONSYS}/cmake/Modules/")

find_package(DooCore)
find_package(DooFit)
find_package(DooSelection)

include_directories(${DooCore_INCLUDE_DIR})
link_directories(${DooCore_LIBRARY_DIR})

include_directories(${DooFit_INCLUDE_DIR})
link_directories(${DooFit_LIBRARY_DIR})

include_directories(${DooSelection_INCLUDE_DIR})
link_directories(${DooSelection_LIBRARY_DIR})

set(LIBRARIES ${LIBRARIES} ${DooCore_LIBRARIES} ${DooFit_LIBRARIES} ${DooSelection_LIBRARIES})

add_subdirectory(src)
add_subdirectory(main)
add_subdirectory(programmierkurs)
# vim:ft=cmake


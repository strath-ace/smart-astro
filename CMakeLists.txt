cmake_minimum_required(VERSION 3.0 FATAL_ERROR)

# Change this section to correspond with your project info.
project (smart-astro
        LANGUAGES C CXX
        VERSION 0.1.0)

set(CMAKE_INCLUDE_SYSTEM_FLAG_C "-isystem")
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_VERBOSE_MAKEFILE OFF)

set(PROJECT_SUMMARY                            "Strathclyde Mechanical and Aerospace Research Toolbox for Astrodynamics (SMART-ASTRO)")
set(${PROJECT_NAME}_VENDOR_NAME                "spaceART")
set(${PROJECT_NAME}_VENDOR_CONTACT             "spaceART@strathclyde.ac.uk")

set(PROJECT_PATH                               "${CMAKE_CURRENT_SOURCE_DIR}")
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH}     "${PROJECT_PATH}/cmake/Modules")
set(INCLUDE_PATH                               "${PROJECT_PATH}/src")
set(SRC_PATH                                   "${PROJECT_PATH}/src")
set(TEST_SRC_PATH                              "${PROJECT_PATH}/test")
if(NOT EXTERNAL_PATH)
    set(EXTERNAL_PATH                            "${PROJECT_PATH}/external")
endif()
if(NOT DOCS_PATH)
    set(DOCS_PATH                                "${PROJECT_PATH}/docs")
endif()
set(LIB_NAME                                   "smart-astro")
set(LIB_NAME_STATIC                            "smart-astro-static")
set(LIB_PATH                                   "${PROJECT_BINARY_DIR}/lib")
set(BIN_PATH                                   "${PROJECT_BINARY_DIR}/bin")
set(MAIN_NAME                                  "main")
set(TEST_PATH                                  "${PROJECT_BINARY_DIR}/test")
set(TEST_NAME                                  "test")
set(EXAMPLES_PATH                              "${PROJECT_BINARY_DIR}/examples")
set(EXAMPLES_NAME                              "examples")
set(IRILIB_PATH                                "${PROJECT_PATH}/irilib")
set(IRILIB_NAME                                "IRILIB")
if (NOT MOID_DIR)
    set(MOID_DIR                                "${EXTERNAL_PATH}/moid_distrib")
endif()

option(BUILD_STATIC                            "Build static library"           OFF)
option(BUILD_MAIN                              "Build main function"            ON)
option(BUILD_DOCS                              "Build docs"                     OFF)
option(BUILD_EXAMPLES                          "Build examples"                 ON)
option(ENABLE_SMARTUQ            	           "UQ methods"                     OFF)
option(ENABLE_CSPICE            	           "JPL CSPICE"                     ON)
option(ENABLE_WORHP            		           "eNLP solver"                    OFF)
option(ENABLE_GRAPHICS                         "X11 libraries (only UNIX systems)" OFF)
option(ENABLE_MONGODB                          "MongoDB database"               OFF)

#if ("$ENV{CI}" STREQUAL "true")
#    set(before_install_output_file before_install_output.log)
#    set(before_install_error_file before_install_error.log)
#
#    message(STATUS "Executing before_install.sh script")
#    execute_process(COMMAND "${CMAKE_CURRENT_SOURCE_DIR}/before_install.sh"
#            OUTPUT_FILE "${before_install_output_file}"
#            ERROR_FILE "${before_install_error_file}"
#            RESULT_VARIABLE before_install_code)
#    if ( NOT ${before_install_code} EQUAL 0 )
#        message(FATAL_ERROR "before_install.sh script failed with error code ${before_install_code}. See ${before_install_output_file} and ${before_install_error_file} logs for more details.")
#    endif()
#
#    mark_as_advanced(before_install_code before_install_output_file before_install_error_file)
#endif()

set(CMAKE_REQUIRED_LIBRARIES gfortran)
include_directories("${IRILIB_PATH}")
add_subdirectory("${IRILIB_PATH}")
file(COPY irilib/data_source DESTINATION ${PROJECT_BINARY_DIR})

include("cmake/Utils.cmake")
build_template ("${INCLUDE_PATH}/smartastro.h.in" "${INCLUDE_PATH}/smartastro.h")

# Initialise (empty) list of libraries to link against.
SET(MANDATORY_LIBRARIES "")

# mandatory dependancy to smartmath
include(ExternalProject)

include(cmake/CMakeLists.txt)
#include(cmake/External/smart-math.cmake)
#include_directories ("${SMART_MATH_INCLUDE_DIR}")
#
#if (ENABLE_SMARTUQ)
#    include(cmake/External/smart-uq.cmake)
#    include_directories ("${SMART_UQ_INCLUDE_DIR}")
#endif()
#
#
#find_package(OpenMP)
#if(OPENMP_FOUND)
#    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
#    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
#endif()
#
#if((CMAKE_COMPILER_IS_GNUCC) AND (CMAKE_BUILD_TYPE STREQUAL "Debug"))
#    # smart-uq requires __gcov_init and __gcov_merge_add symbols
#    set(CMAKE_C_FLAGS  "${CMAKE_C_FLAGS} --coverage")
#    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --coverage")
#    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -lgcov --coverage")
#endif()
#
## Graphics libraries
#if(ENABLE_GRAPHICS)
#    set(MANDATORY_LIBRARIES ${MANDATORY_LIBRARIES} "X11" "Xext")
#endif(ENABLE_GRAPHICS)
#
#if(ENABLE_WORHP)
#    find_library(WORHP_LIBRARY NAMES worhp)
#    find_path(WORHP_INCLUDE_DIR NAMES worhp/worhp.h)
#    if(NOT WORHP_LIBRARY OR NOT WORHP_INCLUDE_DIR)
#        message(FATAL_ERROR "WORHP support was requested, but the worhp library or headers could not be located.")
#    endif(NOT WORHP_LIBRARY OR NOT WORHP_INCLUDE_DIR)
#    message(STATUS "worhp library: ${WORHP_LIBRARY}")
#    message(STATUS "worhp headers: ${WORHP_INCLUDE_DIR}")
#
#    set(MANDATORY_LIBRARIES ${MANDATORY_LIBRARIES} ${WORHP_LIBRARY})
#    include_directories("${WORHP_INCLUDE_DIR}")
#endif()

file(GLOB_RECURSE SRC
        "${SRC_PATH}/*.cpp")

set(MAIN_SRC
        "${SRC_PATH}/main.cpp")

file(GLOB_RECURSE TEST_SRC
        "${TEST_SRC_PATH}/*.cpp")

file(GLOB_RECURSE LibFiles
        "${INCLUDE_PATH}/*.h"
        )


add_custom_target(headers SOURCES ${LibFiles})

# Code for compiling and executing MOID computation (Minimum Orbit Intersection Distance)
# The code is written in Fortran
set(MOID_INPUT ${MOID_DIR}/orbit_data.inp)
set(MOID_EXEC ${MOID_DIR}/CP_comp.exe)
set(MOID_OUTPUT ${MOID_DIR}/crtpts.out)

# Set CMake build-type. If it not supplied by the user, the default built type is "Release".
if(((NOT CMAKE_BUILD_TYPE))
        OR (CMAKE_BUILD_TYPE STREQUAL "Release"))
    message(STATUS "WARNING: building release version!")
    set(CMAKE_BUILD_TYPE Release)
elseif(CMAKE_BUILD_TYPE STREQUAL "Debug")
    message(STATUS "WARNING: building debug version!")
endif()

## Set platform-specific compiler flags.
#if(WIN32)
#    if(MSVC)
#        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /EHsc /Ox /W3 /FC -D_SCL_SECURE_NO_WARNINGS -std=c++11")
#        # We need to force MSVC to use static runtime:
#        # (from http://www.cmake.org/Wiki/CMake_FAQ#Dynamic_Replace).
#        foreach(flag_var
#                CMAKE_CXX_FLAGS CMAKE_CXX_FLAGS_DEBUG CMAKE_CXX_FLAGS_RELEASE
#                CMAKE_CXX_FLAGS_MINSIZEREL CMAKE_CXX_FLAGS_RELWITHDEBINFO)
#            # Find all dynamic runtime (MD) references and replace with static (MT)
#            if(${flag_var} MATCHES "/MD")
#                string(REGEX REPLACE "/MD" "/MT" ${flag_var} "${${flag_var}}")
#            endif()
#        endforeach()
#    else()
#        set(CMAKE_LEGACY_CYGWIN_WIN32 0)
#        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Woverloaded-virtual -Wnon-virtual-dtor")
#    endif()
#else()
#    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Woverloaded-virtual -Wnon-virtual-dtor")
#endif()
#
#if(CMAKE_COMPILER_IS_GNUCXX)
#    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0 -g3")
#    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")
#
#    #Flags for compiling mongocxx driver
#    if(ENABLE_MONGODB)
#        find_package(libmongocxx REQUIRED)
#        find_package(libbsoncxx REQUIRED)
#        include_directories(${LIBMONGOCXX_INCLUDE_DIR})
#        include_directories(${LIBBSONCXX_INCLUDE_DIR})
#        include_directories("${MONGOCXX_INSTALL_PREFIX}/include/mongocxx/v_noabi")
#        include_directories("${MONGOCXX_INSTALL_PREFIX}/include/bsoncxx/v_noabi")
#        include_directories("${MONGOCXX_INSTALL_PREFIX}/include/libmongoc-1.0")
#        include_directories("${MONGOCXX_INSTALL_PREFIX}/include/libbson-1.0")
#        include_directories("${MONGOCXX_INSTALL_PREFIX}/lib")
#        set(MANDATORY_LIBRARIES ${MANDATORY_LIBRARIES} ${LIBMONGOCXX_LIBRARIES} ${LIBBSONCXX_LIBRARIES})
#    endif()
#endif()

#build_template("${INCLUDE_PATH}/config.h.in" "${INCLUDE_PATH}/config.h")

include_directories(src)
include_directories(external)

if(BUILD_STATIC)
    add_library(${LIB_NAME_STATIC} STATIC ${SRC})
    target_link_libraries(${LIB_NAME_STATIC} ${SMART_MATH_STATIC_LIBRARY} ${SMART_UQ_STATIC_LIBRARY} ${MANDATORY_LIBRARIES} ${IRILIB_NAME})
    set_target_properties (${LIB_NAME_STATIC} PROPERTIES
            VERSION ${PROJECT_VERSION}
            SOVERSION ${PROJECT_VERSION}
            OUTPUT_NAME ${LIB_NAME})

    if(SMART_MATH_EXTERNAL)
        add_dependencies(${LIB_NAME_STATIC} smart-math)
    endif()

    if(SMART_UQ_EXTERNAL)
        add_dependencies(${LIB_NAME_STATIC} smart-uq)
    endif()
else()
    set(LIB_NAME_STATIC "")
endif ()

add_library(${LIB_NAME} SHARED ${SRC})
target_link_libraries(${LIB_NAME} ${SMART_MATH_LIBRARY} ${SMART_UQ_LIBRARY} ${MANDATORY_LIBRARIES} ${IRILIB_NAME})
set_target_properties (${LIB_NAME} PROPERTIES
        VERSION ${PROJECT_VERSION}
        SOVERSION ${PROJECT_VERSION}
        OUTPUT_NAME ${LIB_NAME})

if(SMART_MATH_EXTERNAL)
    add_dependencies(${LIB_NAME} smart-math)
endif()

if(SMART_UQ_EXTERNAL)
    add_dependencies(${LIB_NAME} smart-uq)
endif()

if(BUILD_MAIN)
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${BIN_PATH})
    add_executable(${MAIN_NAME} ${MAIN_SRC})
    target_link_libraries(${MAIN_NAME} ${LIB_NAME} ${MANDATORY_LIBRARIES})
endif()

if(BUILD_DOCS)
    find_package(Doxygen)

    if(NOT DOXYGEN_FOUND)
        message(FATAL_ERROR
                "Doxygen is needed to build the documentation. Please make sure it is installed.")
    endif()

    find_package(Boost 1.50.0)
    message(STATUS "Detected Boost version: ${Boost_VERSION}")

    configure_file(${PROJECT_PATH}/Doxyfile.in
            ${PROJECT_BINARY_DIR}/Doxyfile @ONLY IMMEDIATE)

    # Add a custom target to run Doxygen whenever the project is built.
    # If you do NOT want the documentation to be generated EVERY time you build the project
    # then leave out the 'ALL' keyword from the command below.
    add_custom_target(docs ALL COMMAND ${DOXYGEN_EXECUTABLE}
            ${PROJECT_BINARY_DIR}/Doxyfile
            SOURCES ${PROJECT_BINARY_DIR}/Doxyfile)
endif(BUILD_DOCS)

if(BUILD_EXAMPLES)
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${EXAMPLES_PATH})
    add_subdirectory("${PROJECT_PATH}/examples")
endif()


# Install header files and library.
# Destination is set by CMAKE_INSTALL_PREFIX and defaults to usual locations, unless overridden by
# user.
install(DIRECTORY ${INCLUDE_PATH}/ DESTINATION include/${PROJECT_NAME} FILES_MATCHING PATTERN "*.h" PATTERN "*.hpp")
install(DIRECTORY ${EXTERNAL_PATH}/ DESTINATION include/${PROJECT_NAME})

install(TARGETS ${LIB_NAME} ${LIB_NAME_STATIC} ${BIN_NAME}
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib/${PROJECT_NAME}
        ARCHIVE DESTINATION lib/${PROJECT_NAME})

#add_custom_target(dist COMMAND ${CMAKE_MAKE_PROGRAM} package)


# Set up packager.
#if(EXISTS "${CMAKE_ROOT}/Modules/CPack.cmake")
#    include(InstallRequiredSystemLibraries)
#
#    set(CPACK_SET_DESTDIR "on")
#    set(CPACK_PACKAGING_INSTALL_PREFIX "/tmp")
#
#    set(CPACK_PACKAGE_NAME "${PROJECT_NAME}")
#    set(CPACK_PACKAGE_DESCRIPTION "${PROJECT_NAME}")
#    set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "${PROJECT_NAME}")
#    set(CPACK_PACKAGE_VENDOR "${PROJECT_NAME}_VENDOR_NAME")
#
#    set(CPACK_RESOURCE_FILE_LICENSE ${CMAKE_SOURCE_DIR}/LICENSE.md)
#    set(CPACK_RESOURCE_FILE_README ${CMAKE_SOURCE_DIR}/README.md)
#    set(CPACK_PACKAGE_VERSION ${PROJECT_VERSION})
#    set(CPACK_PACKAGE_VERSION_MAJOR ${PROJECT_VERSION_MAJOR})
#    set(CPACK_PACKAGE_VERSION_MINOR ${PROJECT_VERSION_MINOR})
#    set(CPACK_PACKAGE_VERSION_PATCH ${PROJECT_VERSION_PATCH})
#
#    set(CPACK_PACKAGE_DESCRIPTION_FILE "${PROJECT_PATH}/README.md")
#    set(CPACK_RESOURCE_FILE_LICENSE "${PROJECT_PATH}/LICENSE.md")
#
#    set(CPACK_PACKAGE_FILE_NAME "${CMAKE_PROJECT_NAME}_${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}.${PROJECT_VERSION_PATCH}")
#    set(CPACK_SOURCE_PACKAGE_FILE_NAME "${CMAKE_PROJECT_NAME}_${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}.${PROJECT_VERSION_PATCH}")
#    set(CPACK_GENERATOR "ZIP")
#
#    # $LINUX_DISTRO Debian Redhat none
#    set(LINUX_DISTRO "none")
#    if(EXISTS "/etc/debian_version")
#        set ( LINUX_DISTRO "Debian")
#    elseif(EXISTS "/etc/redhat-release")
#        set ( LINUX_DISTRO "Redhat")
#    endif()
#
#    if (UNIX AND (LINUX_DISTRO STREQUAL "Debian"))
#        set(CPACK_GENERATOR "DEB")
#        set(CPACK_DEBIAN_PACKAGE_MAINTAINER "strathclydesmart.slack.com")
#        set(CPACK_DEBIAN_PACKAGE_PRIORITY "optional")
#
#        execute_process(COMMAND dpkg --print-architecture OUTPUT_VARIABLE _CPACK_DEBIAN_PACKAGE_ARCHITECTURE)
#        string(STRIP ${_CPACK_DEBIAN_PACKAGE_ARCHITECTURE} _CPACK_DEBIAN_PACKAGE_ARCHITECTURE)
#        set(CPACK_DEBIAN_PACKAGE_ARCHITECTURE ${_CPACK_DEBIAN_PACKAGE_ARCHITECTURE})
#
#        set(CPACK_DEBIAN_REVISON "a1")
#        set(CPACK_PACKAGE_FILE_NAME "${CPACK_PACKAGE_NAME}-${CPACK_PACKAGE_VERSION}-${CPACK_DEBIAN_REVISON}_${_CPACK_DEBIAN_PACKAGE_ARCHITECTURE}")
#        set(CPACK_DEBIAN_PACKAGE_PRIORITY "optional")
#        set(CPACK_DEBIAN_COMPRESSION_TYPE "gzip")
#        set(CPACK_DEBIAN_PACKAGE_DEPENDS "smart-math (>= 0.1.0); smart-uq (>= 1.0.0)")
#
#        mark_as_advanced(_CPACK_DEBIAN_PACKAGE_ARCHITECTURE)
#    elseif (WIN32)
#        # There is a bug in NSI that does not handle full unix paths properly. Make
#        # sure there is at least one set of four (4) backlasshes.
#        set(CPACK_GENERATOR "${CPACK_GENERATOR};NSIS")
#        set(CPACK_NSIS_DISPLAY_NAME "${CPACK_PACKAGE_INSTALL_DIRECTORY}")
#        set(CPACK_NSIS_CONTACT "${PROJECT_NAME}_VENDOR_CONTACT")
#        set(CPACK_NSIS_MODIFY_PATH ON)
#    endif()
#
#    include(CPack)
#endif()
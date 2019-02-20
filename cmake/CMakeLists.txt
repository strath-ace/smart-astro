#project(utopiae-orbit-dependencies
#        LANGUAGES CXX
#        VERSION 1.0.0)

# Include external libraries
include(ExternalProject)

# Smart-math
include(cmake/External/smart-math.cmake)
list (APPEND MANDATORY_LIBRARIES ${SMART_MATH_LIBRARY} ${SMART_MATH_STATIC_LIBRARY})
include_directories ("${SMART_MATH_INCLUDE_DIR}")

# Smart-uq
if (ENABLE_SMARTUQ)
    include(cmake/External/smart-uq.cmake)
    list (APPEND MANDATORY_LIBRARIES ${SMART_UQ_LIBRARY} ${SMART_UQ_STATIC_LIBRARY})
    include_directories (${SMART_UQ_INCLUDE_DIR})
endif()

# Cspice
if (ENABLE_CSPICE)
    # Add compiler flag
    add_definitions(-D__USE_CSPICE)

    include(cmake/External/cspice.cmake)
    list (APPEND MANDATORY_LIBRARIES ${CSPICE_LIBRARY})
    include_directories(${CSPICE_INCLUDE_DIR})
endif()

# OpenMP
find_package(OpenMP)
if(OPENMP_FOUND)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
endif()

#if((CMAKE_COMPILER_IS_GNUCC) AND (CMAKE_BUILD_TYPE STREQUAL "Debug"))
#    # smart-uq requires __gcov_init and __gcov_merge_add symbols
#    set(CMAKE_C_FLAGS  "${CMAKE_C_FLAGS} --coverage")
#    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --coverage")
#    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -lgcov --coverage")
#endif()

# Graphics libraries
if(ENABLE_GRAPHICS)
    set(MANDATORY_LIBRARIES ${MANDATORY_LIBRARIES} "X11" "Xext")
endif(ENABLE_GRAPHICS)

# WORHP
if(ENABLE_WORHP)
    find_library(WORHP_LIBRARY NAMES worhp)
    find_path(WORHP_INCLUDE_DIR NAMES worhp/worhp.h)
    if(NOT WORHP_LIBRARY OR NOT WORHP_INCLUDE_DIR)
        message(FATAL_ERROR "WORHP support was requested, but the worhp library or headers could not be located.")
    endif(NOT WORHP_LIBRARY OR NOT WORHP_INCLUDE_DIR)
    message(STATUS "worhp library: ${WORHP_LIBRARY}")
    message(STATUS "worhp headers: ${WORHP_INCLUDE_DIR}")

    set(MANDATORY_LIBRARIES ${MANDATORY_LIBRARIES} ${WORHP_LIBRARY})
    include_directories("${WORHP_INCLUDE_DIR}")
endif()

#UI

cmake_minimum_required(VERSION 2.8.11)
project (SMART-ASTRO_GUI)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTOUIC ON)

#find_package(Qt5Core REQUIRED)
find_package(Qt5Widgets REQUIRED)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
#include_directories(${Qt5Widgets_INCLUDE_DIRS})
#include_directories(${Qt5Core_INCLUDE_DIRS})

include_directories(${Qt5Widgets_INCLUDES})

add_definitions(${Qt5Widgets_DEFINITIONS})

set(CMAKE_CXX_FLAGS "${Qt5Widgets_EXECUTABLE_COMPILE_FLAGS}")

set(SMART-ASTRO_GUI_SRCS
    src/SMART-ASTRO_GUI/main.cpp
    src/SMART-ASTRO_GUI/mainwindow.cpp
    src/SMART-ASTRO_GUI/genericfunctiondialog.cpp
    src/SMART-ASTRO_GUI/parameterinputdialog.cpp
    src/SMART-ASTRO_GUI/Controller/parameterinputcontroller.cpp
    src/SMART-ASTRO_GUI/answerdialog.cpp
    src/Astro-Core/conversion_coordinates.cpp
    src/Astro-Core/spice_general_functions.cpp
    src/Astro-Core/conversion_time.cpp
)

set(SMART-ASTRO_GUI_HEADERS
    src/SMART-ASTRO_GUI/mainwindow.h
    src/SMART-ASTRO_GUI/genericfunctiondialog.h
    src/SMART-ASTRO_GUI/parameterinputdialog.h
    src/SMART-ASTRO_GUI/Controller/parameterinputcontroller.h
    src/SMART-ASTRO_GUI/answerdialog.h
    src/Astro-Core/conversion_time.h
    src/Astro-Core/conversion_coordinates.h
)

set(SMART-ASTRO_GUI_UIS
   src/SMART-ASTRO_GUI/mainwindow.ui
   src/SMART-ASTRO_GUI/genericfunctiondialog.ui
   src/SMART-ASTRO_GUI/parameterinputdialog.ui
   src/SMART-ASTRO_GUI/answerdialog.ui
)
add_executable(SMART-ASTRO_GUI ${SMART-ASTRO_GUI_SRCS} ${SMART-ASTRO_GUI_HEADERS})

target_link_libraries(SMART-ASTRO_GUI ${CSPICE_LIBRARY} ${Qt5Widgets_LIBRARIES}  )
#Qt5::Widgets Qt5::Core

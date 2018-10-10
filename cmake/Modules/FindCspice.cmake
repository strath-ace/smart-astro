if (NOT CSPICE_ROOT_DIR)
    set(CSPICE_ROOT_DIR "" CACHE PATH "Folder contains CSPICE library")
endif ()

if (CSPICE_ROOT_DIR)
    set(_CSPICE_INCLUDE_LOCATIONS "${CSPICE_ROOT_DIR}/include")
    set(_CSPICE_LIB_LOCATIONS "${CSPICE_ROOT_DIR}/lib")
else ()
    set(_CSPICE_INCLUDE_LOCATIONS "")
    set(_CSPICE_LIB_LOCATIONS "")
endif ()

find_package(PkgConfig QUIET)
if (PKG_CONFIG_FOUND)
    pkg_check_modules(CSPICE_PKGCONF QUIET Cspice)
endif ()

find_path(CSPICE_INCLUDE_DIR cspice/SpiceUsr.h
#        HINTS "${_CSPICE_INCLUDE_LOCATIONS}"
        PATHS "/usr/local/include/cspice/" "/usr/local/include/" "${CSPICE_PKGCONF_INCLUDE_DIRS}"
        )

find_library(CSPICE_LIBRARY cspice.a
#        HINTS "${_CSPICE_LIB_LOCATIONS}"
        PATHS "/usr/local/lib/cspice/" "${CSPICE_PKGCONF_LIBRARY_DIRS}"
#        PATH_SUFFIXES a
        )

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(CSPICE DEFAULT_MSG CSPICE_INCLUDE_DIR CSPICE_LIBRARY)
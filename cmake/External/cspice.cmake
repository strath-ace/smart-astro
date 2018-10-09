if (NOT __CSPICE_INCLUDED) # guard against multiple includes
    set(__CSPICE_INCLUDED TRUE)

    # use the system-wide cspice if present
    find_package(Cspice)

    if (CSPICE_FOUND)
        set(CSPICE_EXTERNAL FALSE)
    else()
        # build directory
        set(cspice_PREFIX ${CMAKE_BINARY_DIR}/external/src/cspice)
        # install directory
        set(cspice_INSTALL ${CMAKE_BINARY_DIR}/external/install/)

        ExternalProject_Add(cspice
                PREFIX ${cspice_PREFIX}
                DOWNLOAD_COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/cmake/External/cspice_download_install.sh
                UPDATE_COMMAND ""
                INSTALL_DIR ${cspice_INSTALL}
                CONFIGURE_COMMAND ""
                BUILD_COMMAND ""
                INSTALL_COMMAND ""
                #CMAKE_ARGS -DCMAKE_BUILD_TYPE=Release
                #-DCMAKE_INSTALL_PREFIX=${cspice_INSTALL}
                #-DCSPICE_ENABLE_TESTING=OFF
                #-DCMAKE_C_FLAGS=${CMAKE_C_FLAGS}
                #-DCMAKE_CXX_FLAGS=${CMAKE_CXX_FLAGS}
                LOG_DOWNLOAD 1)

        set(CSPICE_FOUND TRUE)
        set(CSPICE_ROOT_DIR ${cspice_INSTALL})
        set(CSPICE_INCLUDE_DIR ${cspice_INSTALL}/include/)
        set(CSPICE_INCLUDE_DIRS ${CSPICE_INCLUDE_DIR})
        set(CSPICE_STATIC_LIBRARY ${cspice_INSTALL}/lib/cspice/cspice.a)
        set(CSPICE_EXTERNAL TRUE)
    endif()

endif()
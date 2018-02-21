if (NOT __SMART_UQ_INCLUDED) # guard against multiple includes
    set(__SMART_UQ_INCLUDED TRUE)

    # use the system-wide smart-math if present
    find_package(smart-uq)
    if (SMART_UQ_FOUND)
        set(SMART_UQ_EXTERNAL FALSE)
    else ()
        # build directory
        set(smart-uq_PREFIX ${CMAKE_BINARY_DIR}/external/src/smart-uq)
        # install directory
        set(smart-uq_INSTALL ${CMAKE_BINARY_DIR}/external/install)

        ExternalProject_Add(smart-uq
                PREFIX ${smart-uq_PREFIX}
                INSTALL_DIR ${smart-uq_INSTALL}
                GIT_REPOSITORY "git@github.com:strath-ace-labs/smart-uq.git"
                GIT_TAG "master"
                GIT_SHALLOW 1
                CMAKE_ARGS -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
                -DCMAKE_INSTALL_PREFIX=${smart-math_INSTALL}
                -DCMAKE_C_FLAGS=${GFLAGS_C_FLAGS}
                -DCMAKE_CXX_FLAGS=${GFLAGS_CXX_FLAGS}
                LOG_DOWNLOAD 1
                LOG_INSTALL 1)

        if(SMART_MATH_EXTERNAL)
            add_dependencies(smart-uq smart-math)
        endif()

        set(SMART_UQ_FOUND TRUE)
        set(SMART_UQ_INCLUDE_DIR ${smart-uq_PREFIX}/src/smart-uq/include)
        set(SMART_UQ_LIBRARY "${smart-uq_INSTALL}/lib/smart-uq/libsmart-uq.so")
        set(SMART_UQ_STATIC_LIBRARY ${smart-uq_INSTALL}/lib/smart-uq/libsmart-uq.a)
        set(SMART_UQ_EXTERNAL TRUE)
    endif ()

endif ()

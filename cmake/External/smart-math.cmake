if (NOT __SMART_MATH_INCLUDED) # guard against multiple includes
    set(__SMART_MATH_INCLUDED TRUE)

    # use the system-wide smart-math if present
    find_package(smart-math)
    if (SMART_MATH_FOUND)
        set(SMART_MATH_EXTERNAL FALSE)
    else ()
        # build directory
        set(smart-math_PREFIX ${CMAKE_BINARY_DIR}/external/src/smart-math)
        # install directory
        set(smart-math_INSTALL ${CMAKE_BINARY_DIR}/external/install)

        ExternalProject_Add(smart-math
                PREFIX ${smart-math_PREFIX}
                INSTALL_DIR ${smart-math_INSTALL}
                GIT_REPOSITORY "git@github.com:strath-ace/smart-math.git"
                GIT_TAG "master"
                #GIT_SHALLOW 1
                CMAKE_ARGS -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
                -DCMAKE_INSTALL_PREFIX=${smart-math_INSTALL}
                -DCMAKE_C_FLAGS=${GFLAGS_C_FLAGS}
                -DCMAKE_CXX_FLAGS=${GFLAGS_CXX_FLAGS}
                LOG_DOWNLOAD 1
                LOG_INSTALL 1)

        set(SMART_MATH_FOUND TRUE)
        set(SMART_MATH_INCLUDE_DIR ${smart-math_PREFIX}/src/smart-math/include)
        if (APPLE)
            set(SMART_MATH_LIBRARY ${smart-math_INSTALL}/lib/smart-math/libsmart-math.dylib)
        elseif(UNIX)
            set(SMART_MATH_LIBRARY ${smart-math_INSTALL}/lib/smart-math/libsmart-math.so)
        endif()
        set(SMART_MATH_STATIC_LIBRARY ${smart-math_INSTALL}/lib/smart-math/libsmart-math.a)
        set(SMART_MATH_EXTERNAL TRUE)
    endif ()

endif ()

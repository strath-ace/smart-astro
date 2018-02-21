# ###############################################################################################
# Builds a file from a template only if the file has changed
# Saving a new copy of a template triggers recompilation of its dependencies. It should be done
# only if required to prevent excessive recompilation.
# Usage:
#   build_template("${INCLUDE_PATH}/header_file.h.in" "${INCLUDE_PATH}/header_file.h")
macro(build_template TEMPLATE_FILE OUTPUT_FILE)
    configure_file(${TEMPLATE_FILE} "${OUTPUT_FILE}.temp")
    get_filename_component(_cacheKeyPrefix ${OUTPUT_FILE} NAME_WE)

    file(SHA256 "${OUTPUT_FILE}.temp" _${_cacheKeyPrefix}.sha256)
    if(NOT (${_cacheKeyPrefix}.sha256 EQUAL _${_cacheKeyPrefix}.sha256 AND EXISTS ${OUTPUT_FILE}))
        set(${_cacheKeyPrefix}.sha256 ${_${_cacheKeyPrefix}.sha256} CACHE INTERNAL "Store value of the previous file hash")
        file(RENAME "${OUTPUT_FILE}.temp" ${OUTPUT_FILE})
    endif()
    file(REMOVE "${OUTPUT_FILE}.temp")
    unset(_${_cacheKeyPrefix}.sha256)
    unset(_cacheKeyPrefix)
endmacro()
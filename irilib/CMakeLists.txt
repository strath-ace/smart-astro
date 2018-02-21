# build wrapper and link to fort lib

set(IRILIBSRCs
    "IRIw"
)
add_library(IRILIB "${IRILIBSRCs}")

set(IRILIBfortSRCs
    "fortran_source/irisub.for"
    "fortran_source/irifun.for"
    "fortran_source/iriflip.for"
    "fortran_source/iritec.for"
    "fortran_source/iridreg.for"
    "fortran_source/igrf.for"
    "fortran_source/cira.for"
)


set(CMAKE_REQUIRED_LIBRARIES gfortran)
enable_language(Fortran)
# build fort lib
add_library(IRILIBfort STATIC ${IRILIBfortSRCs})
set_target_properties(IRILIBfort PROPERTIES LINKER_LANGUAGE Fortran)

target_link_libraries(IRILIB IRILIBfort gfortran)

install(FILES "IRI_wrapper.h" DESTINATION "include")
install(TARGETS IRILIB DESTINATION "lib")

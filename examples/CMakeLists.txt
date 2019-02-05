
if(ENABLE_CSPICE)
    file(COPY "${SRC_PATH}/AstroData/SpiceKernels/lsk/" DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/SpiceKernels/lsk/")
    file(COPY "${SRC_PATH}/AstroData/SpiceKernels/spk/" DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/SpiceKernels/spk/")
    file(COPY "${SRC_PATH}/AstroData/SpiceKernels/fk/" DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/SpiceKernels/fk/")
    file(COPY "${SRC_PATH}/AstroData/SpiceKernels/pck/" DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/SpiceKernels/pck/")
endif()

add_executable(example_cspice example_cspice.cpp)
target_link_libraries(example_cspice ${LIB_NAME} ${MANDATORY_LIBRARIES})
#
#add_executable(example_OPPED_IPA_3dof example_OPPED_IPA_3dof.cpp)
#target_link_libraries(example_OPPED_IPA_3dof ${LIB_NAME} ${MANDATORY_LIBRARIES})
#
#add_executable(example_OPPED_IPA_6dof example_OPPED_IPA_6dof.cpp)
#target_link_libraries(example_OPPED_IPA_6dof ${LIB_NAME} ${MANDATORY_LIBRARIES})
#
#add_executable(example_OPPED_IPA_reentry example_OPPED_IPA_reentry.cpp)
#target_link_libraries(example_OPPED_IPA_reentry ${LIB_NAME} ${MANDATORY_LIBRARIES})
#
#add_executable(example_ephemerides example_ephemerides.cpp)
#target_link_libraries(example_ephemerides ${LIB_NAME} ${MANDATORY_LIBRARIES})
#
#add_executable(example_charged_MEO example_charged_MEO.cpp)
#target_link_libraries(example_charged_MEO ${LIB_NAME} ${MANDATORY_LIBRARIES})
#
#add_executable(example_charged_3dof example_charged_3dof.cpp)
#target_link_libraries(example_charged_3dof ${LIB_NAME} ${MANDATORY_LIBRARIES})
#
#add_executable(example_charged_6dof example_charged_6dof.cpp)
#target_link_libraries(example_charged_6dof ${LIB_NAME} ${MANDATORY_LIBRARIES})
#
#add_executable(example_dearth_orb example_dearth_orb.cpp)
#target_link_libraries(example_dearth_orb ${LIB_NAME} ${MANDATORY_LIBRARIES})
#
#add_executable(example_dearth_analytical example_dearth_analytical.cpp)
#target_link_libraries(example_dearth_analytical ${LIB_NAME} ${MANDATORY_LIBRARIES})
#
#add_executable(example_dearth_semianalytical example_dearth_semianalytical.cpp)
#target_link_libraries(example_dearth_semianalytical ${LIB_NAME} ${MANDATORY_LIBRARIES})
#
#add_executable(example_hamiltonian example_hamiltonian.cpp)
#target_link_libraries(example_hamiltonian ${LIB_NAME} ${MANDATORY_LIBRARIES})
#
#add_executable(example_magnetic_semianalytical example_magnetic_semianalytical.cpp)
#target_link_libraries(example_magnetic_semianalytical ${LIB_NAME} ${MANDATORY_LIBRARIES})

if(ENABLE_SMARTUQ)

#add_executable(example_OPPED_IPA_nonintrusive example_OPPED_IPA_nonintrusive.cpp)
#target_link_libraries(example_OPPED_IPA_nonintrusive ${LIB_NAME} ${MANDATORY_LIBRARIES})
#
#add_executable(example_OPPED_IPA_cheby example_OPPED_IPA_cheby.cpp)
#target_link_libraries(example_OPPED_IPA_cheby ${LIB_NAME} ${MANDATORY_LIBRARIES})
#
#add_executable(example_OPPED_IPA_taylor example_OPPED_IPA_taylor.cpp)
#target_link_libraries(example_OPPED_IPA_taylor ${LIB_NAME} ${MANDATORY_LIBRARIES})

endif(ENABLE_SMARTUQ)



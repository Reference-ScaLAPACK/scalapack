message("Running BLACS TESTS")
message(STATUS "${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ./${TEST_PROG}")
message(STATUS "Output out_${TEST_PROG}.txt")
file(COPY ${RUNTIMEDIR}/${TEST_PROG} DESTINATION ${OUTPUTDIR})

execute_process(COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ./${TEST_PROG}
                OUTPUT_FILE "out_${TEST_PROG}.txt"
                ERROR_FILE "error_${TEST_PROG}.txt"
                RESULT_VARIABLE HAD_ERROR)

if(HAD_ERROR)
    # This is normal to exit in Error (good behaviour)
    # So we are going to check that the output have the last line of the testing : DONE BLACS_GRIDEXIT
    file(READ "out_${TEST_PROG}.txt" TESTSTRING)

    STRING(REPLACE "DONE BLACS_GRIDEXIT" "BLACS OK" tmp ${TESTSTRING})

if("${tmp}" STREQUAL "${TESTSTRING}")
       message( STATUS "Error in error_${TEST_PROG}.txt")
       message(FATAL_ERROR "Test failed - Test did not reach DONE BLACS_GRIDEXIT")
else()
       message( STATUS "Test Passed")
    endif()
endif()

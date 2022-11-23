get_filename_component(OUTPUTDIR ${EXE} DIRECTORY)
get_filename_component(name ${EXE} NAME_WE)

set(cmd ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS} ${EXE})
set(outfile ${OUTPUTDIR}/out_${name}.txt)
set(errfile ${OUTPUTDIR}/err_${name}.txt)
message(STATUS "Output/Error log files: ${outfile} ${errfile}")
set(CMAKE_EXECUTE_PROCESS_COMMAND_ECHO STDOUT)

execute_process(COMMAND ${cmd}
OUTPUT_FILE ${outfile}
ERROR_FILE ${errfile}
RESULT_VARIABLE HAD_ERROR
WORKING_DIRECTORY ${SOURCE_DIR}
)

if(NOT HAD_ERROR EQUAL 0)
    # This is normal to exit in Error (good behaviour)
    # So we are going to check that the output have the last line of the testing : DONE BLACS_GRIDEXIT
    file(READ ${outfile} TESTSTRING)

    STRING(REPLACE "DONE BLACS_GRIDEXIT" "BLACS OK" tmp "${TESTSTRING}")

if("${tmp}" STREQUAL "${TESTSTRING}")
       message( STATUS "Error in ${errfile}")
       message(FATAL_ERROR "Test failed - Test did not reach DONE BLACS_GRIDEXIT")
else()
       message( STATUS "Test Passed")
    endif()
endif()

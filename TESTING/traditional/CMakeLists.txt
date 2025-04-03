foreach(file_dat IN ITEMS
    BLLT.dat
    BLU.dat
    BRD.dat
    EVC.dat
    HRD.dat
    INV.dat
    LLT.dat
    LS.dat
    LU.dat
    NEP.dat
    QR.dat
    SEP.dat
    SEPR.dat
    SVD.dat
    TRD.dat
)
    configure_file(${file_dat} ${file_dat} COPYONLY)
endforeach()

set(traditional_test_workdir ${CMAKE_CURRENT_BINARY_DIR})
function(add_traditional_test test num_proc)
    add_test(NAME ${test}
        COMMAND
            ${MPIEXEC}
                ${MPIEXEC_NUMPROC_FLAG} ${num_proc}
                ${MPIEXEC_PREFLAGS}
                $<TARGET_FILE:${test}>
        WORKING_DIRECTORY ${traditional_test_workdir}
    )
    set_tests_properties(${test} PROPERTIES
        PROCESSORS ${num_proc}
    )
endfunction()

add_subdirectory(LIN)
add_subdirectory(EIG)


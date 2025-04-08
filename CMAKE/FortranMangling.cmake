# Macro that defines variables describing the Fortran name mangling
# convention
#
# Sets the following outputs on success:
#
#  INTFACE
#    Add_
#    NoChange
#    f77IsF2C
#    UpCase
#

include_guard()

block()
    # TODO: This path is hard-coded
    set(BLACS_INSTALL_SRC
        ${CMAKE_CURRENT_LIST_DIR}/../BLACS/INSTALL
    )

    try_run(xintface_res xintface_compile_res
        SOURCES
            ${BLACS_INSTALL_SRC}/Fintface.f
            ${BLACS_INSTALL_SRC}/Cintface.c
        NO_CACHE
        COMPILE_OUTPUT_VARIABLE xintface_compile_output
        RUN_OUTPUT_VARIABLE xintface_output
    )
    if(NOT xintface_compile_res)
        message(FATAL_ERROR
            "Could not compile BLACS/INSTALL:\n"
            "${xintface_compile_output}"
        )
    endif()
    if(NOT ${xintface_res} EQUAL 0)
        message(FATAL_ERROR
            "xintface did not execute properly:\n"
            "${xintface_output}"
        )
    endif()
    string(STRIP "${xintface_output}" xintface_output)
    set(CDEFS ${xintface_output} CACHE STRING "Fortran Mangling")
endblock()

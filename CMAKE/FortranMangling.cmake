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

FUNCTION(COMPILE RESULT)
    MESSAGE(STATUS "=========")
    MESSAGE(STATUS "Compiling and Building BLACS INSTALL Testing to set correct variables")
   
   # Configure: 
    EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND}  
         "-DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}"
         "-DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}" "${PROJECT_SOURCE_DIR}//BLACS/INSTALL"
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/BLACS/INSTALL/        
        RESULT_VARIABLE RESVAR OUTPUT_VARIABLE LOG1_OUT ERROR_VARIABLE LOG1_ERR
    )
    if(RESVAR EQUAL 0)
    MESSAGE(STATUS "Configure in the INSTALL directory successful")
    else()
    MESSAGE(STATUS " Build Output:\n ${LOG1_OUT}")
    MESSAGE(STATUS " Error Output:\n ${LOG1_ERR}")
    MESSAGE(FATAL_ERROR " Configure in the BLACS INSTALL directory FAILED")
    endif()

    # Build:
    EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND} --build
        ${PROJECT_SOURCE_DIR}/BLACS/INSTALL/ 
        RESULT_VARIABLE RESVAR OUTPUT_VARIABLE LOG2_OUT ERROR_VARIABLE LOG2_ERR
    )
    if(RESVAR  EQUAL 0)
    MESSAGE(STATUS "Build in the BLACS INSTALL directory successful")
    else()
    MESSAGE(STATUS " Build Output:\n ${LOG2_OUT}")
    MESSAGE(STATUS " Error Output:\n ${LOG2_ERR}")
    MESSAGE(FATAL_ERROR " Build in the BLACS INSTALL directory FAILED")
    endif()
    # Clean up:
    FILE(REMOVE_RECURSE ${PROJECT_SOURCE_DIR}/BLACS/INSTALL/CMakeCache.txt)
    FILE(REMOVE_RECURSE ${PROJECT_SOURCE_DIR}/BLACS/INSTALL/CMakeFiles )
ENDFUNCTION()


macro(FORTRAN_MANGLING CDEFS)
MESSAGE(STATUS "=========")
MESSAGE(STATUS "Testing FORTRAN_MANGLING")
   
    execute_process ( COMMAND  ${PROJECT_SOURCE_DIR}/BLACS/INSTALL/xintface
                         RESULT_VARIABLE xintface_RES
                         OUTPUT_VARIABLE xintface_OUT
                         ERROR_VARIABLE xintface_ERR)
                         

#    MESSAGE(STATUS "FORTRAN MANGLING:RUN \n${xintface_OUT}")

       if (xintface_RES EQUAL 0)
          STRING(REPLACE "\n" "" xintface_OUT "${xintface_OUT}")
          MESSAGE(STATUS "CDEFS set to ${xintface_OUT}")
          SET(CDEFS ${xintface_OUT} CACHE STRING "Fortran Mangling" FORCE)
      else()
          MESSAGE(STATUS " xintface Output:\n ${xintface_OUT}")
          MESSAGE(FATAL_ERROR "FORTRAN_MANGLING:ERROR ${xintface_ERR}")
      endif() 
endmacro(FORTRAN_MANGLING)

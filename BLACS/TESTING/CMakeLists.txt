if(NOT (BUILD_SINGLE AND BUILD_DOUBLE AND BUILD_COMPLEX AND BUILD_COMPLEX16))
  return()
endif()

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

set(FTestObj
   blacstest.f btprim.f tools.f)

if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU" AND CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER_EQUAL 10)
  set_property(SOURCE blacstest.f PROPERTY COMPILE_FLAGS "-std=legacy")
endif()

add_executable(xFbtest ${FTestObj})
target_link_libraries(xFbtest PRIVATE scalapack MPI::MPI_Fortran)

set(CTestObj
   Cbt.c)

set_property(
   SOURCE Cbt.c
   APPEND PROPERTY COMPILE_DEFINITIONS BTCINTFACE
   )

add_executable(xCbtest ${CTestObj} ${FTestObj})
target_link_libraries(xCbtest PRIVATE scalapack MPI::MPI_Fortran)

# We could run the BLACS TESTING the following way
# But BLACS TESTING are TESTING anormal exit so even if they pass,
# CTest will determine they fail
#add_test(xFbtest0 ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS} ./xFbtest)
#add_test(xCbtest0 ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS} ./xCbtest)

add_test(NAME xCbtest
  COMMAND ${CMAKE_COMMAND}
    -DMPIEXEC_EXECUTABLE=${MPIEXEC_EXECUTABLE}
    -DMPIEXEC_NUMPROC_FLAG=${MPIEXEC_NUMPROC_FLAG}
    -DMPIEXEC_PREFLAGS=${MPIEXEC_PREFLAGS}
    -DSOURCE_DIR=${CMAKE_CURRENT_SOURCE_DIR}
    -DEXE=$<TARGET_FILE:xCbtest>
    -P ${CMAKE_CURRENT_SOURCE_DIR}/runtest.cmake
    )

add_test(NAME xFbtest
  COMMAND ${CMAKE_COMMAND}
    -DMPIEXEC_EXECUTABLE=${MPIEXEC_EXECUTABLE}
    -DMPIEXEC_NUMPROC_FLAG=${MPIEXEC_NUMPROC_FLAG}
    -DMPIEXEC_PREFLAGS=${MPIEXEC_PREFLAGS}
    -DSOURCE_DIR=${CMAKE_CURRENT_SOURCE_DIR}
    -DEXE=$<TARGET_FILE:xFbtest>
    -P ${CMAKE_CURRENT_SOURCE_DIR}/runtest.cmake
    )

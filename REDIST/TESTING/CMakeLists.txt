set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

function(_redist_test name src)
  add_executable(${name} ${src})
  target_link_libraries(${name} PRIVATE scalapack MPI::MPI_Fortran)
  # add_test(NAME ${name}
  # COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS} $<TARGET_FILE:${name}>
  # WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  # )
endfunction()

_redist_test(xigemr pigemrdrv.c)
_redist_test(xitrmr pitrmrdrv.c)

if(BUILD_SINGLE)
  _redist_test(xsgemr psgemrdrv.c)
  _redist_test(xstrmr pstrmrdrv.c)
endif()
if(BUILD_DOUBLE)
  _redist_test(xdgemr pdgemrdrv.c)
  _redist_test(xdtrmr pdtrmrdrv.c)
endif()
if(BUILD_COMPLEX)
  _redist_test(xcgemr pcgemrdrv.c)
  _redist_test(xctrmr pctrmrdrv.c)
endif()
if(BUILD_COMPLEX16)
  _redist_test(xzgemr pzgemrdrv.c)
  _redist_test(xztrmr pztrmrdrv.c)
endif()

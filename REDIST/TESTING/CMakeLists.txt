set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${SCALAPACK_BINARY_DIR}/scalapack/REDIST/TESTING)

file(COPY GEMR2D.dat DESTINATION ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
file(COPY TRMR2D.dat DESTINATION ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})

add_executable(xigemr pigemrdrv.c)
add_executable(xsgemr psgemrdrv.c)
add_executable(xdgemr pdgemrdrv.c)
add_executable(xcgemr pcgemrdrv.c)
add_executable(xzgemr pzgemrdrv.c)

add_executable(xitrmr pitrmrdrv.c)
add_executable(xstrmr pstrmrdrv.c)
add_executable(xdtrmr pdtrmrdrv.c)
add_executable(xctrmr pctrmrdrv.c)
add_executable(xztrmr pztrmrdrv.c)

target_link_libraries(xigemr scalapack MPI::MPI_Fortran )
target_link_libraries(xsgemr scalapack MPI::MPI_Fortran )
target_link_libraries(xdgemr scalapack MPI::MPI_Fortran )
target_link_libraries(xcgemr scalapack MPI::MPI_Fortran )
target_link_libraries(xzgemr scalapack MPI::MPI_Fortran )

target_link_libraries(xitrmr scalapack MPI::MPI_Fortran )
target_link_libraries(xstrmr scalapack MPI::MPI_Fortran )
target_link_libraries(xdtrmr scalapack MPI::MPI_Fortran )
target_link_libraries(xctrmr scalapack MPI::MPI_Fortran )
target_link_libraries(xztrmr scalapack MPI::MPI_Fortran )

#add_test(NAME xigemr COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS} ./xigemr COMMAND_EXPAND_LISTS)
#add_test(NAME xsgemr COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS} ./xsgemr COMMAND_EXPAND_LISTS)
#add_test(NAME xdgemr COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS} ./xdgemr COMMAND_EXPAND_LISTS)
#add_test(NAME xcgemr COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS} ./xcgemr COMMAND_EXPAND_LISTS)
#add_test(NAME xzgemr COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS} ./xzgemr COMMAND_EXPAND_LISTS)

#add_test(NAME xitrmr COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS} ./xitrmr COMMAND_EXPAND_LISTS)
#add_test(NAME xstrmr COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS} ./xstrmr COMMAND_EXPAND_LISTS)
#add_test(NAME xdtrmr COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS} ./xdtrmr COMMAND_EXPAND_LISTS)
#add_test(NAME xctrmr COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS} ./xctrmr COMMAND_EXPAND_LISTS)
#add_test(NAME xztrmr COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS} ./xztrmr COMMAND_EXPAND_LISTS)


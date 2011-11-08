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

target_link_libraries(xigemr scalapack )
target_link_libraries(xsgemr scalapack )
target_link_libraries(xdgemr scalapack )
target_link_libraries(xcgemr scalapack )
target_link_libraries(xzgemr scalapack )

target_link_libraries(xitrmr scalapack )
target_link_libraries(xstrmr scalapack )
target_link_libraries(xdtrmr scalapack )
target_link_libraries(xctrmr scalapack )
target_link_libraries(xztrmr scalapack )

#add_test(xigemr ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ./xigemr)
#add_test(xsgemr ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ./xsgemr)
#add_test(xdgemr ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ./xdgemr)
#add_test(xcgemr ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ./xcgemr)
#add_test(xzgemr ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ./xzgemr)

#add_test(xitrmr ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ./xitrmr)
#add_test(xstrmr ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ./xstrmr)
#add_test(xdtrmr ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ./xdtrmr)
#add_test(xctrmr ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ./xctrmr)
#add_test(xztrmr ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ./xztrmr)


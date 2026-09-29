add_subdirectory(SRC)
if(BUILD_TESTING AND SCALAPACK_BUILD_SINGLE)
  add_library(pblas-test-support OBJECT
    SRC/PTOOLS/PB_Cwarn.c
    SRC/PTOOLS/PB_Cabort.c)
  target_compile_definitions(pblas-test-support PRIVATE TestingPblas)
  # TestingPblas headers are referenced relative to PBLAS/TESTING.
  target_include_directories(pblas-test-support PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}/TESTING")

  add_subdirectory(TESTING)
  add_subdirectory(TIMING)
endif()

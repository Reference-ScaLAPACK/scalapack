add_subdirectory(SRC)
if(${SCALAPACK_BUILD_TESTING})
  add_subdirectory(TESTING)
endif()

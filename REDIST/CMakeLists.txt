add_subdirectory(SRC)
if(${SCALAPACK_BUILD_TESTS})
  add_subdirectory(TESTING)
endif()

cmake_minimum_required(VERSION 3.26...4.0)

## Old test-suite to be modernized
# include only if this is included via `add_subdirectory`

if(DEFINED SCALAPACK_SOURCE_DIR)
    add_subdirectory(traditional)
endif()

## Modern test-suite is an independent project

# TODO: Add the implementations
project(scalapack-tests
    LANGUAGES C Fortran
)
enable_testing()

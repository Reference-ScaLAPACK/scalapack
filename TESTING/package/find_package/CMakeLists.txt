cmake_minimum_required(VERSION 3.26...4.0)

project(test_find_package
    LANGUAGES C Fortran
)

find_package(scalapack REQUIRED CONFIG)

# TODO: Add a proper smoke test
if(NOT TARGET scalapack::scalapack)
    message(FATAL_ERROR "scalapack::scalapack target not present")
endif()

enable_testing()

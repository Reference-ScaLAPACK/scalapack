cmake_minimum_required(VERSION 3.26...4.0)

project(scalapack-tests
    LANGUAGES C Fortran
)

if(NOT TARGET scalapack::scalapack)
    find_package(scalapack REQUIRED CONFIG)
endif()

enable_testing()

# Old test-suite
add_subdirectory(traditional)

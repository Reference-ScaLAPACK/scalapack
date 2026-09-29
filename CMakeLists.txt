cmake_minimum_required(VERSION 3.26...4.0)

project(
  SCALAPACK
  VERSION 2.2.3
  LANGUAGES C Fortran)

list(PREPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/CMAKE")

include(GNUInstallDirs)
include(CMakeDependentOption)
include(CMakePackageConfigHelpers)

set(CMAKE_C_STANDARD 17)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_EXTENSIONS OFF)

# ==============================================================================
# Deprecated option compatibility (TODO: Remove this in a future release)
# ==============================================================================
function(scalapack_deprecate_build_option suffix description)
  set(new_option "SCALAPACK_BUILD_${suffix}")
  set(old_option "BUILD_${suffix}")
  if(DEFINED ${old_option})
    message(WARNING "${old_option} is deprecated; use ${new_option} instead")
    if(NOT DEFINED ${new_option})
      set(${new_option}
          "${${old_option}}"
          CACHE BOOL "${description}")
    endif()
  endif()
endfunction()

scalapack_deprecate_build_option(SINGLE "Build single precision real")
scalapack_deprecate_build_option(DOUBLE "Build double precision real")
scalapack_deprecate_build_option(COMPLEX "Build single precision complex")
scalapack_deprecate_build_option(COMPLEX16 "Build double precision complex")
# ==============================================================================

# Build options
option(BUILD_SHARED_LIBS "Build shared libraries" OFF)
option(BUILD_TESTING "Build tests of the library" ON)

option(SCALAPACK_BUILD_SINGLE "Build single precision real" ON)
option(SCALAPACK_BUILD_DOUBLE "Build double precision real" ON)
cmake_dependent_option(SCALAPACK_BUILD_COMPLEX "Build single precision complex"
                       ON "SCALAPACK_BUILD_SINGLE" OFF)
cmake_dependent_option(SCALAPACK_BUILD_COMPLEX16 "Build double precision complex"
                       ON "SCALAPACK_BUILD_DOUBLE" OFF)

if(UNIX AND CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
  string(APPEND CMAKE_Fortran_FLAGS " -fltconsistency -fp_port")
endif()

# Set default build type
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release")
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release"
                                               "RelWithDebInfo")
endif()

set(DART_TESTING_TIMEOUT 600)

# Keep libraries together and test executables under TESTING.
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/TESTING")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/lib")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/lib")

message(STATUS "Build single precision real: ${SCALAPACK_BUILD_SINGLE}")
message(STATUS "Build double precision real: ${SCALAPACK_BUILD_DOUBLE}")
message(STATUS "Build single precision complex: ${SCALAPACK_BUILD_COMPLEX}")
message(STATUS "Build double precision complex: ${SCALAPACK_BUILD_COMPLEX16}")

# Platform/compiler checks used by the BLACS implementation.
include(CheckBLACSCompilerFlags)
checkblacscompilerflags()
include(FortranMangling)
add_compile_definitions("${CDEFS}")

# External dependencies.  The library itself uses the MPI C interface; the
# Fortran MPI interface is only needed by the test suite (and the MSVC split).
set(_scalapack_mpi_components C)
if(BUILD_TESTING OR MSVC)
  list(APPEND _scalapack_mpi_components Fortran)
endif()
find_package(MPI REQUIRED COMPONENTS ${_scalapack_mpi_components})
find_package(LAPACK REQUIRED)

# Canonical library targets.  MSVC keeps C and Fortran code in separate
# libraries, matching the existing build layout.
add_library(scalapack)
add_library(ScaLAPACK::scalapack ALIAS scalapack)
set(SCALAPACK_FORTRAN_TARGET scalapack)
set(_scalapack_library_targets scalapack)

if(MSVC)
  add_library(scalapack-F)
  add_library(ScaLAPACK::scalapack-F ALIAS scalapack-F)
  set(SCALAPACK_FORTRAN_TARGET scalapack-F)
  list(APPEND _scalapack_library_targets scalapack-F)
endif()

foreach(_scalapack_target IN LISTS _scalapack_library_targets)
  set_target_properties(
    "${_scalapack_target}"
    PROPERTIES POSITION_INDEPENDENT_CODE ON
               VERSION "${PROJECT_VERSION}"
               SOVERSION "${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}")
endforeach()

target_link_libraries(scalapack PRIVATE LAPACK::LAPACK MPI::MPI_C)
if(MSVC)
  target_link_libraries(scalapack-F PRIVATE LAPACK::LAPACK MPI::MPI_Fortran)
endif()

# Each source directory contributes directly to the canonical target(s).
add_subdirectory(BLACS)
add_subdirectory(TOOLS)
add_subdirectory(PBLAS)
add_subdirectory(REDIST)
add_subdirectory(SRC)

if(BUILD_TESTING)
  enable_testing()
  add_subdirectory(TESTING)
endif()

# --------------------------------------------------
# CPACK Packaging

set(CPACK_PACKAGE_NAME "ScaLAPACK")
set(CPACK_PACKAGE_VENDOR
    "University of Tennessee, Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd"
)
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "ScaLAPACK- Linear Algebra Package")
set(CPACK_PACKAGE_VERSION_MAJOR ${SCALAPACK_VERSION_MAJOR})
set(CPACK_PACKAGE_VERSION_MINOR ${SCALAPACK_VERSION_MINOR})
set(CPACK_PACKAGE_VERSION_PATCH ${SCALAPACK_VERSION_PATCH})
set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/LICENSE")
set(CPACK_PACKAGE_INSTALL_DIRECTORY "SCALAPACK")
if(WIN32 AND NOT UNIX)
  # There is a bug in NSI that does not handle full unix paths properly. Make
  # sure there is at least one set of four (4) backlasshes.
  set(CPACK_NSIS_HELP_LINK "http:\\\\\\\\http://icl.cs.utk.edu/lapack-forum")
  set(CPACK_NSIS_URL_INFO_ABOUT "http:\\\\\\\\www.netlib.org/scalapack")
  set(CPACK_NSIS_CONTACT "scalapack@eecs.utk.edu")
  set(CPACK_NSIS_MODIFY_PATH ON)
  set(CPACK_NSIS_DISPLAY_NAME "SCALAPACK-${SCALAPACK_VERSION}")
  set(CPACK_PACKAGE_RELOCATABLE "true")
else(WIN32 AND NOT UNIX)
  set(CPACK_GENERATOR "TGZ")
  set(CPACK_SOURCE_GENERATOR TGZ)
  set(CPACK_SOURCE_PACKAGE_FILE_NAME "scalapack-${SCALAPACK_VERSION}")
  set(CPACK_SOURCE_IGNORE_FILES ~$ .svn ${CPACK_SOURCE_IGNORE_FILES})
endif(WIN32 AND NOT UNIX)
include(CPack)

# --------------------------------------------------
# Installation and package metadata

set(SCALAPACK_INSTALL_CMAKEDIR "${CMAKE_INSTALL_LIBDIR}/cmake/scalapack")

install(TARGETS ${_scalapack_library_targets} EXPORT scalapack-targets)

configure_package_config_file(
  CMAKE/scalapack-config.cmake.in
  "${CMAKE_CURRENT_BINARY_DIR}/scalapack-config.cmake"
  INSTALL_DESTINATION "${SCALAPACK_INSTALL_CMAKEDIR}")
write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/scalapack-config-version.cmake"
  VERSION "${PROJECT_VERSION}"
  COMPATIBILITY SameMinorVersion)

install(
  EXPORT scalapack-targets
  FILE scalapack-targets.cmake
  NAMESPACE ScaLAPACK::
  DESTINATION "${SCALAPACK_INSTALL_CMAKEDIR}")
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/scalapack-config.cmake"
              "${CMAKE_CURRENT_BINARY_DIR}/scalapack-config-version.cmake"
        DESTINATION "${SCALAPACK_INSTALL_CMAKEDIR}")

if(IS_ABSOLUTE "${CMAKE_INSTALL_LIBDIR}")
  set(_scalapack_pc_libdir "${CMAKE_INSTALL_LIBDIR}")
else()
  set(_scalapack_pc_libdir "\${prefix}/${CMAKE_INSTALL_LIBDIR}")
endif()
configure_file(scalapack.pc.in scalapack.pc @ONLY)
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/scalapack.pc"
        DESTINATION "${CMAKE_INSTALL_LIBDIR}/pkgconfig")

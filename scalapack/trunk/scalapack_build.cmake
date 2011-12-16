cmake_minimum_required(VERSION 2.8)
###################################################################
# The values in this section must always be provided
###################################################################
if(UNIX)
  if(NOT compiler)
    set(compiler gcc)
  endif(NOT compiler)
  if(NOT c_compiler)
    set(c_compiler gcc)
  endif(NOT c_compiler)
  if(NOT full_compiler)
    set(full_compiler g++)
  endif(NOT full_compiler)
endif(UNIX)

if(EXISTS "/proc/cpuinfo")
  set(parallel 1)
  file(STRINGS "/proc/cpuinfo" CPUINFO)
  foreach(line ${CPUINFO})
    if("${line}" MATCHES processor)
      math(EXPR parallel "${parallel} + 1")
    endif()
  endforeach(line)
endif()

if(WIN32)
  set(VSLOCATIONS 
    "[HKEY_LOCAL_MACHINE\\SOFTWARE\\Microsoft\\VisualStudio\\6.0\\Setup;VsCommonDir]/MSDev98/Bin"
    "[HKEY_LOCAL_MACHINE\\SOFTWARE\\Microsoft\\VisualStudio\\7.0\\Setup\\VS;EnvironmentDirectory]"
    "[HKEY_LOCAL_MACHINE\\SOFTWARE\\Microsoft\\VisualStudio\\7.1\\Setup\\VS;EnvironmentDirectory]"
    "[HKEY_LOCAL_MACHINE\\SOFTWARE\\Microsoft\\VisualStudio\\8.0;InstallDir]"
    "[HKEY_LOCAL_MACHINE\\SOFTWARE\\Microsoft\\VisualStudio\\8.0\\Setup;Dbghelp_path]"
    "[HKEY_LOCAL_MACHINE\\SOFTWARE\\Microsoft\\VisualStudio\\9.0\\Setup\\VS;EnvironmentDirectory]"
    )
  set(GENERATORS
    "Visual Studio 6"
    "Visual Studio 7"
    "Visual Studio 7 .NET 2003"
    "Visual Studio 8 2005"
    "Visual Studio 8 2005"
    "Visual Studio 9 2008"
    )
  set(vstype 0)
  foreach(p ${VSLOCATIONS})
    get_filename_component(VSPATH ${p} PATH)
    if(NOT "${VSPATH}" STREQUAL "/" AND EXISTS "${VSPATH}")
      message(" found VS install = ${VSPATH}")
      set(genIndex ${vstype})
    endif()
    math(EXPR vstype "${vstype} +1")
  endforeach()
  if(NOT DEFINED genIndex)
    message(FATAL_ERROR "Could not find installed visual stuido")
  endif()
  list(GET GENERATORS ${genIndex} GENERATOR)
  set(CTEST_CMAKE_GENERATOR      "${GENERATOR}")
  message("${CTEST_CMAKE_GENERATOR} - found")
  set(compiler cl)
endif(WIN32)

find_program(HOSTNAME NAMES hostname)
find_program(UNAME NAMES uname)

# Get the build name and hostname
exec_program(${HOSTNAME} ARGS OUTPUT_VARIABLE hostname)
string(REGEX REPLACE "[/\\\\+<> #]" "-" hostname "${hostname}")

message("HOSTNAME: ${hostname}")
# default to parallel 1
if(NOT DEFINED parallel)
  set(parallel 1)
endif(NOT DEFINED parallel)

# find SVN
find_program(SVN svn PATHS $ENV{HOME}/bin /vol/local/bin)
if(NOT SVN)
  message(FATAL_ERROR "SVN not found")
endif()

set(CTEST_UPDATE_COMMAND       ${SVN})
macro(getuname name flag)
  exec_program("${UNAME}" ARGS "${flag}" OUTPUT_VARIABLE "${name}")
  string(REGEX REPLACE "[/\\\\+<> #]" "-" "${name}" "${${name}}")
  string(REGEX REPLACE "^(......|.....|....|...|..|.).*" "\\1" "${name}" "${${name}}")
endmacro(getuname)

getuname(osname -s)
getuname(osver  -v)
getuname(osrel  -r)
getuname(cpu    -m)
if("${osname}" MATCHES Darwin)
  find_program(SW_VER sw_vers)
  execute_process(COMMAND "${SW_VER}" -productVersion OUTPUT_VARIABLE osver)
  string(REPLACE "\n" "" osver "${osver}")
  set(osname "MacOSX")
  set(osrel "")
  if("${cpu}" MATCHES "Power")
    set(cpu "ppc")
  endif("${cpu}" MATCHES "Power")
endif("${osname}" MATCHES Darwin)

if(NOT compiler)
  message(FATAL_ERROR "compiler must be set")
endif(NOT compiler)

  
set(BUILDNAME "${osname}${osver}${osrel}${cpu}-${compiler}")
message("BUILDNAME: ${BUILDNAME}")

# this is the cvs module name that should be checked out
set (CTEST_MODULE_NAME scalapack)
set (CTEST_DIR_NAME "${CTEST_MODULE_NAME}SVN")

# Settings:
message("NOSPACES = ${NOSPACES}")
if(NOSPACES)
  set(CTEST_DASHBOARD_ROOT    "$ENV{HOME}/Dashboards/MyTests-${BUILDNAME}")
else(NOSPACES)
  set(CTEST_DASHBOARD_ROOT    "$ENV{HOME}/Dashboards/My Tests-${BUILDNAME}")
endif(NOSPACES)
set(CTEST_SITE              "${hostname}")
set(CTEST_BUILD_NAME        "${BUILDNAME}")
set(CTEST_TEST_TIMEOUT           "600")

# CVS command and the checkout command
if(NOT EXISTS "${CTEST_DASHBOARD_ROOT}/${CTEST_DIR_NAME}")
  set(CTEST_CHECKOUT_COMMAND     
    "\"${CTEST_UPDATE_COMMAND}\" co https://icl.cs.utk.edu/svn/scalapack-dev/scalapack/trunk ${CTEST_DIR_NAME}")
endif(NOT EXISTS "${CTEST_DASHBOARD_ROOT}/${CTEST_DIR_NAME}")

# Set the generator and build configuration
if(NOT DEFINED CTEST_CMAKE_GENERATOR)
  set(CTEST_CMAKE_GENERATOR      "Unix Makefiles")
endif(NOT DEFINED CTEST_CMAKE_GENERATOR)
set(CTEST_PROJECT_NAME         "ScaLAPACK")
set(CTEST_BUILD_CONFIGURATION  "Release")

# Extra special variables
set(ENV{DISPLAY}             "")
if(CTEST_CMAKE_GENERATOR MATCHES Makefiles)
  set(ENV{CC}                  "${c_compiler}")
  set(ENV{FC}                  "${f_compiler}")
  set(ENV{CXX}                 "${full_compiler}")
endif(CTEST_CMAKE_GENERATOR MATCHES Makefiles)

###################################################################
# Values for the cmake build
###################################################################

set( CACHE_CONTENTS "
# Specific Fortran Compiler (uncomment and add flags directly after = )
#CMAKE_Fortran_COMPILER:STRING=
# Specific Fortran Compiler Flags (uncomment and add flags directly after = )
#CMAKE_Fortran_FLAGS:STRING=
# Use Reference BLAS and LAPACK by default
USE_OPTIMIZED_LAPACK_BLAS:OPTION=OFF
" )

#----------------------------------------------------------------------------------
# Should not need to edit under this line
#----------------------------------------------------------------------------------

# if you do not want to use the default location for a 
# dashboard then set this variable to the directory
# the dashboard should be in
make_directory("${CTEST_DASHBOARD_ROOT}")
# these are the the name of the source and binary directory on disk. 
# They will be appended to DASHBOARD_ROOT
set(CTEST_SOURCE_DIRECTORY  "${CTEST_DASHBOARD_ROOT}/${CTEST_DIR_NAME}")
set(CTEST_BINARY_DIRECTORY  "${CTEST_SOURCE_DIRECTORY}-${CTEST_BUILD_NAME}")
set(CTEST_NOTES_FILES  "${CTEST_NOTES_FILES}"
  "${CMAKE_CURRENT_LIST_FILE}"
  )

# check for parallel
if(parallel GREATER 1)
  if(NOT CTEST_BUILD_COMMAND)
    set(CTEST_BUILD_COMMAND "make -j${parallel} -i")
  endif(NOT CTEST_BUILD_COMMAND)

  message("Use parallel build")
  message("CTEST_BUILD_COMMAND: ${CTEST_BUILD_COMMAND}")
  message("CTEST_CONFIGURE_COMMAND: ${CTEST_CONFIGURE_COMMAND}")
endif(parallel GREATER 1)

##########################################################################
# wipe the binary dir
message("Remove binary directory...")
ctest_empty_binary_directory("${CTEST_BINARY_DIRECTORY}")

message("CTest Directory: ${CTEST_DASHBOARD_ROOT}")
message("Initial checkout: ${CTEST_CVS_CHECKOUT}")
message("Initial cmake: ${CTEST_CMAKE_COMMAND}")
message("CTest command: ${CTEST_COMMAND}")

# this is the initial cache to use for the binary tree, be careful to escape
# any quotes inside of this string if you use it
file(WRITE "${CTEST_BINARY_DIRECTORY}/CMakeCache.txt" "
SITE:STRING=${hostname}
BUILDNAME:STRING=${BUILDNAME}
DART_ROOT:PATH=
SVNCOMMAND:FILEPATH=${CTEST_UPDATE_COMMAND}
DROP_METHOD:STRING=https
DART_TESTING_TIMEOUT:STRING=${CTEST_TEST_TIMEOUT}
")

message("Start dashboard...")
#ctest_start(Nightly)
ctest_start(Experimental)
message("  Update")
ctest_update(SOURCE "${CTEST_SOURCE_DIRECTORY}" RETURN_VALUE res)
message("  Configure")
ctest_configure(BUILD "${CTEST_BINARY_DIRECTORY}" RETURN_VALUE res)
message("read custom files after configure")
ctest_read_custom_files("${CTEST_BINARY_DIRECTORY}")
message("  Build")
ctest_build(BUILD "${CTEST_BINARY_DIRECTORY}" RETURN_VALUE res)
message("  Test")
ctest_test(BUILD "${CTEST_BINARY_DIRECTORY}" RETURN_VALUE res)
message("  Submit")
ctest_submit(RETURN_VALUE res)
message("  All done")



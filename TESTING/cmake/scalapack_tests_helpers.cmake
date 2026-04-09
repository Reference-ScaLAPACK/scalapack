cmake_minimum_required(VERSION 3.26...4.0)
include_guard()

function(scalapack_add_test test)
    #[===[.md
    # scalapack_add_test

    Internal helper for adding functional tests of entire CMake projects

    ## Synopsis
    ```cmake
    scalapack_add_test(<name>
        [TEST_NAME <test_name>]
        [CMAKE_ARGS <arg1> ...]
        [LABELS <label1> ...]
    )
    ```

    ## Options

    `<name>`
    : Test project to be configured

      This is a subfolder relative to `${CMAKE_CURRENT_SOURCE_DIR}` containing the
      test project with `CMakeLists.txt`

    `TEST_NAME` [Default: `<name>`]
    : Unique name of the ctest test name being created

    `CMAKE_ARGS`
    : Additional CMake args passed to the configure step
    ]===]

    set(ARGS_Options)
    set(ARGS_OneValue
        TEST_NAME
    )
    set(ARGS_MultiValue
        CMAKE_ARGS
        LABELS
    )
    cmake_parse_arguments(PARSE_ARGV 1 ARGS "${ARGS_Options}" "${ARGS_OneValue}" "${ARGS_MultiValue}")
    # Check required/optional arguments
    if(NOT DEFINED ARGS_TEST_NAME)
        set(ARGS_TEST_NAME ${test})
    endif()

    if(SCALAPACK_SOURCE_DIR)
        # If scalapack is built as part of the same project, use the build artifacts
        list(APPEND ARGS_CMAKE_ARGS
            -Dscalapack_ROOT=${SCALAPACK_BINARY_DIR}
        )
    else()
        # Otherwise use `scalapack_DIR` which is created after the `find_package`
        list(APPEND ARGS_CMAKE_ARGS
            -Dscalapack_ROOT=${scalapack_DIR}
        )
    endif()
    # Propagate the compilers used
    list(APPEND ARGS_CMAKE_ARGS
        -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
        -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
    )

    add_test(NAME ${ARGS_TEST_NAME}
        COMMAND ${CMAKE_CTEST_COMMAND} --build-and-test ${CMAKE_CURRENT_SOURCE_DIR}/${test}
        ${CMAKE_CURRENT_BINARY_DIR}/${test}
        # Use the same build environment as the current runner
        --build-generator "${CMAKE_GENERATOR}"
        --build-options
            ${ARGS_CMAKE_ARGS}
        --test-command ${CMAKE_CTEST_COMMAND}
            --test-dir ${CMAKE_CURRENT_BINARY_DIR}/${test}
            --no-tests=ignore
            --output-on-failure
    )
endfunction()

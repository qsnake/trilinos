INCLUDE(PackageCopyFilesToBinaryDir)
INCLUDE(PackageAddTest)

MACRO(PyTrilinos_MAKE_TEST TEST_NAME)

  PACKAGE_COPY_FILES_TO_BINARY_DIR(${TEST_NAME}.py
    SOURCE_FILES ${TEST_NAME}.src
    DEST_FILES   ${TEST_NAME}.py)

  PACKAGE_ADD_TEST(
    ${PYTHON_EXECUTABLE}
    NOEXEPREFIX
    NOEXESUFFIX
    NAME ${TEST_NAME}
    ARGS "${TEST_NAME}.py --testharness"
    STANDARD_PASS_OUTPUT
    ${ARGN}
    )

ENDMACRO(PyTrilinos_MAKE_TEST TEST_NAME)


INCLUDE(PackageMacros)
INCLUDE(AddOptionAndDefine)

PACKAGE(MOOCHO)

ASSERT_DEFINED(${PROJECT_NAME}_ENABLE_Fortran)
ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_QPKWIK
  HAVE_MOOCHO_QPKWIK
  "Enable support for QPKWIK"
  ${${PROJECT_NAME}_ENABLE_Fortran}
  )

ADD_SUBDIRECTORY(src)

ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Thyra)
IF (${PACKAGE_NAME}_ENABLE_Thyra)
  ADD_SUBDIRECTORY(thyra)
ENDIF()

PACKAGE_ADD_EXAMPLE_DIRECTORIES(example)

PACKAGE_ADD_TEST_DIRECTORIES(test)

PACKAGE_ARCH_EXCLUDE_AUTOTOOLS_FILES()

PACKAGE_POSTPROCESS()

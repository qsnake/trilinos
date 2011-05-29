INCLUDE(PackageAddExecutableAndTest)
INCLUDE(AddSubdirectories)


INCLUDE_DIRECTORIES( ${CMAKE_CURRENT_SOURCE_DIR}
                     ${CMAKE_CURRENT_SOURCE_DIR}/.. )


PACKAGE_ADD_EXECUTABLE_AND_TEST(
  CTeuchos_any_UnitTests
  SOURCES
    CTeuchos_any_UnitTests.cpp
    ../CTrilinos_UnitTestHelpers.cpp
    ${PACKAGE_SOURCE_DIR}/../teuchos/test/UnitTest/Teuchos_StandardUnitTestMain.cpp
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  CTeuchos_CommandLineProcessor_UnitTests
  SOURCES
    CTeuchos_CommandLineProcessor_UnitTests.cpp
    ../CTrilinos_UnitTestHelpers.cpp
    ${PACKAGE_SOURCE_DIR}/../teuchos/test/UnitTest/Teuchos_StandardUnitTestMain.cpp
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  CTeuchos_ParameterEntry_UnitTests
  SOURCES
    CTeuchos_ParameterEntry_UnitTests.cpp
    ../CTrilinos_UnitTestHelpers.cpp
    ${PACKAGE_SOURCE_DIR}/../teuchos/test/UnitTest/Teuchos_StandardUnitTestMain.cpp
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  CTeuchos_ParameterList_UnitTests
  SOURCES
    CTeuchos_ParameterList_UnitTests.cpp
    ../CTrilinos_UnitTestHelpers.cpp
    ${PACKAGE_SOURCE_DIR}/../teuchos/test/UnitTest/Teuchos_StandardUnitTestMain.cpp
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )


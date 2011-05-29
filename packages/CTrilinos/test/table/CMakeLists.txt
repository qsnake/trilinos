INCLUDE(PackageAddExecutableAndTest)


#INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})


PACKAGE_ADD_EXECUTABLE_AND_TEST(
  CTrilinos_Table_UnitTests
  SOURCES
    CTrilinos_Table_UnitTests.cpp
    ${PACKAGE_SOURCE_DIR}/../teuchos/test/UnitTest/Teuchos_StandardUnitTestMain.cpp
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  CTrilinos_TableRepos_UnitTests
  SOURCES
    CTrilinos_TableRepos_UnitTests.cpp
    ${PACKAGE_SOURCE_DIR}/../teuchos/test/UnitTest/Teuchos_StandardUnitTestMain.cpp
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  CTrilinos_utils_UnitTests
  SOURCES
    CTrilinos_utils_UnitTests.cpp
    ${PACKAGE_SOURCE_DIR}/../teuchos/test/UnitTest/Teuchos_StandardUnitTestMain.cpp
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )

#PACKAGE_ADD_EXECUTABLE_AND_TEST(
#  CTrilinos_Table_memory
#  SOURCES
#    CTrilinos_Table_memory.c
#  ARGS 12
#  COMM serial
#  NUM_MPI_PROCS 1
#  STANDARD_PASS_OUTPUT
#  )


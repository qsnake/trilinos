
INCLUDE(PackageAddExecutableAndTest)
INCLUDE(PackageCopyFilesToBinaryDir)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../base)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../support-Trilinos)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR}/)

PACKAGE_CONFIGURE_FILE(${PACKAGE_NAME}_config.h)

SET(SOURCES "")

APPEND_SET(SOURCES
fei_unit.cpp
fei_unit_Factory.cpp
  )

PACKAGE_ADD_EXECUTABLE_AND_TEST(
    fei_unit
    SOURCES ${SOURCES}
    ARGS ""
    COMM serial mpi
    NUM_MPI_PROCS 1
   STANDARD_PASS_OUTPUT
  )

PACKAGE_ADD_TEST(
    fei_unit
    NAME fei_unit
    ARGS ""
    COMM serial mpi
    NUM_MPI_PROCS 2
   STANDARD_PASS_OUTPUT
  )

PACKAGE_ADD_TEST(
    fei_unit
    NAME fei_unit
    ARGS ""
    COMM serial mpi
    NUM_MPI_PROCS 3
   STANDARD_PASS_OUTPUT
  )

PACKAGE_ADD_TEST(
    fei_unit
    NAME fei_unit
    ARGS ""
    COMM serial mpi
    NUM_MPI_PROCS 4
   STANDARD_PASS_OUTPUT
  )


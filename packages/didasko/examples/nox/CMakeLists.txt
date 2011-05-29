
INCLUDE(PackageAddExecutableAndTest)

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  Nox_ex1
  SOURCES ex1.cpp
  ARGS
  COMM serial mpi
  NUM_MPI_PROCS 1-4
  )

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  Nox_ex2
  SOURCES ex2.cpp
  ARGS
  COMM serial mpi
  NUM_MPI_PROCS 1-4
  )

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  Nox_ex3
  SOURCES ex3.cpp Heq.C ProblemInterface.C
  ARGS
  COMM serial mpi
  NUM_MPI_PROCS 1-4
  )

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  Nox_ex4
  SOURCES ex4.cpp Heq.C ProblemInterface.C
  ARGS
  COMM serial mpi
  NUM_MPI_PROCS 1-4
  )

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  Nox_ex5
  SOURCES ex5.cpp Heq.C ProblemInterface.C
  ARGS
  COMM serial mpi
  NUM_MPI_PROCS 1-4
  )

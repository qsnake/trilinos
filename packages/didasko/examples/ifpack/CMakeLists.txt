
INCLUDE(PackageAddExecutableAndTest)

PACKAGE_ADD_EXECUTABLE(
  Ifpack_ex1
  SOURCES ex1.cpp
#  ARGS
  COMM serial mpi
#  NUM_MPI_PROCS 1-4
  )

PACKAGE_ADD_EXECUTABLE(
  Ifpack_ex2
  SOURCES ex2.cpp
#  ARGS
  COMM serial mpi
#  NUM_MPI_PROCS 1-4
  )


INCLUDE(PackageAddExecutableAndTest)

PACKAGE_ADD_EXECUTABLE(
  Belos_ex1
  SOURCES ex1.cpp
#  ARGS
  COMM serial mpi
#  NUM_MPI_PROCS 1-4
  )


INCLUDE(PackageAddExecutableAndTest)

PACKAGE_ADD_EXECUTABLE(
  Anasazi_ex1
  SOURCES ex1.cpp
#  ARGS
  COMM serial mpi
#  NUM_MPI_PROCS 1-4
  )

PACKAGE_ADD_EXECUTABLE(
  Anasazi_ex2
  SOURCES ex2.cpp
#  ARGS
  COMM serial mpi
#  NUM_MPI_PROCS 1-4
  )

PACKAGE_ADD_EXECUTABLE(
  Anasazi_ex3
  SOURCES ex3.cpp
#  ARGS
  COMM serial mpi
#  NUM_MPI_PROCS 1-4
  )

PACKAGE_ADD_EXECUTABLE(
  Anasazi_ex4
  SOURCES ex4.cpp
#  ARGS
  COMM serial mpi
#  NUM_MPI_PROCS 1-4
  )

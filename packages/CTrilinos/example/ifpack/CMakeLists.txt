
INCLUDE(PackageAddExecutableAndTest)
INCLUDE(AddSubdirectories)

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  Ifpack_example_Amesos
  SOURCES Ifpack_example_Amesos.c
  COMM serial mpi
  NUM_MPI_PROCS 1
  )

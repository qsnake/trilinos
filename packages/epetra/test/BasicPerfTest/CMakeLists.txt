INCLUDE(PackageAddExecutableAndTest)

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  BasicPerfTest_test 
  SOURCES cxx_main.cpp
  ARGS "16 12 1 1 25 -v"
  COMM mpi serial
  NUM_MPI_PROCS 1
  FAIL_REGULAR_EXPRESSION "tests FAILED"
  )

INSTALL(TARGETS Epetra_BasicPerfTest_test DESTINATION bin)


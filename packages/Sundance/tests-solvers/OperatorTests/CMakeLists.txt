# CMake tests specification for Sundance/tests-solvers/OperatorTests


INCLUDE(AddTestBatch)



SET(MPITests 
        LinearCombinationTest 
        CompoundOpTest
        MultiVectorOpTest
        BlockOpTest)

SET(SerialTests )
#  SVDTest)


ADD_TEST_BATCH(
        MPITests
        COMM serial mpi
)

ADD_TEST_BATCH(
        SerialTests
        COMM serial mpi
        NUM_MPI_PROCS 1
)

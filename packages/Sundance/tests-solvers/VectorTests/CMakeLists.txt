# CMake tests specification for Sundance/tests-solvers/VectorTests


INCLUDE(PackageAddExecutableAndTest)

PACKAGE_ADD_EXECUTABLE_AND_TEST(
        VectorTest
        SOURCES VectorTest.cpp
        COMM serial mpi
)

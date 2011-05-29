# CMake tests specification 


INCLUDE(PackageAddExecutableAndTest)

PACKAGE_ADD_EXECUTABLE_AND_TEST(
        FunctorTest
        SOURCES FunctorTest.cpp
        COMM serial
)

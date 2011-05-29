# CMake tests specification 


INCLUDE(PackageAddExecutableAndTest)

PACKAGE_ADD_EXECUTABLE_AND_TEST(
        mapPerfTest
        SOURCES mapPerfTest.cpp
        COMM serial
)

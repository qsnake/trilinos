# CMake tests specification 


INCLUDE(PackageAddExecutableAndTest)

PACKAGE_ADD_EXECUTABLE_AND_TEST(
        CombTest
        SOURCES CombTest.cpp
        COMM serial
)


INCLUDE("${CTEST_SCRIPT_DIRECTORY}/TrilinosCTestDriverCore.zan.gcc.cmake")
INCLUDE("${CTEST_SCRIPT_DIRECTORY}/../../../TrilinosVersion.cmake")

#
# Set the options specific to this build case
#

SET(COMM_TYPE SERIAL)
SET(BUILD_TYPE OPT)
SET(BUILD_DIR_NAME "CUDA_OPT_RELEASE_GCC_${Trilinos_VERSION}")
SET(Trilinos_TRACK ${Trilinos_TESTING_TRACK})
SET(Trilinos_BRANCH ${Trilinos_REPOSITORY_BRANCH})
#SET(CTEST_TEST_TYPE Experimental)

SET(Trilinos_PACKAGES Tpetra Kokkos)

SET( EXTRA_CONFIGURE_OPTIONS
  "-DDART_TESTING_TIMEOUT:STRING=120"
  "-DTrilinos_ENABLE_EXPLICIT_INSTANTIATION:BOOL=ON"
  "-DTrilinos_ENABLE_DEVELOPMENT_MODE:BOOL=OFF"
  "-DTPL_ENABLE_Pthread:BOOL=ON"
  "-DTPL_ENABLE_CUDA:BOOL=ON"
  "-DTPL_ENABLE_Thrust:BOOL=ON"
  "-DTPL_ENABLE_TBB:BOOL=ON"
  "-DTBB_LIBRARY_DIRS=/usr/local/tbb30_018oss/lib/intel64/cc4.1.0_libc2.4_kernel2.6.16.21"
  "-DTBB_INCLUDE_DIRS=/usr/local/tbb30_018oss/include"
  "-DThrust_INCLUDE_DIRS=/usr/local"
  "-DTrilinos_ENABLE_PERFORMANCE_TESTS:BOOL=ON"
  )

#
# Set the rest of the system-specific options and run the dashboard build/test
#

TRILINOS_SYSTEM_SPECIFIC_CTEST_DRIVER()
